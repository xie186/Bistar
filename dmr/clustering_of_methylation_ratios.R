#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Check that packages are available, otherwise raise an error

required_packages = c("argparse", "logging", "bsseq", "ggplot2")

# Look for missing packages
missing_packages <- required_packages[
    !(required_packages %in% installed.packages())]

# Throw an error if any package is missing.
if(length(missing_packages) > 0)
    stop(paste("Missing R package(s):", paste(missing_packages, collapse=" ")))

library(logging)
library(ggplot2)

# -----------------------------------------------------------------------------
# Get command line arguments

# Create parser
parser <- argparse::ArgumentParser()

# Add required arguments
parser$add_argument(
    "--bs-rds", required=T, metavar="<path>", help=paste(
    "Path to the R RDS file containing the 'bs' BSseq object which provides",
    "methylation counts and covariates.")
)
parser$add_argument(
    "--outdir", required=T, metavar="<path>", help="Output directory."
)
parser$add_argument(
    "--test-covariate", required=T, metavar="<name>", help=paste(
    "The name of the covariate to test for association of methylation levels.")
)

# Optional arguments
parser$add_argument(
    "--min-cpgs", required=F, metavar="<int>", default=100, help=paste(
    "Minimum number of CpGs required to run the PCA considering that only",
    "CpGs measured in all samples are used, by default 100.")
)

# Parse
args <- parser$parse_args()

# -----------------------------------------------------------------------------
# Setup logger and log parameters

basicConfig()
loginfo("Parameters:")
loginfo(str(args, nchar.max=200))

# -----------------------------------------------------------------------------
# Load methylation counts
bs <- readRDS(args$bs_rds)

# Log number of CpGs and samples
loginfo("Number of CpGs: %i.", nrow(bs))
loginfo("Number of samples: %i.", ncol(bs))

# -----------------------------------------------------------------------------
# Compute methylation ratios of CpGs measured in all samples

loginfo("Keeping only CpGs measured in all samples before clustering.")

# Get coverage matrix
coverage <- bsseq::getCoverage(bs, type="Cov")

# Identify CpGs that are not measured in all samples
rows_to_remove <- which(DelayedMatrixStats::rowSums2(coverage == 0) > 0)

# Remove those CpGs
loginfo("Removing %i CpGs...", length(rows_to_remove))
if (length(rows_to_remove) > 0) {
    bs <- bs[-rows_to_remove, ]
}

# Check that there are enough CpGs
nb_cpgs <- nrow(bs)
if (nb_cpgs < args$min_cpgs) {
    logerror("Not enough CpGs left: %i, required: %i.", nb_cpgs, args$min_cpgs)
    stop()
} else {
    loginfo("Number of CpGs left for the clustering: %i.", nb_cpgs)
}

# Free memory
rm(coverage, rows_to_remove)

# Compute methylation ratios
ratios <- bsseq::getCoverage(bs, type="M") / bsseq::getCoverage(bs, type="Cov")

# Tranpose to have samples in rows
ratios <- t(ratios)

# -----------------------------------------------------------------------------
# Hierarchical clustering: compute and plot dendrogram

# Number of different colors needed to plot the dendrogram
unique_values <- unique(bsseq::pData(bs)[[args$test_covariate]])
nb_colors <- length(unique_values)

# Check that the test covariate takes at least 2 different values
if (nb_colors < 2)
    stop("All the samples have the same test covariate value.")

# If the test covariate has 2 possible values, use "blue" and "red"
if (nb_colors == 2) {
    # Create vector mapping <test covariate value> -> <color>
    color_of_test_covariate <- c("blue", "red")
    names(color_of_test_covariate) <- unique_values

# If more than 2 colors are needed, assign colors using ColorBrewer package
} else {
    # Use 'Set1' colormap for categorical values and 'Reds' if continuous
    cmap_name <- if (is.numeric(unique_values)) "Set1" else "Reds"
    color_of_test_covariate <- RColorBrewer::brewer.pal(n=nb_colors,
                                                        name=cmap_name)
}

# Compute dissimilarity between samples
loginfo("Computing dissimilarity matrix using euclidean distance...")
distances <- dist(ratios, method="euclidean")

# Hierarchical clustering
loginfo("Clustering using ward.D2 method...")
hcluster <- hclust(distances, method="ward.D2")

# Dendrogram
loginfo("Computing dendrogram...")
dendrogram <- as.dendrogram(hcluster)

## Assign a color to each value of the test covariate

# Function to pass to dendrapply() to set color of leaves according to the
# test covariate value
add_color_to_dendrogram_leaves <- function(node) {
    if (is.leaf(node)) {
        a <- attributes(node)
        test_covariate <- bsseq::pData(bs)[a$label, args$test_covariate]
        color <- color_of_test_covariate[test_covariate]
        attr(node, "nodePar") <- c(a$nodePar, list(lab.col=color))
    }
    return(node)
}

# Colorize dendrogram
dendrogram <- dendrapply(dendrogram, add_color_to_dendrogram_leaves)

# Plot and save dendrogram as an SVG
dendro_svg <- file.path(args$outdir, "dendrogram_of_methylation_ratios.svg")
loginfo("Creating dendrogram plot...")
svg(file=dendro_svg)
dendro_plot <- plot(
    dendrogram,
    main="Clustering of methylation ratios of CpGs measured in all samples")
loginfo("Saving dendrogram plot as %s", dendro_svg)
print(dendro_svg)
dev.off()

# -----------------------------------------------------------------------------
# Compute and plot PCA

loginfo("Computing PCA...")
pca <- stats::prcomp(x=ratios, center=T, scale=F)

# Extract the first 2 components as a dataframe and add a "color" column
# to adapt color to test covariate
df <- as.data.frame(pca$x)[, 1:2]
df$color <- bsseq::pData(bs)[, args$test_covariate]

# Compute percentage of variance explained by each component
percentages <-pca$sdev / sum(pca$sdev) * 100

# Plot PCA and save as an SVG
pca_svg <- file.path(args$outdir, "pca_of_methylation_ratios.svg")
title <- "PCA of methylation ratios of CpGs measured in all samples"
loginfo("Creating PCA plot...")
svg(file=pca_svg)
pca_plot <- ggplot(df, aes(x=PC1, y=PC2, color=color, label=rownames(df))) +
    geom_point() + geom_text(vjust="inward", hjust="inward", size=2) +
    theme_bw() + labs(title=title,
                      x=sprintf("PC1 %.2f %%", percentages[1]),
                      y=sprintf("PC2 %.2f %%", percentages[2]),
                      color=args$test_covariate)
loginfo("Saving PCA plot as %s", pca_svg)
print(pca_plot)
dev.off()
