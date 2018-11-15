#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Check that packages are available, otherwise raise an error

required_packages = c("argparse", "logging", "goldmine")

# Look for missing packages
missing_packages <- required_packages[
    !(required_packages %in% installed.packages())]

# Throw an error if any package is missing.
if(length(missing_packages) > 0)
    stop(paste("Missing R package(s):", paste(missing_packages, collapse=" ")))

library(logging)
library(goldmine)

# -----------------------------------------------------------------------------
# Get command line arguments

# Create parser
parser <- argparse::ArgumentParser()

# Required arguments
parser$add_argument(
    "--dmrs-rds", required=T, metavar="<path>",
    help="Path to the dmrseq results in R RDS format."
)
parser$add_argument(
    "--outdir", required=T, metavar="<path>", help="Output directory."
)
parser$add_argument(
    "--ref-build", required=T, metavar="<name>",
    help="Reference genome build, e.g. hg19."
)

# Optional arguments
parser$add_argument(
    "--goldmine-dir", metavar="<path>", help=paste(
    "Directory to use a cache directory for Goldmine annotation package.",
    "It avoids re-downloading annotation databases from UCSC if you",
    "already have them.")
)
parser$add_argument(
    "--goldmine-sync", action="store_true", help=paste(
    "Check if newer versions of UCSC annotation tables are available and if",
    "so download them in <--goldmine-dir>.")
)

args <- parser$parse_args()

# -----------------------------------------------------------------------------
# Setup logger and log parameters

basicConfig()
loginfo("Parameters:")
loginfo(str(args, nchar.max=200))

# -----------------------------------------------------------------------------
# Checks

# Check existence of input file
if(! file.exists(args$dmrs_rds)) {
    logerror("Input RDS file does not exist: %s", args$dmrs_rds)
    stop()
}

# Check existence of output directory
if(! dir.exists(args$outdir)) {
    logerror("Output directory does not exist: %s", args$outdir)
    stop()
}

# -----------------------------------------------------------------------------
# Load DMRs
dmrs <- readRDS(args$dmrs_rds)

# -----------------------------------------------------------------------------
# Annotate regions

# Get genes from UCSC
ucsc_genes <- goldmine::getGenes(geneset="ucsc",
                                 genome=args$ref_build,
                                 cachedir=args$goldmine_dir,
                                 sync=args$goldmine_sync)

# Annotate with genes
annotated_dmrs <- goldmine::goldmine(query=dmrs,
                                     genome=args$ref_build,
                                     cachedir=args$goldmine_dir,
                                     genes=ucsc_genes,
                                     contextonly=T,
                                     sync=args$goldmine_sync)

# Keep only columns of interest
columns_to_keep <- c(
    "chr", "start", "end", "width", "strand", "L", "area", "beta", "stat",
    "pval", "qval", "call", "call_genes", "overlapped_genes", "nearest_genes",
    "distance_to_nearest_gene")
annotated_dmrs <- annotated_dmrs[, columns_to_keep, with=F]

# Save DMRs as a BED file
annotated_dmrs_bed <- file.path(args$outdir, "annotated_dmrs.bed")
write.table(annotated_dmrs, file=annotated_dmrs_bed, sep="\t",
            col.names=T, row.names=F, quote=F)
