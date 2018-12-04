#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Check that packages are available, otherwise raise an error

required_packages = c("argparse", "logging", "bsseq")

# Look for missing packages
missing_packages <- required_packages[
    !(required_packages %in% installed.packages())]

# Throw an error if any package is missing.
if(length(missing_packages) > 0)
    stop(paste("Missing R package(s):", paste(missing_packages, collapse=" ")))

library(logging)

# -----------------------------------------------------------------------------
# Get command line arguments

# Default values
DEFAULT_OF <- list(
    threads=1
)

# Create parser
parser <- argparse::ArgumentParser()

# Add required arguments
parser$add_argument(
    "--outdir", required=T, metavar="<path>", help="Output directory."
)
parser$add_argument(
    "--samples", required=T, nargs="+", metavar="<ID>", help="The sample IDs."
)
parser$add_argument(
    "--cpg-reports", required=T, nargs="+",  metavar="<path>", help=paste(
    "Filepaths of Bismark CpG reports in the same order as the related sample",
    "IDs provided with --samples.")
)
parser$add_argument(
    "--covariates-tsv", required=T, metavar="<path>", help=paste(
    "Path to the covariate file, providing the test covariate and optionaly",
    "the covariates to be accounted for. The file has to be tab-separated",
    "with a header with the name(s) of the covariate(s). The first column",
    "provides the sample IDs.")
)
parser$add_argument(
    "--test-covariate", required=T, metavar="<name>", help=
    "The name of the covariate to test for association of methylation levels."
)
parser$add_argument(
    "--min-samples", required=T, type="integer", metavar="<int>",
    help="Minimum number of measured samples, otherwise the CpG is removed."
)

# Optional arguments
parser$add_argument(
    "--do-not-merge-complementary-cpgs", action="store_false",
    dest="merge_complementary_cpgs", help=paste(
    "By default a CpG and its complementary CpG on the other strand are",
    "merged as one CpG because we consider methylation as symmetric.",
    "Set this flag to not merge and keep 2 separate CpGs.")
)
parser$add_argument(
    "--threads", type="integer", default=DEFAULT_OF[["threads"]],
    metavar="<int>", help=paste(
    "Number of threads, by default", DEFAULT_OF[["threads"]])
)

args <- parser$parse_args()

# -----------------------------------------------------------------------------
# Setup logger and log parameters

basicConfig()
loginfo("Parameters:")
loginfo(str(args, nchar.max=200))

# -----------------------------------------------------------------------------
# Load covariates and make checks

# Check existence of file
if(! file.exists(args$covariates_tsv)) {
    logerror("Covariate file does not exist: %s", args$covariates_tsv)
    stop()
}

# Check existence of output directory
if(! dir.exists(args$outdir)) {
    logerror("Output directory does not exist: %s", args$outdir)
    stop()
}

# Load table
covariates <- read.table(args$covariates_tsv, sep="\t", header=T)

# Assuming first column provides the sample IDs, rename it to "sampleID"
names(covariates)[1] <- "sampleID"

# Multiple checks on the dataframe
if(ncol(covariates) < 2) {
    logerror("The covariate file should have at least 2 columns: %s",
             args$covariates_tsv)
    stop()
}
    
for(sample_id in args$samples) {
    if(! (sample_id %in% covariates$sampleID)) {
        logerror("Sample '%s' not found in the first column of %s", sample_id,
                 args$covariates_tsv)
        stop()
    }
}

# Check that the test covariate is in the dataframe
if(! (args$test_covariate %in% names(covariates))) {
    logerror("Test covariate '%s' not found in the covariate file.",
             args$test_covariate)
    stop()
}

# Check that the test covariate takes at least 2 different values
unique_values <- unique(covariates[[args$test_covariate]])
if (length(unique_values) < 2) {
    logerror("All the samples have the same test covariate value.")
    stop()
}

# Keep only samples provided with --samples argument
covariates <- covariates[which(covariates$sampleID %in% args$samples), ]

# Reorder the samples to match the order in the BSseq object
covariates <- covariates[match(args$samples, covariates$sampleID), ]

# Make row.names be sample IDs (required by read.bismark colData argument)
row.names(covariates) <- covariates$sampleID

# -----------------------------------------------------------------------------

# Read methylation calls
bs <- bsseq::read.bismark(files=args$cpg_reports,
                          colData=covariates,
                          rmZeroCov=TRUE,
                          strandCollapse=args$merge_complementary_cpgs,
                          nThread=args$threads,
                          verbose=TRUE)

# Log
loginfo("Number of CpGs: %i.", nrow(bs))
loginfo("Number of samples: %i.", ncol(bs))

# -----------------------------------------------------------------------------
# Filter loci without enough measured samples

# Get total read counts
coverage <- bsseq::getCoverage(bs, type="Cov")

# Indexes of CpGs without enough samples
not_enough_samples <- which(
    DelayedMatrixStats::rowSums2(coverage > 0) < args$min_samples)

# Print filtering consequences
loginfo("Removing %i CpGs that are not measured in at least %i samples",
        length(not_enough_samples), args$min_samples)

# Remove low coverage CpGs
if (length(not_enough_samples) > 0) {
    bs <- bs[-not_enough_samples, ]
}

# Free memory
rm(coverage, not_enough_samples)

# Print filtering consequences
nb_cpgs_left <- nrow(bs)
loginfo("%i CpG left after filtering.", nb_cpgs_left)

# If nothing left, stop execution
if(nb_cpgs_left == 0) {
    logerror("0 CpG left after filtering, '--min-samples' might be too high.")
    stop()
}

# For a 2-group test covariate or a categorical test covariate with more than
# 2 levels, remove CpGs that don't have at least one sample per condition.
is_2_group <- length(unique_values) == 2
is_categorical <- !is.numeric(unique_values)
if (is_2_group || is_categorical) {
    for (condition in unique_values) {
        
        # Log
        loginfo(paste("Checking that all CpGs have at least one measured for",
                      "condition '%s'."), condition)
        loginfo("Number of CpGs before filtering: %i.", nrow(bs))

        # Get samples of that condition
        condition_indexes <- which(
            bsseq::pData(bs)[[args$test_covariate]] == condition)
        condition_samples <- bsseq::sampleNames(bs)[condition_indexes]

        # Log samples of that condition
        loginfo("Samples of this condition: %s", condition_samples)

        # Get coverage of samples of this condition
        condition_coverage <- bsseq::getCoverage(
            bs[, condition_indexes], type="Cov")

        # Remove CpGs with 0 coverage for all samples of that condition
        zero_cov_cpgs <- which(
            DelayedMatrixStats::rowSums2(condition_coverage > 0) == 0)

        # If CpGs have been removed print filtering consequences
        nb_zero_cov_cpgs <- length(zero_cov_cpgs)
        if (nb_zero_cov_cpgs > 0) {
            loginfo(paste("Removing %i CpGs that have 0 coverage for all",
                          "samples of condition '%s'"),
                    nb_zero_cov_cpgs, condition)
            bs <- bs[-zero_cov_cpgs, ]
        } else {
            loginfo(paste("All CpGs have at least one measured sample for",
                          "that condition."))
        }
        
        # Log
        loginfo("Number of CpGs after filtering: %i.", nrow(bs))

        # Clean
        rm(condition_indexes, condition_samples, condition_coverage,
           zero_cov_cpgs, nb_zero_cov_cpgs)
    }
}

# Write BSseq object as a RData file
bs_rds <- file.path(args$outdir, "bs.rds")
loginfo("Save methylation counts as RDS file: %s ...", bs_rds)
saveRDS(bs, file=bs_rds)
