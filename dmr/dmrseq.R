#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Check that packages are available, otherwise raise a clear error

required_packages = c("argparse", "logging", "dmrseq")

# Look for missing packages
missing_packages <- required_packages[
    !(required_packages %in% installed.packages())]

# Throw an error if any package is missing.
if(length(missing_packages) > 0)
    stop(paste("Missing R package(s):", paste(missing_packages, collapse=" ")))

library(logging)

# -----------------------------------------------------------------------------

# Default values
DEFAULT_OF <- list(
    cutoff=0.10,
    min_cpgs=5,
    max_gap=100,
    smooth=TRUE,
    bp_span=1000,
    min_in_span=30,
    max_gap_smooth=2500,
    max_perms=20,
    threads=1
)

# Get command line arguments
parser <- argparse::ArgumentParser()

# Required arguments
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
    "dmrseq 'testCovariate' argument: the name of the covariate to test for",
    "association of methylation levels. Two-group comparison is carried out",
    "if the variable has 2 levels, but dmrseq can also work with continuous",
    "or categorical variables that have more than 2 levels.")
)

# Optional arguments
parser$add_argument(
    "--threads", type="integer", default=DEFAULT_OF[["threads"]],
    metavar="<int>", help=paste(
    "Number of threads, by default", DEFAULT_OF[["threads"]])
)
parser$add_argument(
    "--max-gap", type="integer", default=DEFAULT_OF[["max_gap"]],
    metavar="<int>", help=paste(
    "dmrseq 'maxGap' argument: max basepair distance between neighboring CpGs",
    "to be included in the same DMR, by default", DEFAULT_OF[["max_gap"]])
)
parser$add_argument(
    "--min-cpgs", type="integer", default=DEFAULT_OF[["min_cpgs"]],
    metavar="<int>", help=paste(
    "dmrseq 'minNumRegion' argument: minimum number of CpGs in a DMR, by",
    "default", DEFAULT_OF[["min_cpgs"]])
)
parser$add_argument(
    "--cutoff", type="double", default=DEFAULT_OF[["cutoff"]],
    metavar="<float>", help=paste(
    "dmrseq 'cutoff' argument, by default", DEFAULT_OF[["cutoff"]])
)
parser$add_argument(
    "--no-smoothing", action="store_false", dest="smooth",
    help="Set dmrseq 'smooth' argument to False."
)
parser$add_argument(
    "--bp-span", type="integer", default=DEFAULT_OF[["bp_span"]],
    metavar="<int>", help=paste(
    "dmrseq 'bpSpan' argument: smoothing window in basepairs, by default",
    DEFAULT_OF[["bp_span"]])
)
parser$add_argument(
    "--min-in-span", type="integer", default=DEFAULT_OF[["min_in_span"]],
    metavar="<int>", help=paste(
    "dmrseq 'minInSpan' argument: minimum numober of CpGs in a smoothing span",
    "window, by default", DEFAULT_OF[["min_in_span"]])
)
parser$add_argument(
    "--max-gap-smooth", type="integer", default=DEFAULT_OF[["max_gap_smooth"]],
    metavar="<int>", help=paste(
    "dmrseq 'maxGapSmooth' argument: maximum basepair distance between",
    "neighboring CpGs when smoothing, by default",
    DEFAULT_OF[["max_gap_smooth"]])
)
parser$add_argument(
    "--max-perms", type="integer", default=DEFAULT_OF[["max_perms"]],
    metavar="<int>", help=paste(
    "dmrseq 'maxPerms' argument: maximum number of permutations to generate",
    "the global null distribution, by default", DEFAULT_OF[["max_perms"]])
)
parser$add_argument(
    "--adjust-covariates", nargs="+", metavar="<name>", help=paste(
    "dmrseq 'adjusteCovariate' argument: the names of covariates to be",
    "adjusted for when testing for the association of methylation value with",
    "the test covariate.")
)
parser$add_argument(
    "--match-covariate", metavar="<name>", help=paste(
    "dmrseq 'matchCovariate' argument: the name a of the covariate that will",
    "be blocked when permutating to test for the association of methylation ",
    "with the test covariate. Only possible for a two-group comparison and",
    "for one covariate.")
)

args <- parser$parse_args()

# -----------------------------------------------------------------------------
# Setup logger and log parameters

basicConfig()
loginfo("Parameters:")
loginfo(str(args, nchar.max=200))

# -----------------------------------------------------------------------------
# Load BSseq data 
bs <- readRDS(args$bs_rds)

# Print number of CpGs and samples
loginfo("Number of CpGs: %i.", nrow(bs))
loginfo("Number of samples: %i.", ncol(bs))

# Names of the covariates to be used in the analysis
covariates_to_use <- c(args$test_covariate, args$adjust_covariates,
                       args$match_covariate)

# Check that the covariates to be used are in the BSseq object
for(cov_name in covariates_to_use)
    if(! (cov_name %in% names(bsseq::pData(bs)))) {
        logerror("Covariate '%s' not found in BSseq object: 'bs'.", cov_name)
        stop()
    }

# -----------------------------------------------------------------------------

# Detect candidate regions
start_time <- Sys.time()
dmrs <- dmrseq::dmrseq(
    bs=bs,
    testCovariate=args$test_covariate,
    adjustCovariate=args$adjust_covariates,
    matchCovariate=args$match_covariate,
    cutoff=args$cutoff,
    minNumRegion=args$min_cpgs,
    smooth=args$smooth,
    bpSpan=args$bp_span,
    minInSpan=args$min_in_span,
    maxGapSmooth=args$max_gap_smooth,
    maxGap=args$max_gap,
    maxPerms=args$max_perms,
    BPPARAM=BiocParallel::MulticoreParam(workers=args$threads),
    verbose=TRUE)

# Log duraction
loginfo("Execution time: %s", Sys.time() - start_time)

# Remove any row that has a NA value
dmrs <- na.omit(dmrs)

# Write results as a RDS file
dmrs_rds <- file.path(args$outdir, "dmrs.rds")
saveRDS(dmrs, file=dmrs_rds)
