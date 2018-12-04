#!/usr/bin/env python

# Standard import
import os
import argparse
import re
import snakemake
import itertools
import pandas
import yaml

DOC = """
Script that generates the YAML configuration file required to run the
differential methylation analysis pipeline.

The Bismark CpG report files to process along with the sample IDs are
identified using a Snakemake regex. Only the '{sample}' wildcard is required in
the regex. You can add arbitrary wildcards to match variable fields in the
filepath, they will be ignored.

Examples of regex (--cpg-report-regex argument):
    /path/to/data/{sample}/{sample}.CpG_report.txt.gz
    /path/to/data/{sample}/{sample}.{ignore}.CpG_report.txt.gz

All parameters and paths to files to process are stored in the output YAML.
"""

def get_goldmine_default_cache_dir():
    bistar_dir = os.environ.get("BISTAR_DIR", None)
    if bistar_dir is not None:
        return os.path.join(bistar_dir, "Goldmine")
    else:
        return None

# Default values
DEFAULT_OF = dict(
    bsseq_threads=1,
    cutoff=0.10,
    dmrseq_threads=1,
    min_cpgs=5,
    min_samples=2,
    max_gap=100,
    smooth=True,
    bp_span=1000,
    min_in_span=30,
    max_gap_smooth=2500,
    max_perms=10,
    goldmine_dir=get_goldmine_default_cache_dir()
)


def is_valid_regex(snakemake_regex):
    """
    Check that the provided pattern contains the required wildcard.
    """
    if "{sample}" not in snakemake_regex:
        raise argparse.ArgumentTypeError(
            "The regex should contain the wildcard '{sample}'.")
    return snakemake_regex


def is_file(filepath):
    """ Check file's existence - argparse 'type' argument. """
    if not os.path.isfile(filepath):
        raise argparse.ArgumentTypeError("File does not exist: %s" % filepath)
    return filepath


def get_cmd_line_args():
    """
    Create a command line argument parser and return a dict mapping
    <argument name> -> <argument value>.
    """
    parser = argparse.ArgumentParser(
        description=DOC, formatter_class=argparse.RawTextHelpFormatter)

    # Required arguments
    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "--cpg-report-regex", type=is_valid_regex, required=True,
        metavar="<regex>", help=
        "Snakemake regex used to infer the Bismark CpG report files to process"
        " and the related sample IDs. The '{sample}' wildcard is required in "
        "the regex. Arbitrary wildcards can be used to match for variable "
        "fields in the filepath. Examples: "
        "'/path/to/data/{sample}/{sample}.CpG_report.txt.gz' or "
        "'/path/to/data/{sample}_{ignore}.CpG_report.txt.gz'"
    )
    required.add_argument(
        "--outdir", required=True, metavar="<dir>", help="Output directory."
    )
    required.add_argument(
        "--ref-build", metavar="<name>", help=
        "Reference genome build e.g. hg19 or mm10. The Goldmine annotation "
        "will only work if the reference build name is a UCSC genome name."
    )
    required.add_argument(
        "--covariates-tsv", required=True, type=is_file, metavar="<path>",
        help=
        "Path to the covariate file, providing the test covariate and "
        "optionally the covariates to be accounted for. The file has to be "
        "tab-separated with a header with the name(s) of the covariate(s). "
        "The first column provides the sample IDs."
    )
    required.add_argument(
        "--test-covariate", required=True, metavar="<name>",
        help=
        "dmrseq 'testCovariate' argument: the name of the covariate to test "
        "for association of methylation levels. Two-group comparison is car"
        "ried out if the variable has 2 levels, but dmrseq can also work with "
        "continuous or categorical variables that have more than 2 levels.")

    # Optional arguments
    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument(
        "--bsseq-threads", type=int, default=DEFAULT_OF["bsseq_threads"],
        metavar="<int>", help=
        "Number of threads to use when loading methylation counts with BSseq, "
        "by default %i. Increasing the number of threads increases RAM "
        "consumption." % DEFAULT_OF["bsseq_threads"]
    )
    optional.add_argument(
        "--dmrseq-threads", type=int, default=DEFAULT_OF["dmrseq_threads"],
        metavar="<int>", help=
        "dmrseq 'threads' argument, by default %i. Increasing the number of "
        "threads increases RAM consumption." % DEFAULT_OF["dmrseq_threads"]
    )
    optional.add_argument(
        "--do-not-merge-complementary-cpgs", action="store_false",
        dest="merge_complementary_cpgs", help=
        "By default a CpG and its complementary CpG on the other strand are "
        "merged as one CpG because we consider methylation as symmetric. "
        "Set this flag to keep 2 separate CpGs."
    )
    optional.add_argument(
        "--max-gap", type=int, default=DEFAULT_OF["max_gap"], metavar="<int>",
        help=
        "dmrseq 'maxGap' argument: max basepair distance between neighboring "
        "CpGs to be included in the same DMR, by default %i."
        % DEFAULT_OF["max_gap"]
    )
    optional.add_argument(
        "--min-cpgs", type=int, default=DEFAULT_OF["min_cpgs"],
        metavar="<int>", help=
        "dmrseq 'minNumRegion' argument: minimum number of CpGs in a DMR, by "
        "default %i." % DEFAULT_OF["min_cpgs"]
    )
    optional.add_argument(
        "--min-samples", type=int, metavar="<int>",
        default=DEFAULT_OF["min_samples"], help=
        "Minimum number of measured samples, otherwise the CpG is removed, "
        "by default %i." % DEFAULT_OF["min_samples"]
    )
    optional.add_argument(
        "--cutoff", type=float, default=DEFAULT_OF["cutoff"],
        metavar="<float>", help=
        "dmrseq 'cutoff' argument, by default %f." % DEFAULT_OF["cutoff"]
    )
    optional.add_argument(
        "--no-smoothing", action="store_false", dest="smooth",
        help="Set dmrseq 'smooth' argument to False."
    )
    optional.add_argument(
        "--bp-span", type=int, default=DEFAULT_OF["bp_span"], metavar="<int>",
        help="dmrseq 'bpSpan' argument: smoothing window in basepairs, by "
             "default %i." % DEFAULT_OF["bp_span"]
    )
    optional.add_argument(
        "--min-in-span", type=int, default=DEFAULT_OF["min_in_span"],
        metavar="<int>", help=
        "dmrseq 'minInSpan' argument: minimum numober of CpGs in a smoothing "
        "span window, by default %i." % DEFAULT_OF["min_in_span"]
    )
    optional.add_argument(
        "--max-gap-smooth", type=int, default=DEFAULT_OF["max_gap_smooth"],
        metavar="<int>", help=
        "dmrseq 'maxGapSmooth' argument: max basepair distance between "
        "neighboring CpGs when smoothing, by default %i."
        % DEFAULT_OF["max_gap_smooth"]
    )
    optional.add_argument(
        "--max-perms", type=int, default=DEFAULT_OF["max_perms"],
        metavar="<int>", help=
        "dmrseq 'maxPerms' argument: maximum number of permutations to "
        "generate the global null distribution, by default %i."
        % DEFAULT_OF["max_perms"]
    )
    optional.add_argument(
        "--adjust-covariates", nargs="+", metavar="<name>", help=
        "dmrseq 'adjusteCovariate' argument: the names of the covariates to "
        "be adjusted for when testing for the association of methylation "
        "value with the test covariate."
    )
    optional.add_argument(
        "--match-covariate", metavar="<name>", help=
        "dmrseq 'matchCovariate' argument: a covariate that will be blocked "
        "when permutating to test for the association of methylation with the "
        "test covariate. Only permutations with balanced composition will be "
        "used. Only possible for a two-group comparison and for one covariate."
    )
    optional.add_argument(
        "--goldmine-dir", metavar="<dir>", default=DEFAULT_OF["goldmine_dir"],
        help=
        "Directory to use a cache directory for Goldmine annotation package. "
        "It avoids re-downloading annotation databases from UCSC if you "
        "already have them. By default %s ." % DEFAULT_OF["goldmine_dir"]
    )
    optional.add_argument(
        "--goldmine-sync",  action="store_true", help=
        "Check if newer versions of UCSC annotation tables are available and "
        "if so download them in <--goldmine-dir>. It requires a web "
        "connection which is not always the case if you run on a cluster."
    )

    # Parse the command line
    args = parser.parse_args()

    # Check that the --min-samples argument value is > 0
    if not args.min_samples > 0:
        raise ValueError("the --min-samples argument should be > 0.")

    # Convert the argparse object to a dict mapping <arg name> -> <value>
    kwargs = vars(args)

    return kwargs


# Function copied from snakemake.utils.listfiles but modified to make it work
# with symbolic links:
# 'os.walk(dirname)' replaced by 'os.walk(dirname, followlinks=True)'
# Since the 2 parameters 'restriction' and 'omit_value' are not used the code
# was also reduced.
def listfiles(pattern):
    """
    Yield a tuple of existing filepaths for the given pattern.
    Wildcard values are yielded as the second tuple item.

    Args:
        pattern (str): a filepattern. Wildcards are specified in snakemake
                       syntax, e.g. "{id}.txt"
    Yields:
        tuple: The next file matching the pattern, and the corresponding
               wildcards object
    """
    pattern = os.path.normpath(pattern)
    first_wildcard = re.search("{[^{]", pattern)
    if first_wildcard:
        dirname = os.path.dirname(pattern[:first_wildcard.start()])
        if not dirname:
            dirname = "."
    else:
        dirname = os.path.dirname(pattern)
    pattern = re.compile(snakemake.io.regex(pattern))
    for dirpath, dirnames, filenames in os.walk(dirname, followlinks=True):
        for f in itertools.chain(filenames, dirnames):
            if dirpath != ".":
                f = os.path.normpath(os.path.join(dirpath, f))
            match = re.match(pattern, f)
            if match:
                wildcards = snakemake.io.Namedlist(fromdict=match.groupdict())
                yield f, wildcards


def main(cpg_report_regex, outdir, covariates_tsv, test_covariate,
         merge_complementary_cpgs, bsseq_threads, dmrseq_threads, max_gap,
         min_cpgs, min_samples, cutoff, smooth, bp_span, min_in_span,
         max_gap_smooth, max_perms, adjust_covariates, match_covariate,
         ref_build, goldmine_dir, goldmine_sync):
    """
    Create configuration file: documentation at the top of the module.
    """
    # Load covariates
    df_covariates = pandas.read_csv(covariates_tsv, header=0, sep="\t")

    # Check number of columns
    if df_covariates.shape[1] < 2:
        raise ValueError(
            "The covariate file should have at least 2 columns (sample IDs "
            "and test covariate) tab-separated: %s" % covariates_tsv)

    # List names of the covariates to use in the analysis
    covariates_to_use = [test_covariate] + (adjust_covariates or [])
    if match_covariate is not None:
        covariates_to_use += [match_covariate]

    # Check that these covariates are defined in the covariate file
    # Excluding the first column which is assumed to provide the sample IDs
    tsv_covariates = set(df_covariates.columns[1:])
    for cov in covariates_to_use:
        if cov not in tsv_covariates:
            raise ValueError("Covariate %s not found in the covariate file: %s"
                             % (cov, covariates_tsv))

    # Sample IDs in the covariate file
    tsv_samples = set(df_covariates.iloc[:, 0].tolist())

    # Dict mapping <sample> -> <Bismark coverage filepath>
    cpg_report_of_sample = dict()

    # Get couples (<file>, <wildcards>)
    couples = list(listfiles(cpg_report_regex))
    if len(couples) == 0:
        raise ValueError(
            "Found 0 Bismark CpG report file using Snakemake regex: %s"
            % cpg_report_regex)

    # For each couple (Bismark CpG report file, wildcards)
    for cpg_report, wildcards in couples:
        sample = wildcards.sample

        # Check for sample appearing twice
        if sample in cpg_report_of_sample:
            raise ValueError(
                "Sample %s has 2 Bismark CpG report files: %s and %s"
                % (sample, cpg_report, cpg_report_of_sample[sample]))

        # Check if the sample is present in the covariate file, if not
        # print a message indicating that this sample and related Bismark
        # file are ignored in the analysis.
        if sample not in tsv_samples:
            print("Sample ignored: %s, not in covariate file: %s. Related "
                  "Bismark file: %s" % (sample, covariates_tsv, cpg_report))
            continue

        cpg_report_of_sample[sample] = cpg_report

    # Create <outdir> if not existing
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Replace relative paths by absolute paths
    outdir = os.path.abspath(outdir)

    # Create dict to be stored in the YAML config file
    config = dict(cpg_report_regex=cpg_report_regex,
                  cpg_report_of_sample=cpg_report_of_sample,
                  outdir=outdir,
                  min_samples=min_samples,
                  covariates_tsv=covariates_tsv,
                  test_covariate=test_covariate,
                  bsseq_threads=bsseq_threads,
                  dmrseq_threads=dmrseq_threads,
                  max_gap=max_gap,
                  merge_complementary_cpgs=merge_complementary_cpgs,
                  min_cpgs=min_cpgs,
                  cutoff=cutoff,
                  smooth=smooth,
                  bp_span=bp_span,
                  min_in_span=min_in_span,
                  max_gap_smooth=max_gap_smooth,
                  max_perms=max_perms,
                  adjust_covariates=adjust_covariates,
                  match_covariate=match_covariate,
                  ref_build=ref_build,
                  goldmine_dir=goldmine_dir,
                  goldmine_sync=goldmine_sync)

    # Create "<outdir>/cluster_logs" directory if not existing
    # This directory has to exist before we run on the cluster.
    cluster_logs_dir = os.path.join(outdir, "cluster_logs")
    if not os.path.isdir(cluster_logs_dir):
        os.mkdir(cluster_logs_dir)

    # Create cluster config dict
    cluster_out = os.path.join(cluster_logs_dir, "{rule}.{wildcards}.stdout")
    cluster_err = os.path.join(cluster_logs_dir, "{rule}.{wildcards}.stderr")
    cluster_config = {
        "__default__": dict(stdout=cluster_out, stderr=cluster_err),
        "create_bsseq_object": dict(
            walltime="02:00:00", mem_gb=40, cpus=bsseq_threads),
        "clustering_of_methylation_ratios": dict(
            walltime="01:00:00", mem_gb=20, cpus=1),
        "dmrseq": dict(
            walltime="01-00:00:00", mem_gb=40, cpus=dmrseq_threads),
        "annotate_dmrs": dict(
            walltime="01:00:00", mem_gb=4, cpus=1),
        "dmr_html_report": dict(
            walltime="01:00:00", mem_gb=2, cpus=1)
    }

    # Output config paths
    config_yaml = os.path.join(outdir, "dmr_config.yaml")
    cluster_yaml = os.path.join(outdir, "dmr_cluster_config.yaml")

    # Write YAML config file
    with open(config_yaml, "w") as f:
        yaml.dump(config, f, default_flow_style=False, indent=4)

    # Write cluster config file
    with open(cluster_yaml, "w") as f_cluster:
        yaml.dump(cluster_config, f_cluster, default_flow_style=False,
                  indent=4)


if __name__ == "__main__":
    kwargs = get_cmd_line_args()
    main(**kwargs)
