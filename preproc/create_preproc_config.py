#!/usr/bin/env python

# Standard import
import sys
import os
import argparse
import logging
import shutil
import subprocess
import yaml

# Third party import
import snakemake

DOC = """
Script that generates the YAML configuration file required to run the
bisulfite preprocessing pipeline. It requires Python >= 3.3 .

The FASTQ files to process along with the <sample ID>, <lane ID> and
<endedness> (e.g. R1, R2) are identified using a Snakemake regex.

'{sample}' is always required. '{end}' is required for paired-end sequencing.
'{lane}' is optional, it is set to 'L1' in output files if not provided
(meaning that there is one unique lane). You can add arbitrary wildcards to
match variable fields in the filepath, they will be ignored.

In paired-end sequencing, the endedness ('{end}' wildcard) is assumed to be
R1/R2 or r1/r2, if other IDs are used like 'forward' and 'reverse', set the
arguments --r1-id and --r2-id. These arguments are case-insensitive.

Examples of regex (--fastq-regex argument):
    /path/to/study/data/{sample}/{sample}_{ignore}_{lane}_{end}_001.fastq.gz
    /path/to/study/data/{sample}_{lane}_{end}.fastq.gz
    /path/to/study/data/{sample}.fq.gz

All parameters and paths to files to process are stored in the output
'preproc_config.yaml' file.
"""

# Exit if Python version is not high enough
if sys.version_info < (3, 3):
    raise Exception("This script requires Python >= 3.3")


# Default values
DEFAULT_OF = dict(
    bismark_bowtie2=dict(
        seed_mismatch=1,
        threads=4
    ),
    bismark_meth_extract=dict(
        threads=4
    ),
    fastqc=dict(
        threads=4
    ),
    phred=33,
    r1_id="R1",
    r2_id="R2",
    picard=dict(
        jvm_args="-Xmx6g -Xms1g"
    ),
    samtools=dict(
        threads=4
    ),
    trim_galore=dict(
        error_rate=0.1,
        min_length=50,
        quality=20,
        stringency=1
    )
)


def is_file(filepath):
    """ Check existence of file. """
    if not os.path.isfile(filepath):
        raise argparse.ArgumentTypeError("File does not exist: %s" % filepath)
    return filepath


def is_dir(dirpath):
    """ Check existence of directory. """
    if not os.path.isdir(dirpath):
        raise argparse.ArgumentTypeError(
            "Directory does not exist: %s" % dirpath)
    return dirpath


def is_fastq_regex(fastq_regex):
    """
    Check that the provided pattern contains the required wildcards.
    """
    if "{sample}" not in fastq_regex:
        raise argparse.ArgumentTypeError(
            "The regex should contain at least the wildcard '{sample}'.")
    return fastq_regex


def is_ref_genome_dir(dirpath, ref_build):
    """ Check that the directory is a valid reference folder. """
    # Paths relative to directory that should exist
    required_files = [
        "%s.fa" % ref_build, "%s.fa.dict" % ref_build, "%s.fa.fai" % ref_build,
        "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
        "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"]
    for relative_path in required_files:
        path = os.path.join(dirpath, relative_path)
        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(
                "Reference genome folder is invalid: %s, missing file %s"
                % (dirpath, path))
    return dirpath


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
        "--fastq-regex", required=True, type=is_fastq_regex, metavar="<regex>",
        help=
        "Snakemake regex used to infer the FASTQ files to process and the "
        "related wildcards: {sample} (mandatory), {lane} (optional) and "
        "{end} (mandatory if paired-end sequencing), e.g. "
        "/path/to/data/{sample}/{sample}_{ignore1}_{lane}_{end}_001.fastq.gz"
    )
    required.add_argument(
        "--outdir", required=True, metavar="<dir>", help="Output directory."
    )
    required.add_argument(
        "--ref-build", required=True, metavar="<version>", help=
        "Reference genome build, e.g. hg19 or mm10. Assuming the existence of "
        "the 3 following files in <--ref-genome-dir>: <ref build>.fa "
        "<ref build>.fa.fai and <ref build>.fa.dict"
    )
    required.add_argument(
        "--ref-genome-dir", metavar="<dir>", help=
        "Bisulfite reference genome directory, including '<ref build>.fa', "
        "'<ref build>.fa.fai', '<ref build>.fa.dict' and the "
        "'Bisulfite_Genome' directory created by running the "
        "'bismark_genome_preparation' script. See README.md documentation."
    )


    # Optional general arguments
    optional = parser.add_argument_group("Optional")
    optional.add_argument(
        "--single-end", action="store_false", dest="paired_end", help=
        "By default paired-end sequencing is assumed, for single-end set this "
        "flag."
    )
    optional.add_argument(
        "--rrbs", action="store_true", help=
        "For Reduced Representation Bisulfite Sequencing (RRBS) set this flag."
    )
    optional.add_argument(
        "--no-deduplication", action="store_false",
        dest="use_bismark_deduplicate", help=
        "Set this flag to not apply Bismark BAM deduplication. The deduplica"
        "tion removes reads with similar start/end positions on a given chromo"
        "some. It is not a valid PCR correction for RRBS or amplicon data. "
        "The deduplication is not applied if the --rrbs flag is set."
    )
    optional.add_argument(
        "--non-directional-library", action="store_false",
        dest="directional_library", help=
        "By default the library is assumed to be directional, if not set this "
        "flag. See Bismark documentation for more information."
    )
    optional.add_argument(
        "--target-bed", type=is_file, metavar="<path>", help=
        "For targeted sequencing, the path to the BED file listing the regions"
        " targeted. Used only for read coverage computation. If no BED is pro"
        "vided the coverage will be computed on the whole reference genome."
    )
    optional.add_argument(
        "--target-kit", metavar="<name>", help=
        "For targeted sequencing, the name of the kit used to target to be "
        "reported in the preprocessing report. Does not affect processing."
    )
    optional.add_argument(
        "--phred", type=int, choices={33, 64}, default=DEFAULT_OF["phred"],
        metavar="<33|64>", help=
        "Base quality encoding of input FASTQ files: 33|64, by default %i."
        % DEFAULT_OF["phred"]
    )
    optional.add_argument(
        "--r1-id", default=DEFAULT_OF["r1_id"], metavar="<ID>", help=
        "Case-insensitive ID used to identify R1 (forward) reads in paired-end"
        " sequencing, by default '%s'." % DEFAULT_OF["r1_id"]
    )
    optional.add_argument(
        "--r2-id", default=DEFAULT_OF["r2_id"], metavar="<ID>", help=
        "Case-insensitive ID used to identify R2 (reverse) reads in paired-end"
        " sequencing, by default '%s'." % DEFAULT_OF["r2_id"]
    )
    optional.add_argument(
        "--read-length", type=int, metavar="<int>", help=
        "Length of reads (e.g. 150) to write in the HTML report. "
        "Does not affect the processing."
    )

    # Optional FastQC arguments
    fastqc = parser.add_argument_group("FastQC optional")
    optional.add_argument(
        "--fastqc-threads", type=int, metavar="<int>",
        default=DEFAULT_OF["fastqc"]["threads"], help=
        "FastQC '--threads' argument, by default %i."
        % DEFAULT_OF["fastqc"]["threads"]
    )

    # Optional Trim Galore arguments
    trim_galore = parser.add_argument_group("Trim Galore optional")
    ADAPTERS_URL = (
        "https://support.illumina.com/bulletins/2016/12/what-sequences-do-i"
        "-use-for-adapter-trimming.html")
    trim_galore.add_argument(
        "--adapter-r1", metavar="<sequence>", help=
        "Trim Galore '--adapter' argument: adapter sequence to be trimmed off "
        "read 1. Common sequences: %s" % ADAPTERS_URL
    )
    trim_galore.add_argument(
        "--adapter-r2", metavar="<sequence>", help=
        "Trim Galore '--adapter2' argument: adapter sequence to be trimmed "
        "off read 2. Common sequences: %s" % ADAPTERS_URL
    )
    trim_galore.add_argument(
        "--quality", type=int, default=DEFAULT_OF["trim_galore"]["quality"],
        metavar="<int>", help=
        "Trim Galore '--quality' argument, by default %i."
         % DEFAULT_OF["trim_galore"]["quality"]
    )
    trim_galore.add_argument(
        "--stringency", type=int, metavar="<int>",
        default=DEFAULT_OF["trim_galore"]["stringency"], help=
        "Trim Galore '--stringency' argument: overlap with adapter sequence "
        "required to trim, by default %i (very stringent)."
        % DEFAULT_OF["trim_galore"]["stringency"]
    )
    trim_galore.add_argument(
        "--min-length", type=int, metavar="<int>",
        default=DEFAULT_OF["trim_galore"]["min_length"], help=
        "Trim Galore '--length' argument: minimum read length after trimming "
        "otherwise removed, by default %i."
        % DEFAULT_OF["trim_galore"]["min_length"]
    )
    trim_galore.add_argument(
        "--error-rate", type=float, metavar="<float>",
        default=DEFAULT_OF["trim_galore"]["error_rate"], help=
        "Trim Galore '-e' argument: maximum allowed error rate with the "
        "matching region, by default {}"
        .format(DEFAULT_OF["trim_galore"]["error_rate"])
    )
    trim_galore.add_argument(
        "--max-n", type=int, metavar="<int>", help=
        "Trim Galore '--max_n' argument: Maximum number of 'N's in a read "
        "otherwise removed. By default not applied."
    )
    trim_galore.add_argument(
        "--trim-n", action="store_true", help=
        "Trim Galore '--trim-n' argument: remove 'N's from ends of the read."
    )
    trim_galore.add_argument(
        "--clip-r1-5p", type=int, metavar="<int>", help=
        "Trim Galore '--clip_R1' argument: remove basepairs from 5' end of "
        "read 1. Useful if there is a methylation bias at this end."
    )
    trim_galore.add_argument(
        "--clip-r2-5p", type=int, metavar="<int>", help=
        "Trim Galore '--clip_R2' argument: remove basepairs from 5' end of "
        "read 2. Useful if there is a methylation bias at this end."
    )
    trim_galore.add_argument(
        "--clip-r1-3p", type=int, metavar="<int>", help=
        "Trim Galore '--three_prime_clip_R1' argument: remove basepairs from "
        "3' end of read 1. Useful if there is a methylation bias at this end."
    )
    trim_galore.add_argument(
        "--clip-r2-3p", type=int, metavar="<int>", help=
        "Trim Galore '--three_prime_clip_R2' argument: remove basepairs from "
        "3' end of read 2. Useful if there is a methylation bias at this end."
    )

    # Optional Bismark tools arguments
    bismark = parser.add_argument_group("Bismark optional")
    bismark.add_argument(
        "--seed-mismatch", type=int, choices=[0, 1], metavar="<0|1>",
        default=DEFAULT_OF["bismark_bowtie2"]["seed_mismatch"], help=
        "Maximum number of mismatch allowed in a seed alignment: 0|1, "
        "by default %i." % DEFAULT_OF["bismark_bowtie2"]["seed_mismatch"]
    )
    bismark.add_argument(
        "--bowtie2-threads", type=int, metavar="<int>",
        default=DEFAULT_OF["bismark_bowtie2"]["threads"], help=
        "Bowtie2 '--threads' argument, by default %i"
        % DEFAULT_OF["bismark_bowtie2"]["threads"])
    bismark.add_argument(
        "--meth-extract-threads", type=int, metavar="<int>",
        default=DEFAULT_OF["bismark_meth_extract"]["threads"], help=
        "bismark_methylation_extractor '--multicore' argument, by default %i."
        % DEFAULT_OF["bismark_meth_extract"]["threads"]
    )

    # Optional Picard arguments
    picard = parser.add_argument_group("Picard optional")
    picard.add_argument(
        "--picard-jvm-args", default=DEFAULT_OF["picard"]["jvm_args"],
        metavar="<args>", help=
        "Java virtual machine arguments, e.g. to control starting and maximum "
        "heap size when running Picard, by default '%s'."
        % DEFAULT_OF["picard"]["jvm_args"]
    )

    # Optional Samtools arguments
    samtools = parser.add_argument_group("Samtools optional")
    samtools.add_argument(
        "--samtools-threads", type=int, metavar="<int>",
        default=DEFAULT_OF["samtools"]["threads"], help=
        "Samtools '--threads' argument, by default %i"
        % DEFAULT_OF["samtools"]["threads"])

    # Parse the command line
    args = parser.parse_args()

    # For paired-end sequencing, check that the '{end}' wildcard is provided
    if args.paired_end is True and "{end}" not in args.fastq_regex:
        raise ValueError(
            "The wildcard '{end}' is required in --fastq-regex argument when "
            "working with paired-end sequencing.")

    # Set 'use_bismark_deduplicate' to False if RRBS data
    if args.rrbs is True:
        args.use_bismark_deduplicate = False

    # Check reference genome directory
    is_ref_genome_dir(args.ref_genome_dir, args.ref_build)

    # Convert the argparse object to a dict mapping <arg name> -> <val>
    kwargs = vars(args)

    return kwargs


# Function copied from snakemake.utils.listfiles but modified to make it work
# with symbolic links:
# os.walk(dirname) -> os.walk(dirname, followlinks=True)
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


def main(fastq_regex, outdir, ref_build, ref_genome_dir, paired_end, rrbs,
         use_bismark_deduplicate, target_kit, target_bed, read_length, phred,
         r1_id, r2_id, directional_library, adapter_r1, adapter_r2,
         quality, stringency, min_length, error_rate, max_n, trim_n,
         clip_r1_5p, clip_r2_5p, clip_r1_3p, clip_r2_3p, seed_mismatch,
         fastqc_threads, bowtie2_threads, meth_extract_threads,
         samtools_threads, picard_jvm_args):
    """
    Create configuration file: documentation at the top of the module.
    """
    # Dict mapping [<sample>][<lane>][<end>] -> FASTQ file path
    fastq_of = dict()

    # Get couples (<file>, <wildcards>)
    couples = list(snakemake.utils.listfiles(fastq_regex))
    if len(couples) == 0:
        raise ValueError(
            "Found 0 FASTQ files using the Snakemake regex: %s" % fastq_regex)

    # For each couple (FASTQ, wildcards)
    for fastq_path, wildcards in couples:
        sample = wildcards.sample

        # If paired-end sequencing
        if paired_end is True:
            end = wildcards.end

            # Check 'end' ID and normalize as R1/R2
            if end.lower() == r1_id.lower():
                end = "R1"
            elif end.lower() == r2_id.lower():
                end = "R2"
            else:
                raise ValueError(
                    "Bad 'end' ID: %s for FASTQ: %s, should be %s or %s"
                    % (end, fastq_path, r1_id, r2_id))

        # If single-end sequencing
        else:
            end = "R1"

        # If no lane, assume a unique lane called "L1"
        lane = wildcards.get("lane", "L1")

        # Update 'fastq_of' dict
        if sample not in fastq_of:
            fastq_of[sample] = {lane: {end: fastq_path}}
        elif lane not in fastq_of[sample]:
            fastq_of[sample][lane] = {end: fastq_path}
        elif end not in fastq_of[sample][lane]:
            fastq_of[sample][lane][end] = fastq_path
        else:
            raise ValueError(
                "The combination <sample> and <end> (and <lane> if provided) "
                "should be unique. Violated by %s and %s"
                % (fast_path, fastq_of[sample][lane][end]))

    # If paired-end sequencing check for unpaired FASTQ files
    if paired_end is True:
        for sample in fastq_of:
            for lane in fastq_of[sample]:
                nb_fastq_files = len(fastq_of[sample][lane])
                if nb_fastq_files == 1:
                    unpaired_fastq = fastq_of[sample][lane].values()[0]
                    raise Exception(
                        "Unpaired FASTQ file: %s . If you are using single-end"
                        " sequencing set the '--single-end' flag."
                        % unpaired_fastq)
                elif nb_fastq_files > 2:
                    raise Exception(
                        "Found {} FASTQ files when expecting 2, for the couple"
                        " (<sample>, <lane>): ({}, {}), FASTQ files: {}"
                        .format(nb_fastq_files, sample, lane,
                                fastq_of[sample][lane].values()))

    # Create <outdir> if not existing
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Write absolute paths in the configuration file
    outdir = os.path.abspath(outdir)

    # Create config dict
    config = dict(
        bismark_bowtie2=dict(
            seed_mismatch=seed_mismatch,
            threads=bowtie2_threads,
        ),
        bismark_meth_extract=dict(
            threads=meth_extract_threads
        ),
        target_kit=target_kit,
        directional_library=directional_library,
        fastq_regex=fastq_regex,
        fastq_of=fastq_of,
        fastqc=dict(
            threads=fastqc_threads
        ),
        outdir=outdir,
        paired_end=paired_end,
        picard=dict(
            jvm_args=picard_jvm_args
        ),
        phred=phred,
        rrbs=rrbs,
        read_length=read_length,
        ref_build=ref_build,
        ref_genome_dir=ref_genome_dir,
        samtools=dict(
            threads=samtools_threads
        ),
        target_bed=target_bed,
        trim_galore=dict(
            adapter_r1=adapter_r1,
            adapter_r2=adapter_r2,
            quality=quality,
            stringency=stringency,
            min_length=min_length,
            error_rate=error_rate,
            max_n=max_n,
            trim_n=trim_n,
            clip_r1_5p=clip_r1_5p,
            clip_r2_5p=clip_r2_5p,
            clip_r1_3p=clip_r1_3p,
            clip_r2_3p=clip_r2_3p
        ),
        use_bismark_deduplicate=use_bismark_deduplicate
    )

    # Create "<outdir>/cluster_logs" directory if not existing
    cluster_logs_dir = os.path.join(outdir, "cluster_logs")
    if not os.path.isdir(cluster_logs_dir):
        os.mkdir(cluster_logs_dir)

    # Create cluster config dict
    cluster_out = os.path.join(cluster_logs_dir, "{rule}.{wildcards}.stdout")
    cluster_err = os.path.join(cluster_logs_dir, "{rule}.{wildcards}.stderr")
    cluster_config = {
        "__default__": dict(stdout=cluster_out, stderr=cluster_err),
        "fastqc": dict(
            walltime="01:00:00", mem_gb=4, cpus=fastqc_threads),
        "trim_galore": dict(
            walltime="08:00:00", mem_gb=1, cpus=2),
        "fastqc_post_trimming": dict(
            walltime="01:00:00", mem_gb=4, cpus=fastqc_threads),
        "bismark_bowtie2": dict(
            walltime="24:00:00", mem_gb=20, cpus=bowtie2_threads),
        "samtools_merge_lanes": dict(
            walltime="01:00:00", mem_gb=4, cpus=samtools_threads),
        "bismark_deduplicate": dict(
            walltime="02:00:00", mem_gb=12, cpus=1),
        "bismark_methylation_extractor": dict(
            walltime="08:00:00", mem_gb=8, cpus=meth_extract_threads),
        "samtools_sort_by_position": dict(
            walltime="01:00:00", mem_gb=12, cpus=samtools_threads),
        "format_target_bed_for_picard_hs_metrics": dict(
            walltime="01:00:00", mem_gb=8, cpus=1),
        "picard_hs_metrics": dict(
            walltime="01:00:00", mem_gb=8, cpus=1),
        "picard_wgs_metrics": dict(
            walltime="01:00:00", mem_gb=8, cpus=1),
        "multiqc": dict(
            walltime="00:30:00", mem_gb=8, cpus=1),
        "preproc_html_report": dict(
            walltime="00:10:00", mem_gb=1, cpus=1)
    }

    # Output config paths
    config_yaml = os.path.join(outdir, "preproc_config.yaml")
    cluster_yaml = os.path.join(outdir, "preproc_cluster_config.yaml")

    # Write config file
    with open(config_yaml, "w") as f_config:
        yaml.dump(config, f_config, default_flow_style=False, indent=4)

    # Write cluster config file
    with open(cluster_yaml, "w") as f_cluster:
        yaml.dump(cluster_config, f_cluster, default_flow_style=False,
                  indent=4)


if __name__ == "__main__":
    kwargs = get_cmd_line_args()
    main(**kwargs)
