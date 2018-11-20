# Bistar: bisulfite-sequencing preprocessing pipeline

## Table of contents
* [Description](#description)
    * [Steps of the pipeline](#steps-of-the-pipeline)
* [Genome of reference](#genome-of-reference)
* [Creation of configuration files](#creation-of-configuration-files)
    * [Options](#options)
    * [Example of configuration file creation](#example-of-configuration-file-creation)
* [How to run](#how-to-run)
    * [Running on a workstation](#running-on-a-workstation)
    * [Running on a SLURM cluster](#running-on-a-slurm-cluster)
* [Report](#report)

## Description

The first pipeline preprocesses bisulfite sequencing data from FASTQ files to
methylation calls at each CpG. The pipeline should be applicable to most approaches
to bisulfite sequencing as long as it is Next-Generation Sequencing (NGS) data
including whole genome sequencing, targeted sequencing (by capture or amplicon-based)
and reduced representation bisulfite sequencing (RRBS), and should work whether
it is single-end or paired-end sequencing. The pipeline is not meant to work with
single cell data.

Note that the PCR bias correction step removes reads with similar genomic positions
(i.e. with the same chromosome, start and end positions) and should not be applied
if the start/end positions of reads are not random. Meaning that the PCR bias
cannot be corrected for RRBS or amplicon data using this pipeline.

### Files

The pipeline is provided as 4 files:

| File | Description |
| ---- | ----------- |
| [README.md](README.md) | This documentation. |
| [create_preproc_config.py](create_preproc_config.py) | Script that helps creating the pipeline's configuration files so that the user can specify what to process and how. |
| [Snakefile_preproc](Snakefile_preproc) | the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow. |
| [preproc_report.Rmd](preproc_report.Rmd) | [R Markdown](https://rmarkdown.rstudio.com/) file used to generate a report at the end of the pipeline. |

### Steps of the pipeline:

1. QC raw reads with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
2. Run [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
   to remove a few basepairs at the read ends (if requested), trim low quality
   bases at both ends and remove adapters.
3. QC trimmed reads with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
4. Align trimmed reads to the reference genome with [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
  and [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/).
5. Merge lanes after alignment with [Samtools](http://www.htslib.org/) if there
   are multiple lanes.
6. Remove duplicates (PCR bias correction) with [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
   if requested (i.e. if --no-deduplication or --rrbs flags were not used when
   calling create_preproc_config.py).
7. Call methylation with [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
   (methylation counts at each CpG).
8. Sort BAM files by genomic position with [Samtools](http://www.htslib.org/)
   (sorting is done here because Bismark's steps require the BAM to be sorted by
   read tag).
9. Compute read coverage with [Picard](https://broadinstitute.github.io/picard/)
   using [CollectHsMetrics](http://broadinstitute.github.io/picard/picard-metric-definitions.html#HsMetrics)
   if a target BED is provided otherwise with
   [CollectWgsMetrics](http://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetrics.WgsMetrics).
10. Create plots with [MultiQC](http://multiqc.info/).
11. Create a preprocessing report using [R Markdown](https://rmarkdown.rstudio.com/)
    and [knitr](https://yihui.name/knitr/).


## Genome of reference

The alignment of reads is done with Bismark and requires a 'bisulfited' reference
genome. You can create one by providing your reference FASTA file to the
*bismark_genome_preparation* script.

### Example of creation for hg19 from scratch

The proposed example is for hg19 but could be easily translated to other organisms
or human builds, in particular if the reference FASTA files are available on the
[UCSC FTP server](http://hgdownload.cse.ucsc.edu/downloads.html).

Assuming that you are in the *Bistar* Conda environment, the reference genome is
created in the current directory.

1. Download the TAR with all the FASTA files of the organism, here from the UCSC
   FTP server:
```bash
TAR_URL=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
wget --no-clobber ${TAR_URL} --tries=5 --continue --timeout=30
```

2. Extract the FASTA files of interest from archive (here *chromFa.tar.gz*) and
   merge them, keeping the order, as one unique FASTA file: *&lt;ref build&gt;.fa*
```bash
echo chr{1..22}.fa chrX.fa chrY.fa chrM.fa | xargs -n1 tar -zxvf chromFa.tar.gz --to-stdout > hg19.fa
```

3. Create the index of your FASTA file with Samtools: *&lt;ref build&gt;.fa.fai*
```bash
samtools faidx hg19.fa
```

4. Create the sequence dictionary with Picard of your FASTA file:
   *&lt;ref build&gt;.fa.dict*
```bash
picard CreateSequenceDictionary REFERENCE=hg19.fa OUTPUT=hg19.fa.dict
```

5. Run the Bismark script: *bismark_genome_preparation*. Note that this step can
   take more than an hour and consumes a lot of RAM (can be >20GB depending on
   the size of your FASTA file).
```bash
bismark_genome_preparation --bowtie2 .
```

The pipeline will expect to find the 3 files named *&lt;ref build&gt;.fa*,
*&lt;ref build&gt;.fa.fai* and *&lt;ref build.fa.dict&gt;* in your reference genome
directory and a sub-directory named *Bisulfite_Genome* created by the
*bismark_genome_preparation* script.

## Creation of configuration files

The pipeline works with 2 YAML configuration files that are created by the provided
Python script (*create_preproc_config.py*):
- *preproc\_config.yaml*: general configuration file with the preprocessing
  parameters, the files to process and the associated metadata like sample ID,
  lane ID and endedness (forward or reverse reads if paired-end sequencing).
- *preproc\_cluster\_config.yaml*: used if you run on a cluster to adapt the amount
  of ressources that are requested to the cluster for each step of the pipeline.
  You can edit it to adapt the number of CPUs, wall time or memory consumption.

### List the files to process

In Snakemake you generally don't explicitly provide a list of files to process
to run the pipeline, instead you assume a regularity in the filesystem's structure
and you identify files and their metadata using a Snakemake regular expression.
It is used to infer the list of files to process and in our case the sample ID,
lane ID (optional) and endedness (for paired-end sequencing) of the reads from
the filepath/filename. You can name the FASTQ as you wish <as long as these
metadata explicitly appear in the filename/filepath.

#### Example 1

Assuming that you want to process these 8 files (paired-end sequencing example):

    /dir1/dir2/S1_L001_R1_001.fq.gz
    /dir1/dir2/S1_L001_R2_001.fq.gz
    /dir1/dir2/S1_L002_R1_001.fq.gz
    /dir1/dir2/S1_L002_R2_001.fq.gz
    /dir1/dir2/S2_L001_R1_001.fq.gz
    /dir1/dir2/S2_L001_R2_001.fq.gz
    /dir1/dir2/S2_L002_R1_001.fq.gz
    /dir1/dir2/S2_L002_R2_001.fq.gz

the regular expression would be:

    /dir1/dir2/{sample}_{lane}_{end}_001.fq.gz

#### Example 2

Assuming that you want to process these 4 files (single-end sequencing example):

    /dir1/dir2/S1/S1.fastq.gz
    /dir1/dir2/S1/S1.fastq.gz
    /dir1/dir2/S2/S2.fastq.gz
    /dir1/dir2/S2/S2.fastq.gz

In that case the sample ID appears twice in the filepath and there is no lane ID
or end ID. The lane is then assumed to be unique and will identified as "L1" in
output files. '{end}' is only mandatory for paired-end sequencing and will be set
to R1 in output files for in single-end sequencing. '{sample}' is always required.
The regular expression would be:

    /dir1/dir2/{sample}/{sample}.fastq.gz

#### Example 3

Assuming that you want to process these 4 files:

    /dir1/dir2/S1_739_R1_999.fq
    /dir1/dir2/S1_739_R2_999.fq
    /dir1/dir2/S2_555_R1_222.fq
    /dir1/dir2/S2_555_R2_222.fq

If there are other variable fields in the filepaths you can define arbitrary wildcards
to match them, as long as they are not named '{sample}', '{lane}' or '{end}'. Here
there are 2 variable fields in the filepaths that are not metadata of interest,
we define 2 wildcards to match them '{ignore1}' and '{ignore2}':

    /dir1/dir2/{sample}_{ignore1}_{end}_{ignore2}.fq

#### Example 4

Assuming that you want to process these 4 files:

    /dir1/dir2/S1_L001_FORWARD.fastq.gz
    /dir1/dir2/S1_L001_REVERSE.fastq.gz
    /dir1/dir2/S2_L001_FORWARD.fastq.gz
    /dir1/dir2/S2_L001_REVERSE.fastq.gz

Sample IDs and lane IDs can be anything but by default the endedness (in paired-end
sequencing) has to be 'R1'/'R2' or 'r1'/'r2'. If not, you have to specify what
IDs identify forward and reverse reads by using the --r1-id and --r2-id arguments
when calling the Python script. The regular expression would be:

    /dir1/dir2/{sample}_{lane}_{end}.fastq.gz

and you would use the following arguments when calling the Python script:

```bash
create_preproc_config.py \
    --fastq-regex /dir1/dir2/{sample}_{lane}_{end}.fq.gz \
    ...
    --r1-id FORWARD \
    --r2-id REVERSE \
    ...
```

### Required settings

| Argument | Description |
| -------- | ----------- |
| --fastq-regex | A Snakemake regular expression see [list files to process](#list-the-files-to-process) |
| --outdir | Output directory |
| --ref-build | Reference genome build, e.g. hg19 or mm10. |
| --ref-genome-dir | Bisulfite reference genome directory, including '&lt;ref build&gt;.fa', '&lt;ref build&gt;.fa.fai', '&lt;ref build&gt;.fa.dict' and the 'Bisulfite_Genome' directory created by running the 'bismark_genome_preparation' script. See README.md documentation. |

### Optional settings

| Argument | Description |
| -------- | ----------- |
| --adapter-r1 | Trim Galore '--adapter' argument: adapter sequence to be trimmed off read 1. [Common sequences](https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html). |
| --adapter-r2 | Trim Galore '--adapter2' argument: adapter sequence to be trimmed off read 2. |
| --quality | Trim Galore '--quality' argument, by default 20. |
| --stringency | Trim Galore '--stringency' argument: overlap with adapter sequence required to trim, by default 1 (very stringent). |
| --min-length | Trim Galore '--length' argument: minimum read length after trimming otherwise removed, by default 50. |
| --error-rate | Trim Galore '-e' argument: maximum allowed error rate with the matching region, by default 0.1 |
| --max-n | Trim Galore '--max_n' argument: Maximum number of 'N's in a read otherwise removed. By default not applied. |
| --trim-n | Trim Galore '--trim-n' argument: remove 'N's from ends of the read. | 
| --clip-r1-5p, --clip-r2-5p, --clip-r1-3p or --clip-r2-3p | Trim Galore '--clip_R1', '--clip_R2', '--three_prime_clip_R1' or '--three_prime_clip_R2' argument: remove basepairs from read ends. Useful if there is a methylation bias at any end. For more information on methylation bias, see: https://sequencing.qcfail.com/articles/library-end-repair-reaction-introduces-methylation-biases-in-paired-end-pe-bisulfite-seq-applications/ . |
| --single-end | By default paired-end sequencing is assumed, for single-end set this flag. |
| --rrbs | For Reduced Representation Bisulfite Sequencing (RRBS) set this flag. |
| --non-directional-library | By default the library is assumed to be directional, if not set this flag. See [Bismark's documentation](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html) for more information. |
| --target-bed | For targeted sequencing, the path to the BED file listing the regions targeted. Used only for read coverage computation. If no BED is provided the coverage will be computed on the whole reference genome. |
| --target-kit | For targeted sequencing, the name of the kit used to target to be reported in the preprocessing report. Does not affect processing. |
| --threads | Number of threads per job, by default 4. |
| --phred | Base quality encoding of input FASTQ files: 33 (default) or 64. |
| --r1-id | Case-insensitive ID used to identify R1 (forward) reads in paired-end sequencing, by default 'R1'. |
| --r2-id | Case-insensitive ID used to identify R2 (reverse) reads in paired-end sequencing, by default 'R2'. |
| --read-length | Length of reads (e.g. 150) to report in the preprocessing report. Does not affect the processing. |
| --max-mismatch | Bismark bowtie2 alignment: maximum number of mismatch allowed in a seed alignment (0 or 1), by default 1. |
| --picard-jvm-args | Java virtual machine arguments, e.g. to control starting and maximum heap size when running Picard, by default '-Xmx6g -Xms1g'. |
| --help | List all the available options that can be set and for optional parameters their default value. |

### Example of configuration file creation

Assuming that you are in the *Bistar* environment.

```bash
# Using a shell variable or adapt the following commands
PREPROC_OUTDIR=/path/to/preproc/outdir

create_preproc_config.py \
    --fastq-regex /path/to/data/{sample}/{sample}_{ignore}_{lane}_{end}_001.fastq.gz \
    --outdir ${PREPROC_OUTDIR} \
    --target-bed ${BISTAR_DIR}/capture_kit/Nimblegen_Epi_SeqCap_CpGiant/130912_HG19_CpGiant_4M_EPI.bed \
    --target-kit "Nimblegen Epi SeqCap CpGiant" \
    --ref-genome-dir /path/to/reference/genome/folder \
    --ref-build hg19 \
    --read-length 150 \
    --threads 4 \
    --phred 33 \
    --r1-id R1 \
    --r2-id R2 \
    --quality 20 \
    --stringency 1 \
    --min-length 50 \
    --error-rate 0.1 \
    --trim-n \
    --max-mismatch 1 \
    --clip-r1-5p 2 \
    --clip-r2-5p 2 \
    --clip-r2-3p 1
```

Here no adapter sequence is provided with *--adapter-r1* (and *--adapter-r2* in
paired-end sequencing), Trim Galore will look for the Illumina universal adapter
sequence with should be enough for most projects sequenced with Illumina. See
[Trim Galore's documentation](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md).

If you don't provide a target BED file the coverage will be computed on the whole
reference genome.

Once the configuration file *preproc_config.yaml* is created, open it to check
that all the files to be processed are listed.

## How to run

Assuming that you are in the *Bistar* environment, that you have the 'bisulfited'
reference genome and the configuration files have been created (at least *preproc_config.yaml*
if you don't run on a cluster). The pipeline is run using Snakemake. Here the
'--jobs' argument sets the maximum number of parallel jobs, could be 2 on a
workstation or 20 on a cluster. Each job computes a step of a the pipeline for
a given input using the number of threads that is specified in the *preproc_config.yaml*
configuration file.

The $BISTAR_DIR environment variable is created automatically when you activate
the Bistar environment.

#### Running on a workstation

```bash
# Dry run to check, remove --dryrun to run
snakemake \
    --configfile ${PREPROC_OUTDIR}/preproc_config.yaml \
    --snakefile ${BISTAR_DIR}/preproc/Snakefile_preproc \
    --jobs 6 \
    --printshellcmds \
    --latency-wait 30 \
    --keep-going \
    --dryrun
```

#### Running on a cluster

The example bellow is for a SLURM cluster, if you are using another scheduler you
have to adapt the '--cluster' argument with the command used on your cluster to
create a job for a shellscript. Snakemake will automatically replace the variables
'{cluster.<variable>}' with the value found in *preproc_cluster_config.yaml* for
each step.

You can edit the *preproc_cluster_config.yaml* file to adapt the number of CPUs,
wall time or RAM memory that are requested to the cluster for each step of the
pipeline. Note that if you change the number of CPUs requested for a given step
it should be consistent with the number of threads used in the command which is
provided by the *preproc_config.yaml* configuration file. The default values for
memory and walltime were set to be largely enough when using 8 threads for all
steps on our internal projects: targeted sequencing with ~80 million basepairs.

Run the complete pipeline
```bash
# Dry run to check, remove --dryrun to run
snakemake \
    --configfile ${PREPROC_OUTDIR}/preproc_config.yaml \
    --cluster-config ${PREPROC_OUTDIR}/preproc_cluster_config.yaml \
    --snakefile ${BISTAR_DIR}/preproc/Snakefile_preproc \
    --jobs 6 \
    --cluster "sbatch --verbose \
                      --account=bioinfo \
                      --partition=normal \
                      --time={cluster.walltime} \
                      --mem={cluster.mem_gb}G \
                      --cpus-per-task={cluster.cpus} \
                      --output={cluster.stdout} \
                      --error={cluster.stderr}" \
    --printshellcmds \
    --latency-wait 60 \
    --keep-going \
    --dryrun
```

## Report
At the end of processing an HTML report file is created (*preprocessing_report.html*)
to help qualifying the results with plots and tables. They may indicate that you
should rerun the preprocessing with different parameters before using the other
pipeline e.g. if there is a methylation bias at read ends that should be trimmed.

## Differential methylation analysis pipeline

[README](../dmr/README.md)
