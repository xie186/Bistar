
## Bistar

### Description

This repository provides 2 Snakemake pipelines for identification of CpG differentially
methylated regions from bisulfite next-generation sequencing data. The first pipeline
preprocesses raw FASTQ files and extracts methylation counts at each CpG. The second
pipeline computes and annotates the differentially methylated CpG regions from
the CpG methylation counts.

![flowchart](flowchart.svg)

The pipelines are built around [Conda](https://conda.io/docs/) and
[Snakemake](https://snakemake.readthedocs.io/en/stable/) and can only be used
with Linux whether you are working with a workstation or a cluster.

### Conda

[Conda](https://conda.io/docs/) is used to make installation of pipelines and their
dependencies automatical, with controlled versions and without admin privileges.
If you don't have Conda installed, you can install [Miniconda3](https://conda.io/miniconda.html),
a lightweight version of Conda with [Python3](https://www.python.org/). Miniconda3
can be installed by downloading and running an installation script available at
<https://conda.io/miniconda.html>.

### Snakemake

[Snakemake](https://snakemake.readthedocs.io/en/stable/) is a workflow manager
that makes the pipeline scalable on both a single workstation or a cluster. It
provides other functionalities like allowing to run only selected steps of the
pipeline, or not re-running steps that had successfully run with the same parameters.

### How to install
1. Activate Conda's base environment if not already the case.
```bash
# Replace '~/miniconda3/' with your Conda's installation directory
source ~/miniconda3/bin/activate
```

2. Create the processing environment, here called Bistar, and install the pipelines
with all their dependencies. It can take a while.
```bash
conda create -n Bistar -c icm-iconics bistar
```

3. Activate the environment
```bash
conda activate Bistar
```

The $BISTAR_DIR environment variable is defined and points to the installation
directory.

If you want to leave an environment
```bash
conda deactivate
```

[Preprocessing pipeline](preproc/README.md)

[Differential methylation analysis pipeline](dmr/README.md)
