# mpileup_call_pipeline

A snakemake pipeline for calling SNP genotypes at known positions from short-read data using `bcftools mpileup`.

Here are the main steps:

1. Trim adapters with `trim_galore`, and including running `fastqc` on the output.
2. Collate `fastqc` results into a single report with `multiqc`.
3. Align reads to a FASTA genome file with `bwa-mem`.
4. Identify and remove reads that are likely PCR duplicates with `samtools`.
5. Call genotypes at known SNP positions using `bcftools mpileup` and `bcftools call`. This is done separately for each chromosome.
6. Merge VCF files from each chromosome.
7. Rename samples in the VCF file to match names from the sample sheet.

## Contents

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [mpileup\_call\_pipeline](#mpileup_call_pipeline)
  - [Contents](#contents)
  - [Installation](#installation)
  - [Clone this repo](#clone-this-repo)
    - [Dependencies](#dependencies)
      - [As a conda environment](#as-a-conda-environment)
      - [With `module load`](#with-module-load)
  - [Input data](#input-data)
    - [Config file](#config-file)
    - [Raw sequence data](#raw-sequence-data)
    - [Reference genome](#reference-genome)
    - [Positions of known variable sites](#positions-of-known-variable-sites)
      - [Ready-made files](#ready-made-files)
      - [Create a targets file yourself](#create-a-targets-file-yourself)
  - [Usage](#usage)
    - [Define paths](#define-paths)
    - [Run the pipeline](#run-the-pipeline)
  - [Output files](#output-files)
  - [Acknowledgements](#acknowledgements)

<!-- TOC end -->

## Installation

## Clone this repo

Clone the repo to your project folder:
```sh
git clone https://github.com/ellisztamas/mpileup_call_pipeline.git
```

### Dependencies

#### As a conda environment

A conda environment is provided to install the necessary dependencies.
This will create a conda environment called `mpileup_call_pipeline`.

```sh
cd mpileup_call_pipeline
mamba env create -f environment.yml
```

I recommend using mamba instead of conda to install the environment.

#### With `module load`

On the VBC Clip cluster you can load most of the dependencies with `module load`:

```sh
module load  build-env/f2022
module load bcftools/1.17-gcc-12.2.0
module load bwa/0.7.17-gcccore-12.3.0
module load multiqc/1.14-foss-2021b
module load samtools/1.18-gcc-12.3.0
module load snakemake/9.5.1-foss-2023b
module load trim_galore/0.6.10-gcccore-12.2.0
```

The [snakemake-executor-plugin-slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) dependency is needed (probably?) to run
the pipeline with SLURM, which is not available with `module load`.
You can install it directly:

```sh
pip install snakemake-executor-plugin-slurm
```

## Input data

### Config file

There is a config file `config.yaml` in the same directory as the snakefile.
It specifies the paths to the data files.
It also defines a "project_name" which will be the name of the final VCF file.
If you run the pipeline more than once, copy the config file for each analysis and save the changes.

### Raw sequence data

The pipeline uses a comma-separated sample sheet giving:

1. sample name
2. paths to fastq files for mate pair 1.
3. paths to fastq files for mate pair 2.

The path should be defined in the config file.

For example:
```
sample,fastq1,fastq2
bob,/rawdata/rawdata_001_R1.fastq,/rawdata/rawdata_001_R2.fastq
alice,/rawdata/rawdata_002_R1.fastq,/rawdata/rawdata_002_R2.fastq
steve,/rawdata/rawdata_003_R1.fastq,/rawdata/rawdata_003_R2.fastq
```
It is probably best if paths are given as absolute paths.

If you're using data from the VBC NGS facility, you may be able to use [this tool](https://methlab.readthedocs.io/en/latest/modules/align_plate_positions.html) to create that sample sheet.

### Reference genome

Give the path to a reference genome to which reads will be aligned in FASTA format in the config file.
This should be indexed with `samtools faidx`.

### Positions of known variable sites

You also need to provide a 'targets file' giving the positions of known SNP positions in the reference genome.
Define the path to the targets file in the config file.

#### Ready-made files

For internal use in the Nordborg group, there are two ready-made target files you may wish to use if you work on *Arabidopsis thaliana* and are aligning reads to the TAIR10 assembly:
```
/groups/nordborg/common_data/fernandos_SNP_calls/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.targets_file.tsv.gz
/groups/nordborg/common_data/1001genomes/1135g_SNP_BIALLELIC.tsv.gz
```
You probably want to use the first one, and this is the default path in the config file.

#### Create a targets file yourself

This should be formatted as a file that can be passed to the argument `-T` or ``--targets-file` in `bcftools` (see the [help page](https://samtools.github.io/bcftools/bcftools.html#common_options)).
It should be tab separated with a row for each marker, and columns giving chromosome/contig, position, and a comma-separated list of alleles starting with the reference allele.
It should also be gzipped and have no header.
Here is an example of how that might look:

```ts
Chr1	167	T,A
Chr1	276	T,G
Chr1	314	C,T
Chr1	322	G,A
Chr1	323	G,A
```

Using the example from the `bcftools` documentation you can easily create a file like this using a command similar to the following:
```sh
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' file.vcf | bgzip -c > als.tsv.gz && tabix -s1 -b2 -e2 als.tsv.gz
```

## Usage

The pipeline is designed to submit individual steps to the CLIP cluster via the SLURM scheduler, so you only need to run the command to start the pipeline in a terminal (see below).

However, as the pipeline will take a long time, you should probably run this inside a `tmux` window so that the window stays active, even if your local machine goes to sleep.
([Here](https://www.howtogeek.com/671422/how-to-use-tmux-on-linux-and-why-its-better-than-screen/) is a tutorial on getting started with `tmux`).

### Define paths

Activate the conda environment.
```sh
conda activate mpileup_call_pipeline
```

Define paths to snakemake file, config file and where to store the output.
Change these paths as necessary.

```sh
# Input files
pipeline=$PWD/02_library/mpileup_call_pipeline/snakefile.smk
# Config file giving 
config_file=$PWD/02_library/mpileup_call_pipeline/config.yaml
# Directory to store the output
outdir=/scratch-cbe/users/$(whoami)/my_genotyping
```

On the VBC cluster I strongly recommend setting an output directory on `scratch-cbe`.
The pipeline will create a lot of intermediate files that you probably don't need, but do take up a lot of space.
When the pipeline is created, you can just copy the final VCF file out and remove the rest.
The `whoami` command will insert your username and ensure the path works.

### Run the pipeline

Run the pipeline via SLURM:

```sh
snakemake \
    --snakefile "$pipeline" \
    --configfile "$config_file" \
    --directory "$outdir" \
    --cores all \
    --executor slurm \
    --rerun-incomplete \
    --restart-times 2 \
    -j5 
```

Some explanations:

* The convoluted arguments to `--config` pass the variables defined above to snakemake as wildcards. You could also put these in a config file, but this could get messy if you need to run the pipeline several times with different settings and aren't careful.
* The argument `-j N` tells the pipeline to run N jobs in parallel, assuming there are N input files. I used five as a place to start because the genotype calls are done by chromosome, and *Arabidopsis thaliana* has five chromosomes. Increase the number if you want to run more samples in parallel (for example, `-j 96` for a 96-well plate).
* `--rerun-incomplete` tells the pipeline to start any failed steps from scratch if you encounter and error and need to come back to it.
* `--restart-times 2` tells the pipeline to try failed jobs again. For a couple of resource intensive steps, it wil try again with additional memory and time.

Change these inputs as necessary.

## Output files

The pipeline outputs one directory per rule.
It should be fairly clear what those files and logs are, but you can probably ignore them unless you need to debug something.

The main output file is a VCF file in the root folder of the output directory, which will be called the project name you passed to the pipeline, followed by `.vcf.gz`.
There should be an accompanying index file.
If the pipeline worked, these files are probably the only ones you need to keep.

You may also want to save the `multiqc` report which is run on fastq files after trimmin, found in the directory `trim_galore`.

## Acknowledgements

This was adapted from shell scripts by Pieter Clauw.
I used GPT 5.1 and Claude Sonnet for drafting rules.
I created the TOC for this README using https://bitdowntoc.derlin.ch/.