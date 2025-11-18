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

- [mpileup_call_pipeline](#mpileup_call_pipeline)
   * [Contents](#contents)
   * [Installation](#installation)
   * [Input data](#input-data)
      + [Reference genome](#reference-genome)
      + [Raw sequence data](#raw-sequence-data)
      + [Positions of known variable sites](#positions-of-known-variable-sites)
      + [Config file](#config-file)
   * [Usage](#usage)
      + [Define paths](#define-paths)
      + [Run the pipeline](#run-the-pipeline)
   * [Output files](#output-files)
   * [Acknowledgements](#acknowledgements)

## Installation

Clone the repo to your project folder:
```sh
git clone https://github.com/ellisztamas/mpileup_call_pipeline.git
```

A conda environment is provided to install the necessary dependencies.
This will create a conda environment called `mpileup_call_pipeline`.
I recommend using mamba instead of conda to install the environment.
```sh
cd mpileup_call_pipeline
mamba env create -f environment.yml
```

## Input data

### Reference genome

Give the path to a reference genome to which reads will be aligned in FASTA format.
This should be indexed with `samtools faidx`.

### Raw sequence data

The pipeline uses a comma-separated sample sheet giving

1. sample name
2. paths to fastq files for mate pair 1.
3. paths to fastq files for mate pair 2.

```
sample,fastq1,fastq2
bob,/rawdata/rawdata_001_R1.fastq,/rawdata/rawdata_001_R2.fastq
alice,/rawdata/rawdata_002_R1.fastq,/rawdata/rawdata_002_R2.fastq
steve,/rawdata/rawdata_003_R1.fastq,/rawdata/rawdata_003_R2.fastq
```
It is probably best if paths are given as absolute paths.

If you're using data from the VBC NGS facility, you may be able to use this tool to create that sample sheet.

### Positions of known variable sites

You also need to provide a 'targets file' giving the positions of known SNP positions.
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

Using the example from `bcftools` you can easily create a file like this using a command similar to the following:
```sh
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' file.vcf | bgzip -c > als.tsv.gz && tabix -s1 -b2 -e2 als.tsv.gz
```

### Config file

There is a config file in the same directory as the snakefile.
Currently it only has some options to specify options to trim_galore.
You could edit this if you needed to, but it is probably fine to leave it as is.

## Usage

The pipeline is designed to submit individual steps to the CLIP cluster via the SLURM scheduler, so you only need to run the command to start the pipeline in a terminal (see below).

However, as the pipeline will take a long time, you should probably run this inside a `tmux` window so that the window stays active, even if your local machine goes to sleep.
([Here](https://www.howtogeek.com/671422/how-to-use-tmux-on-linux-and-why-its-better-than-screen/) is a tutorial on getting started with `tmux`).

### Define paths

Activate the conda environment.
```sh
conda activate mpileup_call_pipeline
```

Define paths to input files, pipeline, and where to store the output.
This also gives a `project_name`; the final VCF file will have this name.
```sh
# Input files
sample_sheet=$PWD/02_library/mpileup_call_pipeline/test_sample_sheet.csv
project_name=test_data
targets=$PWD/data/targets_file.tsv.gz
genome=$PWD/data/reference_genome.fas
# Pipeline files
pipeline=$PWD/02_library/mpileup_call_pipeline/snakefile.smk
config_file=$PWD/02_library/mpileup_call_pipeline/config.yaml
# Directory to store the output
outdir=scratchdir/mpileup_call_pipeline
```
I am assuming you are working from the root of a project directory, but are saving to a temporary working directory called `scratchdir` - edit as required.
If it is available on your system, using a separate output directory can be useful because most of the files created will be intermediate files that you probably don't need to keep, but do take up a lot of space.
When the pipeline is created, you can just copy the final VCF file out and remove the rest.
Note the `$PWD`, which ensures paths are defined relative to the current (project) working directory, but allows relative paths to be passed flexibly to snakemake as absolute paths.

### Run the pipeline

Run the pipeline via SLURM:

```sh
snakemake \
    --snakefile $pipeline \
    --configfile $config_file \
    --config sample_sheet=$sample_sheet project_name=$project_name targets=$targets fasta=$genome \
    --directory $outdir \
    --cores all \
    --rerun-incomplete \
    --restart-times 2 \
    --executor slurm \
    -j5 
```

Some explanations:

* The convoluted arguments to `--config` pass the variables defined above to snakemake as wildcards. You could also put these in a config file, but this could get messy if you need to run the pipeline several times with different settings and aren't careful.
* The argument `-j N` tells the pipeline to run N jobs in parallel, assuming there are N input files. I used five as a place to start because the genotype calls are done by chromosome, and *Arabidopsis thaliana* has five chromosomes.
* `--rerun-incomplete` tells the pipeline to start any failed steps from scratch if you encounter and error and need to come back to it.
* `--restart-times 2` tells the pipeline to attempt the assembly or scaffolding steps with twice the memory if these fail initially.

Change these inputs as necessary.

## Output files

The pipeline outputs one dircetory per rule.
Some, but not all, rules also output stout and sterr to `logs/`.
It should be fairly clear what those files and logs are, but you can probably ignore them unless you need to debug something.

The main output file is a VCF file in the root folder of the output directory, which will be called the project name you passed to the pipeline, followed by `.vcf.gz`.
There should be an accompanying index file.
If the pipeline worked, these files are probably the only ones you need to keep.

You may also want to save the `multiqc` report which is run on fastq files after trimmin, found in the directory `trim_galore`.

## Acknowledgements

This was adapted from shell scripts by Pieter Clauw.
I used GPT 5.1 and Claude Sonnet for drafting rules.
I created the TOC for this README using https://bitdowntoc.derlin.ch/