import pandas as pd

# Name for the VCF file.
project_name=config.get("project_name")

# Get sample sheet path from CLI/config; can provide a default if you like
SAMPLESHEET = config.get("sample_sheet")
samples_df = pd.read_csv(SAMPLESHEET)
# Check whether any samples appear more than once.
counts = samples_df['sample'].value_counts()
if any(counts>1):
  raise ValueError("One or more samples have identical names. Either merge the fastq files, or give them unique identifiers.")
# Dict of samples
sample_dict = samples_df.set_index("sample", drop=False).to_dict(orient="index")
SAMPLES = samples_df["sample"].tolist()

# Path to the reference genome
fasta = config.get("fasta")
# List of chromosome names.
chromosomes = [line.split()[0] for line in open(f"{fasta}.fai")]

# Path to file with known SNP positions
targets=config.get("targets")


# ---------- Include rule modules ----------
include: "rules/trim_galore.smk"
include: "rules/multiqc.smk"
include: "rules/align_reads.smk"
include: "rules/remove_duplicate_reads.smk"
include: "rules/generate_gvcf.smk"
include: "rules/merge_gvcfs.smk"
include: "rules/call_genotypes.smk"
include: "rules/merge_chromosomes.smk"
include: "rules/reheader_vcf.smk"



# ---------- Final output rule ----------
rule all:
    input:
        f"{project_name}.vcf.gz",
        f"trimmed_reads/{project_name}_multiqc_report_post_trimming.html"