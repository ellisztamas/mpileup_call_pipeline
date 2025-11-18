import pandas as pd

# Project name to give to the multiqc and VCF files.
project_name = config.get("project_name", "genotype_calls")
        
# Path to the file listing variable sites.
variable_sites = config.get("variable_sites")

# Get sample sheet path from CLI/config; can provide a default if you like
SAMPLESHEET = config.get("sample_sheet")
samples_df = pd.read_csv(SAMPLESHEET)
sample_dict = samples_df.set_index("sample", drop=False).to_dict(orient="index")
SAMPLES = samples_df["sample"].tolist()

# List of chromosome names.
chromosomes = [line.split()[0] for line in open(f"{config["fasta"]}.fai")]


# ---------- Include rule modules ----------
include: "rules/trim_galore.smk"
include: "rules/multiqc.smk"
include: "rules/align_reads.smk"
include: "rules/remove_duplicate_reads.smk"
include: "rules/mpileup_by_chromosome.smk"
include: "rules/merge_chromosomes.smk"
include: "rules/reheader_vcf.smk"



# ---------- Final output rule ----------
rule all:
    input:
        f"{project_name}.vcf.gz",
        f"trimmed_reads/{project_name}_multiqc_report_post_trimming.html"