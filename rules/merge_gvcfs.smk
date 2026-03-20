import math

BATCH_SIZE = 100  # ~100 files open per job

# Create batch indices for each chromosome
def get_batches(samples, batch_size):
    n = math.ceil(len(samples) / batch_size)
    return [str(i) for i in range(n)]

BATCHES = get_batches(SAMPLES, BATCH_SIZE)

def get_batch_samples(wildcards):
    """Return the samples for a given batch index."""
    i = int(wildcards.batch)
    return SAMPLES[i * BATCH_SIZE : (i + 1) * BATCH_SIZE]


# --- Step 1: Merge each batch ---
rule merge_gvcfs_batch:
    input:
        gvcfs = lambda wildcards: expand(
            "gvcfs/{sample}.{chrom}.g.vcf.gz",
            sample=get_batch_samples(wildcards),
            chrom=wildcards.chrom
        ),
        tbis = lambda wildcards: expand(
            "gvcfs/{sample}.{chrom}.g.vcf.gz.tbi",
            sample=get_batch_samples(wildcards),
            chrom=wildcards.chrom
        ),
        fasta = config["fasta"]
    output:
        merged_gvcf = temp("merged_gvcfs/batches/{chrom}.batch{batch}.g.vcf.gz"),
        merged_tbi  = temp("merged_gvcfs/batches/{chrom}.batch{batch}.g.vcf.gz.tbi")
    threads: 4
    resources:
        open_files = BATCH_SIZE * 2 + 10,
        qos        = 'rapid',
        mem_mb     = lambda wildcards, attempt: 1024 * 4 * attempt,
        runtime    = lambda wildcards, attempt: 60 * attempt
    log:
        out = "logs/merge_gvcfs_batch/{chrom}.batch{batch}.out",
        err = "logs/merge_gvcfs_batch/{chrom}.batch{batch}.err"
    benchmark:
        "benchmarks/merge_gvcfs_batch/{chrom}.batch{batch}.tsv"
    shell:
        """
        bcftools merge \
            --threads {threads} \
            --gvcf {input.fasta} \
            --merge none \
            --output-type z \
            --output {output.merged_gvcf} \
            {input.gvcfs} \
        > {log.out} 2> {log.err} && \
        bcftools index --tbi {output.merged_gvcf}
        """


# --- Step 2: Merge all batch GVCFs into the final output ---
rule merge_gvcfs_final:
    input:
        gvcfs = expand(
            "merged_gvcfs/batches/{{chrom}}.batch{batch}.g.vcf.gz",
            batch=BATCHES
        ),
        tbis = expand(
            "merged_gvcfs/batches/{{chrom}}.batch{batch}.g.vcf.gz.tbi",
            batch=BATCHES
        ),
        fasta = config["fasta"]
    output:
        merged_gvcf = "merged_gvcfs/{chrom}.merged.g.vcf.gz",
        merged_tbi  = "merged_gvcfs/{chrom}.merged.g.vcf.gz.tbi"
    threads: 4
    resources:
        open_files = len(BATCHES) * 2 + 10,
        qos        = 'rapid',
        mem_mb     = lambda wildcards, attempt: 1024 * 4 * attempt,
        runtime    = lambda wildcards, attempt: 60 * attempt
    log:
        out = "logs/merge_gvcfs_final/{chrom}.out",
        err = "logs/merge_gvcfs_final/{chrom}.err"
    benchmark:
        "benchmarks/merge_gvcfs_final/{chrom}.tsv"
    shell:
        """
        if [ $(echo {input.gvcfs} | wc -w) -eq 1 ]; then
            cp {input.gvcfs} {output.merged_gvcf}
            cp {input.tbis} {output.merged_tbi}
        else
            bcftools merge \
                --threads {threads} \
                --gvcf {input.fasta} \
                --merge none \
                --output-type z \
                --output {output.merged_gvcf} \
                {input.gvcfs} \
            > {log.out} 2> {log.err} && \
            bcftools index --tbi {output.merged_gvcf}
        fi
        """