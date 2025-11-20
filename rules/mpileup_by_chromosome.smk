
MAX_MEM_MB = 64000

def scaled_mem(wildcards, attempt):
    mem = 8*1024 * (2 ** (attempt - 1))
    return min(mem, MAX_MEM_MB)

rule mpileup_call_by_chrom:
    """
    Joint genotyping of all samples at known sites, per chromosome.
    Uses bcftools mpileup | bcftools call pipeline.
    """
    input:
        bams = expand("deduplicated_bams/{sample}.bam", sample=SAMPLES),
        bais = expand("deduplicated_bams/{sample}.bam.bai", sample=SAMPLES),
        ref = config["fasta"],
        targets = config["targets"]
    output:
        vcf = temp("mpileup_call_by_chrom/{chrom}.vcf.gz"),
        indices = temp("mpileup_call_by_chrom/{chrom}.vcf.gz.csi")
    params:
        chrom = "{chrom}"
    threads: 4
    resources:
        qos='short',
        mem_mb=mem_for_attempt,
        runtime=60 * 4 * attempt
    log:
        out="logs/mpileup_call_by_chrom/{chrom}.out",
        err="logs/mpileup_call_by_chrom/{chrom}.err"
    shell:
        """
        bcftools mpileup \
            --threads {threads} \
            --fasta-ref {input.ref} \
            --regions {params.chrom} \
            --targets-file {input.targets} \
            --annotate FORMAT/AD,FORMAT/DP \
            --max-depth 1000 \
            --min-MQ 20 \
            --min-BQ 20 \
            --output-type u \
            {input.bams} \
        | bcftools call \
            --threads {threads} \
            --multiallelic-caller \
            --variants-only \
            --targets-file {input.targets} \
            --output-type z \
            --output {output.vcf} \
            > {log.out} 2> {log.err}
        
        bcftools index --threads {threads} {output.vcf} 2>> {log}
        """