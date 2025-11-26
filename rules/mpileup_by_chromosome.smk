
rule mpileup_call_by_chrom:
    """
    Joint genotyping of all samples at known sites, per chromosome.
    Uses bcftools mpileup | bcftools call pipeline.
    """
    input:
        bam = expand("deduplicated_bams/{sample}.bam", sample=SAMPLES),
        bai = expand("deduplicated_bams/{sample}.bam.bai", sample=SAMPLES),
    output:
        vcf = "mpileup_call_by_chrom/{chrom}.vcf.gz",
        indices = "mpileup_call_by_chrom/{chrom}.vcf.gz.csi"
    params:
        chrom = "{chrom}"
    threads: 16
    resources:
        qos='short',
        mem_mb  = lambda wildcards, attempt: 8*1024 * (2**(attempt-1)),
        runtime = lambda wildcards, attempt: 60 * 4 * attempt
    log:
        out="logs/mpileup_call_by_chrom/{chrom}.out",
        err="logs/mpileup_call_by_chrom/{chrom}.err"
    shell:
        """
        bcftools mpileup \
            --threads {threads} \
            --fasta-ref {fasta} \
            --regions {params.chrom} \
            --targets-file {targets} \
            --annotate FORMAT/AD,FORMAT/DP \
            --max-depth 1000 \
            --min-MQ 20 \
            --min-BQ 20 \
            --output-type u \
            {input.bam} \
        | bcftools call \
            --threads {threads} \
            --multiallelic-caller \
            --variants-only \
            --targets-file {targets} \
            --output-type z \
            --output {output.vcf} \
            > {log.out} 2> {log.err}
        
        bcftools index --threads {threads} {output.vcf} 2>> {log}
        """