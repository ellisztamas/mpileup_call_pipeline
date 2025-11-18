
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
        mem_mb=8*1024,
        runtime=240
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