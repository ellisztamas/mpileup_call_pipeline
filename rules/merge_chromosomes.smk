rule merge_chromosomes:
    """
    Concatenate per-chromosome VCFs into a single cohort VCF.
    """
    input:
        vcfs = expand("mpileup_call_by_chrom/{chrom}.vcf.gz", chrom=chromosomes),
        indices = expand("mpileup_call_by_chrom/{chrom}.vcf.gz.csi", chrom=chromosomes)
    output:
        vcf = f"merge_chromosomes/{project_name}.vcf.gz"
    threads: 4
    resources:
        qos='medium',
        mem_mb = 4000,
        runtime = 8*60
    log:
        out="logs/merge_chromosomes.out",
        err="logs/merge_chromosomes.err"
    shell:
        """
        bcftools concat \
            --threads {threads} \
            --output-type z \
            --output {output.vcf} \
            {input.vcfs} \
            > {log.out} 2> {log.err}
        
        bcftools index --threads {threads} {output.vcf} 2>> {log}
        """