rule merge_chromosomes:
    """
    Concatenate per-chromosome VCFs into a single cohort VCF.
    """
    input:
        vcf = expand("call_genotypes/{chrom}.vcf", chrom=chromosomes),
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
    benchmark:
        "benchmarks/merge_chromosomes.tsv"
    shell:
        """
        bcftools concat \
            --threads {threads} \
            --output-type z \
            --output {output.vcf} \
            {input.vcf} \
            > {log.out} 2> {log.err}
        
        bcftools index --threads {threads} {output.vcf} 2>> {log}
        """