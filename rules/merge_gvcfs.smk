rule merge_gvcfs:
    input:
        gvcfs = expand("gvcfs/{sample}.{{chrom}}.g.vcf.gz", sample=SAMPLES),
        tbis = expand("gvcfs/{sample}.{{chrom}}.g.vcf.gz.tbi", sample=SAMPLES)
    output:
        merged_gvcf = "merged_gvcfs/{chrom}.merged.g.vcf.gz",
        merged_tbi = "merged_gvcfs/{chrom}.merged.g.vcf.gz.tbi"
    threads: 4 # Merging many files benefits from a few more threads for I/O
    resources:
        mem_mb = 8*1024,
        runtime = 60
    log:
        out = "logs/merge_gvcfs/{chrom}.out",
        err = "logs/merge_gvcfs/{chrom}.err"
    benchmark:
        "benchmarks/merge_gvcfs/{chrom}.tsv"
    shell:
        """
        (
            bcftools merge \
                --threads {threads} \
                --output-type z \
                {input.gvcfs} > {output.merged_gvcf}
            
            bcftools index --tbi {output.merged_gvcf}
        ) > {log.out} 2> {log.err}
        """