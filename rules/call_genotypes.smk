rule call_genotypes:
    input:
        merged_gvcf = "merged_gvcfs/{chrom}.merged.g.vcf.gz",
        merged_tbi  = "merged_gvcfs/{chrom}.merged.g.vcf.gz.tbi"
    output:
        vcf = "call_genotypes/{chrom}.vcf.gz",
        tbi = "call_genotypes/{chrom}.vcf.gz.tbi"
    threads: 2
    resources:
        qos     = 'short',
        mem_mb  = lambda wildcards, attempt: 1024 * 4 * attempt,
        runtime = lambda wildcards, attempt: 60 * attempt
    log:
        out = "logs/call_genotypes/{chrom}.out",
        err = "logs/call_genotypes/{chrom}.err"
    benchmark:
        "benchmarks/call_genotypes/{chrom}.tsv"
    shell:
        """
        bcftools call \
            --threads {threads} \
            --ploidy 2 \
            --multiallelic-caller \
            --keep-alts \
            --format-fields GQ \
            --output-type z \
            --output {output.vcf} \
            {input.merged_gvcf} \
        > {log.out} 2> {log.err} && \
        bcftools index --tbi {output.vcf}
        """