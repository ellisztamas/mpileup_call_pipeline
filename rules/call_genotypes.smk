rule call_genotypes:
    input:
        merged_gvcf = "merged_gvcfs/{chrom}.merged.g.vcf.gz",
        merged_tbi = "merged_gvcfs/{chrom}.merged.g.vcf.gz.tbi"
    output:
        vcf = "call_genotypes/{chrom}.vcf"
    threads: 2
    resources:
        mem_mb = 4096,
        runtime = 60
    log:
        out = "logs/call_genotypes/{chrom}.out",
        err = "logs/call_genotypes/{chrom}.err"
    benchmark:
        "benchmarks/call_genotypes/{chrom}.tsv"
    shell:
        """
        (
            bcftools call \
                --threads {threads} \
                --ploidy 2 \
                --multiallelic-caller \
                --variants-only \
                --format-fields GQ \
                --output-type v \
                {input.merged_gvcf} \
            | sed 's/##INFO=<ID=MQ,Number=1,Type=Integer/##INFO=<ID=MQ,Number=1,Type=Float/' \
            > {output.vcf}
        ) > {log.out} 2> {log.err}
        """