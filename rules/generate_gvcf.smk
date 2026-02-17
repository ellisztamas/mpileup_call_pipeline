rule generate_gvcf:
    input:
        fasta = config["fasta"],
        bam = "deduplicated_bams/{sample}.bam",
        bai = "deduplicated_bams/{sample}.bam.bai",
        targets = config["targets"]
    output:
        gvcf = "gvcfs/{sample}.{chrom}.g.vcf.gz",
        tbi = "gvcfs/{sample}.{chrom}.g.vcf.gz.tbi"
    params:
        chrom = "{chrom}",
        max_depth = config.get("max_depth", 10000)
        min_MQ = config.get("min_MQ", 15)
        min_BQ = config.get("min_BQ", 20)
    threads: 2
    resources:
        mem_mb = 4*1024,
        runtime = 120 # 2 hours is usually plenty for a single sample/chrom
    log:
        out = "logs/generate_gvcf/{sample}.{chrom}.out",
        err = "logs/generate_gvcf/{sample}.{chrom}.err"
    benchmark:
        "benchmarks/generate_gvcf/{sample}.{chrom}.tsv"
    shell:
        """
        (
            bcftools mpileup \
                --threads {threads} \
                --fasta-ref {input.fasta} \
                --regions {params.chrom} \
                --targets-file {input.targets} \
                --annotate FORMAT/AD,FORMAT/DP \
                --skip-indels \
                --adjust-MQ 50 \
                --max-depth {params.max_depth} \
                --min-MQ {params.min_MQ} \
                --min-BQ {params.min_BQ} \
                --gvcf 0 \
                --output-type z \
                {input.bam} > {output.gvcf}
            
            bcftools index --tbi {output.gvcf}
        ) > {log.out} 2> {log.err}
        """