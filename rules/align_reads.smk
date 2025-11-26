rule align_reads:
    input:
        genome = fasta,
        r1 = "trimmed_reads/{sample}_val_1.fq.gz",
        r2 = "trimmed_reads/{sample}_val_2.fq.gz"
    output:
        bam="aligned_bams/{sample}.bam",
        # bai="aligned_bams/{sample}.bam.bai"
    resources:
        qos='short',
        mem_mb=lambda wildcards, attempt: 1024*10 * attempt,
        runtime= lambda wildcards, attempt: 60 * attempt,
    log:
        out = "logs/align_reads/{sample}.out",
        err = "logs/align_reads/{sample}.err"
    threads: 10
    shell:
        """
        bwa mem \
            -o {output.bam} \
            {input.genome} \
            -t {threads} \
            {input.r1} {input.r2} \
        > {log.out} 2> {log.err}

        """