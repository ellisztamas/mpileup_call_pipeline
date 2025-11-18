rule align_reads:
    input:
        genome = lambda wildcards: sample_dict[wildcards.sample]["genome"],
        r1 = "trimmed_reads/{sample}_val_1.fq.gz",
        r2 = "trimmed_reads/{sample}_val_2.fq.gz"
    output:
        bam="aligned_bams/{sample}.bam",
        bai="aligned_bams/{sample}.bam.bai"
    resources:
        qos='short',
        mem_mb=1024,
        runtime=2*60,
    log:
        out = "logs/align_reads/{sample}.out",
        err = "logs/align_reads/{sample}.err"
    shell:
        """
        bwa mem \
            -o {output.bam} \
            {input.genome} \
            {input.r1} {input.r2} \
            > {log.out} 2> {log.err}
        samtools index {output.bam}
        """