# Set the numbers of bases to clip from 5` end of reads 
clip_R1 = config.get("clip_R1", 15)  # default 15 if not set
clip_R2 = config.get("clip_R2", 15)  # default 15 if not set
quality = config.get("quality", 20)
length = config.get("length", 20)

rule trim_galore:
    input:
        r1 = lambda wildcards: sample_dict[wildcards.sample]["fastq1"],
        r2 = lambda wildcards: sample_dict[wildcards.sample]["fastq2"]
    output:
        r1 = "trimmed_reads/{sample}_val_1.fq.gz",
        r2 = "trimmed_reads/{sample}_val_2.fq.gz"
    params:
        outdir = "trimmed_reads"
    resources:
        qos='short',
        mem_mb=lambda wildcards, attempt: 1024 * (2**(1-attempt)),
        runtime=lambda wildcards, attempt: 60 * attempt,
    log:
        out = "logs/trim_galore/{sample}.out",
        err = "logs/trim_galore/{sample}.err"
    threads: 10
    shell:
        """
        trim_galore \
            --paired \
            --quality {quality} \
            --length {length} \
            --clip_R1 {clip_R1} --clip_R2 {clip_R2} \
            --basename {wildcards.sample} \
            --fastqc \
            --trim-n \
            --cores {threads} \
            --output_dir {params.outdir} \
            {input.r1} {input.r2} \
            > {log.out} 2> {log.err}
        """
