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
        qos='rapid',
        mem_mb=1024,
        runtime=60,
    log:
        out = "logs/trim_galore/{sample}.out",
        err = "logs/trim_galore/{sample}.err"
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
            --output_dir {params.outdir} \
            {input.r1} {input.r2} \
            > {log.out} 2> {log.err}
        """
