import os

# ---------- Helper function for read groups ----------
def get_read_group(sample):
    """
    Generate read group string for bwa mem -R flag
    Args:
        sample: The sample name (unique identifier)
    Returns:
        Read group string formatted for bwa mem
    """
    # Get the BioSample (individual) for this sample
    biosample = sample_dict[sample]['BioSample']
    
    # Construct read group string
    # ID and LB use the unique sample name
    # SM uses the BioSample (individual)
    rg_string = f"@RG\\tID:{sample}\\tSM:{biosample}\\tLB:{sample}\\tPL:ILLUMINA"
    
    return rg_string


rule align_reads:
    input:
        genome = fasta,
        r1 = "trimmed_reads/{sample}_val_1.fq.gz",
        r2 = "trimmed_reads/{sample}_val_2.fq.gz"
    output:
        bam="aligned_bams/{sample}.bam",
    params:
        rg = lambda wildcards: get_read_group(wildcards.sample)
    resources:
        qos='medium',
        mem_mb=lambda wildcards, attempt: 1024*10 * attempt,
        runtime= lambda wildcards, attempt: 60*8 * attempt,
    log:
        out = "logs/align_reads/{sample}.out",
        err = "logs/align_reads/{sample}.err"
    benchmark:
        "benchmarks/align_reads/{sample}.tsv"
    threads: 20
    shell:
        """
        bwa mem \
            -t {threads} \
            -R '{params.rg}' \
            {input.genome} \
            {input.r1} {input.r2} \
            -o {output.bam} \
        > {log.out} 2> {log.err}
        """