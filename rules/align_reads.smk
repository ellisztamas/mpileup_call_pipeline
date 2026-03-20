import os

def calculate_runtime(wildcards, attempt, input):
    """Calculate runtime based on input file sizes with 30-minute minimum"""
    # Get file sizes in bytes and convert to GB
    r1_size_gb = os.path.getsize(input.r1) / (1024**3)
    r2_size_gb = os.path.getsize(input.r2) / (1024**3)
    
    # Calculate runtime: sum of sizes * 7 hours
    calculated_runtime = (r1_size_gb + r2_size_gb) * 7 * 60
    
    # Apply minimum
    runtime_minutes = max(calculated_runtime, 30)
    
    return int(runtime_minutes * attempt)

rule align_and_deduplicate:
    input:
        genome = fasta,
        r1 = "trimmed_reads/{sample}_val_1.fq.gz",
        r2 = "trimmed_reads/{sample}_val_2.fq.gz"
    output:
        bam = "deduplicated_bams/{sample}.bam",
        bai = "deduplicated_bams/{sample}.bam.bai"
    resources:
        qos = 'medium',
        mem_mb = lambda wildcards, attempt: 1024*24 * attempt,
        runtime = calculate_runtime
    log:
        err = "logs/align_reads/{sample}.err"
    benchmark:
        "benchmarks/align_reads/{sample}.tsv"
    threads: 12
    shell:
        r"""
        set -o pipefail

        minimap2 -ax sr -t {threads} \
            -R "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}" \
            {input.genome} {input.r1} {input.r2} 2> {log.err} | \
        samblaster -r 2>> {log.err} | \
        samtools sort -@ {threads} -m 2G -o {output.bam} - 2>> {log.err}

        samtools index -@ {threads} {output.bam} 2>> {log.err}
        """