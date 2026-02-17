"""
Identify and mark reads that are likely PCR duplicates.

Steps 

- collate: group reads by name
- fixmate Add ms and MC tags for markdup to use later.
    - Flag -m adds mate score tag
- sort: Sort reads by position
- markdup: Identify and remove (-r) duplicate reads.
- index: Index the output .bam file.
"""

rule remove_duplicate_reads:
    input:
        bam="aligned_bams/{sample}.bam",
    output:
        collated=temp("deduplicated_bams/{sample}_collated.bam"),
        fixmate=temp("deduplicated_bams/{sample}_fixmate.bam"),
        sorted=temp("deduplicated_bams/{sample}_sorted.bam"),
        markdup="deduplicated_bams/{sample}.bam",
        bai="deduplicated_bams/{sample}.bam.bai"
    resources:
        qos='short',
        mem_mb =  lambda wildcards, attempt: 10*1024 * (2**(attempt-1)),
        runtime = lambda wildcards, attempt: 60*2*attempt,
    threads:10
    log:
        out="logs/remove_duplicate_reads/{sample}.out",
        err="logs/remove_duplicate_reads/{sample}.err"
    benchmark:
        "benchmarks/remove_duplicate_reads/{sample}.tsv"
    shell:
        """
        (
            samtools collate -@ {threads} -o {output.collated} {input.bam}
            samtools fixmate -@ {threads} -m {output.collated} {output.fixmate}
            samtools sort -@ {threads} -o {output.sorted} {output.fixmate}
            samtools markdup -@ {threads} -r {output.sorted} {output.markdup}
            samtools index {output.markdup}
        ) > {log.out} 2> {log.err}
        """