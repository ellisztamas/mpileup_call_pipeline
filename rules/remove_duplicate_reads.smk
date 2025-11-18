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
        bai="aligned_bams/{sample}.bam.bai"
    output:
        bam="deduplicated_bams/{sample}.bam",
        bai="deduplicated_bams/{sample}.bam.bai"
    resources:
        qos='short',
        mem_mb=4*1024,
        runtime=2*60,
    shell:
        """
        samtools collate     -O -u {input.bam} | \
            samtools fixmate -m -u - - | \
            samtools sort    -u - | \
            samtools markdup -r - {output.bam}
        samtools index {output.bai}
        """