"""
Collate fastqc results into a single report.

Note that the name for the output file is hard-coded as 
trimmed_reads/multiqc_report_post_trimming.html.
This will not automatically be updated if you run the pipeline again.
"""
rule multiqc:
    input:
        "trimmed_reads"
    output:
        f"trimmed_reads/{config['project_name']}_multiqc_report_post_trimming.html"
    resources:
        qos='rapid',
        mem_mb=4*1024,
        runtime=60,
    benchmark:
        "benchmarks/multiqc.tsv"
    shell:
        "multiqc --filename {output} {input}"