rule summarise_vcf:
    """
    Create a summary of the output VCF file using bcftools stats.
    """
    input:
        f"{project_name}.vcf.gz"
    output:
        stats=f"summarise_vcf/{project_name}.vchk",
        pdf="summarise_vcf/summary.pdf"
    shell:
        """
        bcftools stats {input} > {output.stats}
        plot-vcfstats -p summarise_vcf {output.stats}
        """