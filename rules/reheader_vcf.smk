rule reheader_vcf:
    input:
        vcf=f"merge_chromosomes/{project_name}.vcf.gz",
    output:
        f"{project_name}.vcf.gz",
    resources:
        qos='rapid',
        mem_mb=1024,
        runtime=10,
    log:
        "logs/reheader_vcf.log"
    shell:
        """
        # Extract current sample names from VCF
        bcftools query -l {input.vcf} > {input.vcf}.old_names.txt
        
        # Generate new sample names (basename without .bam)
        sed 's|.*/||; s|\.bam$||' {input.vcf}.old_names.txt > {input.vcf}.new_names.txt
        
        # Create reheader file (old_name new_name pairs)
        paste {input.vcf}.old_names.txt {input.vcf}.new_names.txt > {input.vcf}.rename.txt
        
        # Reheader the VCF
        bcftools reheader \
            --samples {input.vcf}.rename.txt \
            --output {output} \
            {input.vcf} \
        2> {log}
        
        # Index the new VCF
        bcftools index {output} 2>> {log}
        
        # Clean up temporary files
        rm {input.vcf}.old_names.txt {input.vcf}.new_names.txt {input.vcf}.rename.txt
        """