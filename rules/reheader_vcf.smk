rule reheader_vcf:
    input:
        f"merge_chromosomes/{project_name}.vcf.gz",
    output:
        f"{project_name}.vcf.gz",
    resources:
        qos='rapid',
        mem_mb=1024,
        runtime=10,
    shell:
        """
        # Extract current sample names from VCF
        bcftools query -l {input} > {input}.old_names.txt
        
        # Generate new sample names (basename without .bam)
        sed 's|.*/||; s|\.bam$||' {input}.old_names.txt > {input}.new_names.txt
        
        # Create reheader file (old_name new_name pairs)
        paste {input}.old_names.txt {input}.new_names.txt > {input}.rename.txt
        
        # Reheader the VCF
        bcftools reheader \
            --samples {input}.rename.txt \
            --output {output} \
            {input}
        
        # Index the new VCF
        bcftools index {output}
        
        # Clean up temporary files
        rm {input}.old_names.txt {input}.new_names.txt {input}.rename.txt
        """