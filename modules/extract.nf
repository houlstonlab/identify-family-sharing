process EXTRACT {
    tag "${pheno}:${famid}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/extracted", mode: 'copy')

    input:
    tuple val(famid), path(ped), val(pheno), val(category),
          path(file), path(index), path(variants)

    output:
    tuple val(famid), val(pheno), val(category),
          path("${famid}.${pheno}.${category}.annotation.tsv")

    script:
    """
    #!/bin/bash
    # Create header 
    echo -e "SNP\t\$(bcftools +split-vep -l ${file} | cut -f 2 | tr '\n' '\t')" > "${famid}.${pheno}.${category}.annotation.tsv"

    # Extract variant annotations
    bcftools +split-vep \
	    -s worst \
		-f '%CHROM:%POS:%REF:%ALT\t%CSQ\n' \
		-d -A tab \
		${file} \
		>> "${famid}.${pheno}.${category}.annotation.tsv"
    """
}
