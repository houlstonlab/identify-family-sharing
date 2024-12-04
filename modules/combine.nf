process COMBINE {
    tag "${famid}:${pheno}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/combined/", mode: 'copy')

    input:
    tuple val(famid), val(pheno), val(category),
          path(file), path(index),
          path(variants)

    output:
    tuple val(famid), val(pheno),
          path("${famid}.${pheno}.vcf.gz"),
          path("${famid}.${pheno}.vcf.gz.tbi"),
          path("${famid}.${pheno}.tsv")
        
    script:
    """
    #!/bin/bash
    # Combine vcfs
    bcftools concat \
        -a -D \
        ${file} \
        --threads ${task.cpu} \
        -Oz -o ${famid}.${pheno}.vcf.gz
    tabix ${famid}.${pheno}.vcf.gz

    # Combine variants
    cat ${variants} | sort -u > ${famid}.${pheno}.tsv
    """
}

        