process SHARING {
    tag "${famid}:${pheno}:${category}:${type}"

    label 'simple'

    container = params.rvs

    publishDir("${params.output_dir}/sharing", mode: 'copy')

    input:
    tuple val(famid), val(pheno), val(category),
          path(ped), path(phenotypes), path(variants),
          val(type)

    output:
    tuple val(famid), val(pheno), val(category), val(type),
          path(ped), path(phenotypes),
          path("${famid}.${pheno}.${category}.${type}.tsv")

    script:
    """
    #!/bin/bash
    sharing.R ${famid} ${pheno} ${category} ${type} ${ped} ${phenotypes} ${variants}
    """
}
