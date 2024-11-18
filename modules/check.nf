process CHECK {
    tag "${pheno}:${famid}"

    label 'simple'

    container params.rvs

    publishDir("${params.output_dir}/family/${famid}", mode: 'copy')

    input:
    tuple val(famid), val(pheno), 
          path(genotype), path(phenotype), path(pedigree)

    output:
    tuple val(famid), val(pheno),
          path("${famid}.${pheno}.mendelian_errors.txt")
 
    script:
    """
    #!/bin/bash
    check.R ${famid} ${pheno} ${pedigree}
    """
}
