process REPORT {
    tag "${pheno}:${famid}:${report}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/reprots", mode: 'copy')

    input:
    tuple val(famid), val(pheno),
          path(bim), path(bed), path(fam), path(nosex), path(log)
    each report

    output:
    tuple val(famid), val(pheno),
          path("${famid}.${pheno}.*")

    script:
    """
    #!/bin/bash        
    plink \
        --bfile ${bim.baseName} \
        --${report} \
        --out ${famid}.${pheno}
    """
}
