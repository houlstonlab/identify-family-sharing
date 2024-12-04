process TEST {
    tag "${pheno}:${famid}:${test}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/tests", mode: 'copy')

    input:
    tuple val(famid), val(pheno),
          path(bim), path(bed), path(fam), path(nosex), path(log)
    each test

    output:
    tuple val(famid), val(pheno), val(test),
          path("${famid}.${pheno}.${test}")

    script:
    """
    #!/bin/bash        
    plink \
        --bfile ${bim.baseName} \
        --${test} \
        --out ${famid}.${pheno}
    # TODO: implement othe tests
    """
}
