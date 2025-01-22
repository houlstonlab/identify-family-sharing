process ATTACH {
    tag "${famid}:${pheno}:${category}"

    label 'simple'

    container params.rvs

    publishDir("${params.output_dir}/markers/", mode: 'copy')

    input:
    tuple val(famid), val(pheno), val(category), path(ped),
          path(bim), path(bed), path(fam), path(nosex), path(log)

    output:
    tuple val(famid), val(pheno), val(category),
          path("${famid}.${pheno}.${category}.marked.ped")

    script:
    """
    #!/bin/bash
    attach.R ${bim} ${bed} ${fam} ${ped} ${famid}.${pheno}.${category}.marked
    """
}
