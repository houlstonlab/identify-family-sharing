process ATTACH {
    tag "${pheno}:${famid}"

    label 'simple'

    container params.rvs

    publishDir("${params.output_dir}/markers/", mode: 'copy')

    input:
    tuple val(famid), val(pheno),
          path(bim), path(bed), path(fam), path(nosex), path(log),
          path(pedigree)

    output:
    tuple val(famid), val(pheno),
          path("${famid}.${pheno}.marked.ped")

    script:
    """
    #!/bin/bash
    attach.R ${bim} ${bed} ${fam} ${pedigree} ${famid}.${pheno}.marked
    """
}
