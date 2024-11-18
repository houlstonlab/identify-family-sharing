process DRAW {
    tag "${pheno}:${famid}"

    label 'simple'

    container params.rvs

    publishDir("${params.output_dir}/sharing/${famid}", mode: 'copy')

    input:
    tuple val(famid), val(pheno), 
          path(genotype), path(phenotype), path(pedigree),
          val(which), path(markers)

    output:
    tuple val(famid), val(pheno), val(which),
          path("${famid}.${pheno}.${which}.*.png")
 
    script:
    """
    #!/bin/bash
    draw.R ${famid} ${pheno} ${genotype} ${phenotype} ${pedigree} ${which} ${markers}
    """
}
