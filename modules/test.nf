process TEST {
    tag "${pheno}:${famid}"

    label 'simple'

    container params.rvs

    publishDir("${params.output_dir}/sharing/${famid}", mode: 'copy')

    input:
    tuple val(famid), val(pheno), 
          path(genotype), path(phenotype), path(pedigree),
          val(which)

    output:
    tuple val(famid), val(pheno), val(which),
          path("${famid}.${pheno}.${which}.sharing.tsv")
 
    script:
    """
    #!/bin/bash
    test.R ${famid} ${pheno} ${genotype} ${phenotype} ${pedigree} ${which}
    """
}
