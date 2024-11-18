process ATTACH {
    tag "${pheno}:${famid}"

    label 'simple'

    container params.rvs

    publishDir("${params.output_dir}/family/${famid}", mode: 'copy')

    input:
    tuple val(famid), path(pedigree),
          val(pheno), val(category), val(variable),
          path(genotypes)

    output:
    tuple val(famid), val(pheno), 
          path("${famid}.${pheno}.genotype.tsv"),
          path("${famid}.${pheno}.phenotype.tsv"),
          path("${famid}.${pheno}.ped")

    script:
    """
    #!/bin/bash
    attach.R ${famid} ${pedigree} ${pheno} ${category.join(',')} ${genotypes.join(',')} 
    """
}
