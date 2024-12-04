process DRAW {
    tag "${pheno}:${famid}:${test}"

    label 'simple'

    container params.rvs

    publishDir("${params.output_dir}/plots", mode: 'copy')

    input:
    tuple val(famid), val(pheno),
          path(pedigree),
          val(test), path(test_file)

    output:
    tuple val(famid), val(pheno), val(test),
          path("${famid}.${pheno}.${test}.*.png")

    script:
    """
    #!/bin/bash
    draw.R ${pedigree} ${test_file} ${famid}.${pheno}.${test}
    """
}
