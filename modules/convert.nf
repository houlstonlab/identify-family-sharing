process CONVERT {
    tag "${pheno}:${famid}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/plinked", mode: 'copy')

    input:
    tuple val(famid), val(pheno),
          path(vcf_in), path(index_in), path(variants),
          path(parents), path(sex), path(aff)

    output:
    tuple val(famid), val(pheno),
          path("${famid}.${pheno}.bim"),
          path("${famid}.${pheno}.bed"),
          path("${famid}.${pheno}.fam"),
          path("${famid}.${pheno}.nosex"),
          path("${famid}.${pheno}.log")

    script:
    """
    #!/bin/bash
    plink \
        --vcf ${vcf_in} \
        --make-bed \
        --const-fid ${famid} \
        --update-parents ${parents} \
        --update-sex ${sex} \
        --pheno ${aff} \
        --out ${famid}.${pheno}
    """
}
