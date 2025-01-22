process CONVERT {
    tag "${famid}:${pheno}:${category}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/plinked", mode: 'copy')

    input:
    tuple val(famid), path(ped),
          val(pheno), val(category),
          path(vcf_in), path(index_in), path(variants)

    output:
    tuple val(famid), val(pheno), val(category),
          path("${famid}.${pheno}.${category}.ped"),
          path("${famid}.${pheno}.${category}.bim"),
          path("${famid}.${pheno}.${category}.bed"),
          path("${famid}.${pheno}.${category}.fam"),
          path("${famid}.${pheno}.${category}.nosex"),
          path("${famid}.${pheno}.${category}.log")

    script:
    """
    #!/bin/bash
    # Supbset pedigree columns
    cp ${ped} ${famid}.${pheno}.${category}.ped
    cat ${ped} | awk '{print \$6,\$1,\$2,\$3}' > parents.txt
    cat ${ped} | awk '{print \$6,\$1,\$4}'     > sex.txt
    cat ${ped} | awk '{print \$6,\$1,\$5}'     > aff.txt

    plink \
        --vcf ${vcf_in} \
        --make-bed \
        --const-fid ${famid} \
        --update-parents parents.txt \
        --update-sex sex.txt \
        --pheno aff.txt \
        --out ${famid}.${pheno}.${category}
    """
}
