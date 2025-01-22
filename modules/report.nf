process REPORT {
    tag "${famid}:${pheno}"

    label 'simple'

    container = params.plink

    publishDir("${params.output_dir}/reports", mode: 'copy')

    input:
    tuple val(famid), path(ped), val(pheno), 
          path(vcf_in), path(vcf_index)
          
    output:
    tuple val(famid), val(pheno), 
          path("${famid}.${pheno}.ped"),
          path("${famid}.${pheno}.bim"),
          path("${famid}.${pheno}.bed"),
          path("${famid}.${pheno}.fam"),
          path("${famid}.${pheno}.nosex"),
          path("${famid}.${pheno}.log"),
          path("${famid}.${pheno}.*")

    script:
    """
    #!/bin/bash 
    # Supbset pedigree columns
    cp ${ped} ${famid}.${pheno}.ped

    cat ${ped} | awk '{print \$6,\$1,\$2,\$3}' > parents.txt
    cat ${ped} | awk '{print \$6,\$1,\$4}'     > sex.txt
    cat ${ped} | awk '{print \$6,\$1,\$5}'     > aff.txt

    # Convert to plink
    plink \
        --vcf ${vcf_in} \
        --make-bed \
        --const-fid ${famid} \
        --update-parents parents.txt \
        --update-sex sex.txt \
        --pheno aff.txt \
        --out ${famid}.${pheno}
    
    # Run mendel test
    plink \
        --bfile ${famid}.${pheno} \
        --mendel-multigen \
        --nonfounders \
        --out tmp
        
    cat tmp.mendel  | sed 's/  */\t/g' | sed 's/^[ \t]*//' > ${famid}.${pheno}.mendel
    cat tmp.imendel | sed 's/  */\t/g' | sed 's/^[ \t]*//' > ${famid}.${pheno}.imendel
    cat tmp.lmendel | sed 's/  */\t/g' | sed 's/^[ \t]*//' > ${famid}.${pheno}.lmendel
    cat tmp.fmendel | sed 's/  */\t/g' | sed 's/^[ \t]*//' > ${famid}.${pheno}.fmendel
    
    # Run kinship test
    if [ \$(wc -l < ${famid}.${pheno}.fam) -gt 1 ]; then
        plink \
            --bfile ${famid}.${pheno} \
            --make-rel square \
            --nonfounders \
            --out tmp

        cat tmp.rel > ${famid}.${pheno}.rel
        cat tmp.rel.id > ${famid}.${pheno}.rel.id
    else
        touch ${famid}.${pheno}.rel
        touch ${famid}.${pheno}.rel.id
    fi
    """
}
