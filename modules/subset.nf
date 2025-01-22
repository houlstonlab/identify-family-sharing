process SUBSET {
    tag "${pheno}:${famid}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/subsets", mode: 'copy')

    input:
    tuple val(famid), path(ped),
          val(pheno), path(file), path(index)

    output:
    tuple val(famid), path(ped),
          val(pheno), 
          path("${famid}.${pheno}.vcf.gz"),
          path("${famid}.${pheno}.vcf.gz.tbi")

    script:
    """
    #!/bin/bash
    # Subset pheno
    cat ${ped} | awk '{ print \$1 }' > samples.txt
    bcftools view --force-samples -g het -S samples.txt ${file} | \
    bcftools view -i 'FILTER="PASS"' | \
    bcftools norm -m -any | \
    bcftools +fill-tags -- -t all | \
    bcftools +setGT -- -t . -n 0 | \
    bcftools +setGT -- -t q -n 0 -i 'FMT/GQ < ${params.GQ} | FMT/DP < ${params.DP} | VAF < ${params.VAF}' | \
    bcftools +fill-tags -- -t all | \
    bcftools view -g het --threads ${task.cpu} -Oz -o ${famid}.${pheno}.vcf.gz

    tabix ${famid}.${pheno}.vcf.gz
    """
}
