process FILTER {
    tag "${pheno}:${famid}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/filtered", mode: 'copy')

    input:
    tuple val(famid), val(pheno), 
          path(file), path(index),
          val(category)

    output:
    tuple val(famid), val(pheno), val(category),
          path("${famid}.${pheno}.${category}.vcf.gz"),
          path("${famid}.${pheno}.${category}.vcf.gz.tbi"),
          path("${famid}.${pheno}.${category}.tsv")
        
    script:
    if ( category == 'Pathogenic' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,gnomADe_AF:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'MAF > ${params.MAF} || HWE < ${params.HWE} || ExcHet < ${params.ExcHet}' | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'CLIN_SIG ~ "pathogenic" || CLIN_SIG ~ "likely_pathogenic"' | \
        bcftools view -e 'gnomADe_AF > ${params.PAF} || MAX_AF > ${params.PAF}' | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${famid}.${pheno}.${category}.vcf.gz

        tabix ${famid}.${pheno}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${famid}.${pheno}.${category}.vcf.gz > ${famid}.${pheno}.${category}.tsv
        """
    } else if ( category == 'Damaging' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,CADD_PHRED:Float,gnomADe_AF:Float,MAX_AF:Float ${file} | \
         bcftools view -e 'MAF > ${params.MAF} || HWE < ${params.HWE} || ExcHet < ${params.ExcHet}' | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'IMPACT="HIGH"' | \
        bcftools view -i 'CADD_PHRED > ${params.CADD}' | \
        bcftools view -e 'gnomADe_AF > ${params.gnomADe_AF} || MAX_AF > ${params.MAX_AF}' | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${famid}.${pheno}.${category}.vcf.gz

        tabix ${famid}.${pheno}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${famid}.${pheno}.${category}.vcf.gz > ${famid}.${pheno}.${category}.tsv
        """
    } else if ( category == 'Splicing' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,SpliceAI_pred_DS_AG:Float,SpliceAI_pred_DS_AL:Float,SpliceAI_pred_DS_DG:Float,SpliceAI_pred_DS_DL:Float,gnomADe_AF:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'MAF > ${params.MAF} || HWE < ${params.HWE} || ExcHet < ${params.ExcHet}' | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'SpliceAI_pred_DS_AG > ${params.DS} || SpliceAI_pred_DS_AL > ${params.DS} || SpliceAI_pred_DS_DG > ${params.DS} || SpliceAI_pred_DS_DL > ${params.DS}' | \
        bcftools view -e 'gnomADe_AF > ${params.gnomADe_AF} || MAX_AF > ${params.MAX_AF}' | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${famid}.${pheno}.${category}.vcf.gz

        tabix ${famid}.${pheno}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${famid}.${pheno}.${category}.vcf.gz > ${famid}.${pheno}.${category}.tsv
        """
    }
}
