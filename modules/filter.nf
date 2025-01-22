process FILTER {
    tag "${pheno}:${famid}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/filtered", mode: 'copy')

    input:
    tuple val(famid), path(ped),
          val(pheno), path(file), path(index),
          val(category)

    output:
    tuple val(famid), path(ped),
          val(pheno), val(category),
          path("${famid}.${pheno}.${category}.vcf.gz"),
          path("${famid}.${pheno}.${category}.vcf.gz.tbi"),
          path("${famid}.${pheno}.${category}.tsv")
        
    script:
    if ( category == 'Pathogenic' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,${params.AF_COL}:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'CLIN_SIG ~ "pathogenic" || CLIN_SIG ~ "likely_pathogenic"' | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${famid}.${pheno}.${category}.vcf.gz

        tabix ${famid}.${pheno}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${famid}.${pheno}.${category}.vcf.gz > ${famid}.${pheno}.${category}.tsv
        """
    } else if ( category == 'Damaging' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,CADD_PHRED:Float,${params.AF_COL}:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'IMPACT="HIGH"' | \
        bcftools view -i 'CADD_PHRED > ${params.CADD}' | \
        bcftools view -e "${params.AF_COL} > ${params.AF} || MAX_AF > ${params.AF}" | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${famid}.${pheno}.${category}.vcf.gz

        tabix ${famid}.${pheno}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${famid}.${pheno}.${category}.vcf.gz > ${famid}.${pheno}.${category}.tsv
        """
    } else if ( category == 'High' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,CADD_PHRED:Float,${params.AF_COL}:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'IMPACT="HIGH"' | \
        bcftools view -e "${params.AF_COL} > ${params.AF} || MAX_AF > ${params.AF}" | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${famid}.${pheno}.${category}.vcf.gz

        tabix ${famid}.${pheno}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${famid}.${pheno}.${category}.vcf.gz > ${famid}.${pheno}.${category}.tsv
        """
    } else if ( category == 'Rare' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,CADD_PHRED:Float,${params.AF_COL}:Float,MAX_AF:Float ${file} | \
        bcftools view -e "${params.AF_COL} > ${params.AF} || MAX_AF > ${params.AF}" | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${famid}.${pheno}.${category}.vcf.gz

        tabix ${famid}.${pheno}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${famid}.${pheno}.${category}.vcf.gz > ${famid}.${pheno}.${category}.tsv
        """
    } else if ( category == 'Stop' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,Consequence,CADD_PHRED:Float,${params.AF_COL}:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'Consequence~"stop_gained"' | \
        bcftools view -e "${params.AF_COL} > ${params.AF} || MAX_AF > ${params.AF}" | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${famid}.${pheno}.${category}.vcf.gz

        tabix ${famid}.${pheno}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${famid}.${pheno}.${category}.vcf.gz > ${famid}.${pheno}.${category}.tsv
        """
    } else if ( category == 'PTV' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,Consequence,CADD_PHRED:Float,${params.AF_COL}:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'Consequence~"stop_gained" || Consequence~"frameshift_variant" || Consequence~"splice_acceptor_variant"' | \
        bcftools view -e "${params.AF_COL} > ${params.AF} || MAX_AF > ${params.AF}" | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${famid}.${pheno}.${category}.vcf.gz

        tabix ${famid}.${pheno}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${famid}.${pheno}.${category}.vcf.gz > ${famid}.${pheno}.${category}.tsv
        """
    } else if ( category == 'Scored' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,CADD_PHRED:Float,${params.AF_COL}:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'IMPACT="HIGH" || IMPACT="MODERATE"' | \
        bcftools view -i 'CADD_PHRED > ${params.CADD}' | \
        bcftools view -e "${params.AF_COL} > ${params.AF} || MAX_AF > ${params.AF}" | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${famid}.${pheno}.${category}.vcf.gz

        tabix ${famid}.${pheno}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${famid}.${pheno}.${category}.vcf.gz > ${famid}.${pheno}.${category}.tsv
        """
    } else if ( category == 'NotScored' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,IMPACT,CADD_PHRED:Float,${params.AF_COL}:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'IMPACT="HIGH" || IMPACT="MODERATE"' | \
        bcftools view -e "${params.AF_COL} > ${params.AF} || MAX_AF > ${params.AF}" | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${famid}.${pheno}.${category}.vcf.gz

        tabix ${famid}.${pheno}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${famid}.${pheno}.${category}.vcf.gz > ${famid}.${pheno}.${category}.tsv
        """
    } else if ( category == 'Splicing' ) {
        """
        #!/bin/bash
        # Filter cohort
        bcftools +split-vep -s worst -c CLIN_SIG,SpliceAI_pred_DS_AG:Float,SpliceAI_pred_DS_AL:Float,SpliceAI_pred_DS_DG:Float,SpliceAI_pred_DS_DL:Float,${params.AF_COL}:Float,MAX_AF:Float ${file} | \
        bcftools view -e 'CLIN_SIG ~ "conflicting" || CLIN_SIG ~ "benign"' | \
        bcftools view -i 'SpliceAI_pred_DS_AG > ${params.DS} || SpliceAI_pred_DS_AL > ${params.DS} || SpliceAI_pred_DS_DG > ${params.DS} || SpliceAI_pred_DS_DL > ${params.DS}' | \
        bcftools view -e "${params.AF_COL} > ${params.AF} || MAX_AF > ${params.AF}" | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view --threads ${task.cpu} -Oz -o ${famid}.${pheno}.${category}.vcf.gz

        tabix ${famid}.${pheno}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${famid}.${pheno}.${category}.vcf.gz > ${famid}.${pheno}.${category}.tsv
        """
    }
}
