process EXTRACT {
    tag "${pheno}:${famid}:${category}:${variable}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/extracted", mode: 'copy')

    input:
    tuple val(famid), val(pheno), val(category),
          path(file), path(index), path(variants),
          val(variable)

    output:
    tuple val(famid), val(pheno), val(category), val(variable),
	      path("${famid}.${pheno}.${category}.${variable}.tsv")

    script:
    if ( variable == 'variants' ) {
		"""
		#!/bin/bash
		bcftools +split-vep \
			-s worst \
			-c SYMBOL \
			-f '%SYMBOL\t%CHROM:%POS:%REF:%ALT\n' \
			${file} \
			> "${famid}.${pheno}.${category}.${variable}.tsv"
		"""
    } else if ( variable == 'annotations' ) {
		"""
		#!/bin/bash
		bcftools +split-vep \
			-s worst \
			-c SYMBOL \
			-f '%SYMBOL\t%CHROM:%POS:%REF:%ALT\t%CSQ\n' \
			-d -A tab \
			${file} \
			> "${famid}.${pheno}.${category}.${variable}.tsv"
		"""
    } else if ( variable == 'genotypes' ) {
		"""
		#!/bin/bash
		bcftools +split-vep \
			-s worst \
			-c SYMBOL \
			-f '%SYMBOL\t%CHROM:%POS:%REF:%ALT[\t%SAMPLE=%GT]\n' \
			${file} \
			> "${famid}.${pheno}.${category}.${variable}.tsv"
		"""
    } else if ( variable == 'frequency' ) {
		"""
		#!/bin/bash
		bcftools +split-vep \
			-s worst \
			-c SYMBOL \
			-f '%SYMBOL\t%CHROM:%POS:%REF:%ALT[\t%SAMPLE=%GT]\n' \
			${file} | \
			summarize_genotypes.awk \
			> "${famid}.${pheno}.${category}.${variable}.tsv"
		"""
    } 
}
