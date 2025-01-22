#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PEDIGREE }    from './modules/pedigree.nf'
include { SUBSET }      from './modules/subset.nf'
include { FILTER }      from './modules/filter.nf'
include { EXTRACT }     from './modules/extract.nf'
include { CONVERT }     from './modules/convert.nf'
include { REPORT }      from './modules/report.nf'
include { ATTACH }      from './modules/attach.nf'
include { SHARING }     from './modules/sharing.nf'
include { DRAW }        from './modules/draw.nf'

// Define input channels
families_ch = Channel.fromPath(params.pedigrees)
    | splitCsv(header: true, sep: '\t')
    | map { row -> [row.famid, row.id, row.fid, row.mid, row.sex, row.aff, row.famid] }
    | groupTuple(by: 0)

phenotypes_ch = Channel.fromPath(params.phenotypes)
variants_ch = Channel.fromFilePairs(params.variants, flat: true)
category_ch = Channel.of(params.categories.split(','))
type_ch     = Channel.of(params.type.split(','))

workflow {
    // Extract families, Subset and Filter
    families_ch
        | PEDIGREE
        | combine(variants_ch)
        | SUBSET
        | combine(category_ch)
        | FILTER
    
    // Reports
    SUBSET.out | REPORT

    // Extract variant annotations
    FILTER.out | EXTRACT

    // Extract genotypes, and indentify sharing
    FILTER.out
        | CONVERT 
        | ATTACH
        | combine(phenotypes_ch)
        | combine(EXTRACT.out, by: [0,1,2])
        | combine(type_ch)
        | SHARING
    
    // Export annotate pedigrees
    SHARING.out | DRAW

    // Export sharing tables
    SHARING.out
        | collectFile (
            keepHeader: true,
            storeDir: "${params.output_dir}/summary",
        )
        { it -> [ "${it[1]}.${it[3]}.tsv", it.last() ] }
}
