#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SUBSET }      from './modules/subset.nf'
include { FILTER }      from './modules/filter.nf'
include { EXTRACT }     from './modules/extract.nf'
include { ATTACH }      from './modules/attach.nf'
include { DRAW }        from './modules/draw.nf'
include { CHECK }       from './modules/check.nf'
include { TEST }        from './modules/test.nf'

// Define input channels
pedigrees_file  = Channel.fromPath(params.pedigrees)
pedigrees_ch    = pedigrees_file
    | splitCsv(header: true, sep: '\t')
    | map { row -> [row.famid, row.id, row.fid, row.mid, row.sex, row.aff] }

variants_ch = Channel.fromFilePairs(params.variants, flat: true)
category_ch = Channel.of( 'Pathogenic', 'Damaging', 'Splicing' )
variable_ch = Channel.of( 'variants', 'annotations', 'genotypes', 'frequency' )
which_ch    = Channel.of('variant', 'gene')

workflow {
    // Identify: Subset, filter and extract variants
    pedigrees_ch 
        | groupTuple(by: 0)
        | map { it[0,1] }
        | combine(variants_ch)
        | SUBSET
        | combine(category_ch)
        | FILTER
        | combine(variable_ch)
        | EXTRACT
        | groupTuple(by: [0,1,3])
        | filter { it[3] == 'genotypes' }
        | set { genotypes_ch }

    // Family: Attach genotypes to pedigrees
    pedigrees_ch 
        | map { it.first() }
        | unique
        | combine(pedigrees_file)
        | combine(genotypes_ch, by: 0)
        | ATTACH
    
    // Sharing: Check and test sharing
    ATTACH.out | CHECK
    ATTACH.out | combine(which_ch) | TEST
    ATTACH.out | combine(TEST.out, by: [0,1]) | DRAW
}
