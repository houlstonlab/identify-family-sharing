#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PEDIGREE }    from './modules/pedigree.nf'
include { SUBSET }      from './modules/subset.nf'
include { FILTER }      from './modules/filter.nf'
include { COMBINE }     from './modules/combine.nf'
include { CONVERT }     from './modules/convert.nf'
include { REPORT }      from './modules/report.nf'
include { TEST }        from './modules/test.nf'
include { ATTACH }      from './modules/attach.nf'
include { DRAW }        from './modules/draw.nf'

// Define input channels
families_ch = Channel.fromPath(params.pedigrees)
    | splitCsv(header: true, sep: '\t')
    | map { row -> [row.famid, row.id, row.fid, row.mid, row.sex, row.aff] }
    | groupTuple(by: 0)

variants_ch = Channel.fromFilePairs(params.variants, flat: true)

category_ch = Channel.of(params.categories.split(','))
report_ch = Channel.of(params.reports.split(','))
test_ch = Channel.of(params.tests.split(','))

workflow {
    // Load Pedigrees
    families_ch
        | PEDIGREE
        | multiMap {
            samples : it[0,1]
            info    : it[0,2,3,4]
            ped     : it[0,5]
        }
        | set {pedigree }

    // Subset, Filter, Combine and Convert
    pedigree.samples
        | combine(variants_ch)
        | SUBSET
        | combine(category_ch)
        | FILTER
        | groupTuple(by: [0,1])
        | COMBINE
        | combine(pedigree.info, by: 0)
        | CONVERT

    // Reports
    REPORT(CONVERT.out, report_ch)

    // Tests
    TEST(CONVERT.out, test_ch)

    // Attach genotypes and draw pedigree
    CONVERT.out
        | combine(pedigree.ped, by: 0)
        | ATTACH
        | combine(TEST.out, by: [0,1])
        | DRAW
}
