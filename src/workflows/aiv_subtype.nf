#! /usr/bin/env nextflow

include {
    BWAMem ;
    BWAIndex
} from '../modules/bwa.nf'
include {
    GenomeCov
} from '../modules/bedtools.nf'
include {
    FindSubtypes
} from '../modules/python.nf'


workflow AIVSubtype {
    take:
    samples_ch

    main:
    Channel.fromFilePairs("${projectDir}/resources/aiv/*_subtypes.fa", size: 1)
        | set { references_ch }

    BWAIndex(references_ch)

    references_ch
        | join(BWAIndex.out)
        | set { indexed_references_ch }

    samples_ch
        | combine(indexed_references_ch)
        | map { [it[0], it[1], it[3], it[4]] }
        | BWAMem

    BWAMem.out
        | GenomeCov
        | FindSubtypes
        | map { [it[0], it[1]] }
        | groupTuple
        | set { result_ch }

    emit:
    result_ch
}
