#! /usr/bin/env nextflow

include {
    FastQC as FastQCRaw;
    FastQC as FastQCClean;
} from './modules/fastqc.nf'
include {
    Cutadapt;
} from './modules/cutadapt.nf'
include {
    BWAmem;
    BWAIndex;
} from './modules/bwa.nf'
include{
    GenomeCov;
} from './modules/bedtools.nf'
include{
    CoveragePlotter;
} from './modules/coverplotter.nf'

workflow {

    Channel.fromFilePairs("raw_reads/*_R{1,2}*fastq.gz")
    | map{row -> [row[0].split("_S")[0], row[1]]}  
    | set { raw_reads }

    Channel.fromPath('samplesheet.csv')
    | splitCsv( header:true, sep:'\t' )
    | map { row -> [row.Sample, [sample:row.Sample, name:row.Name, reference:row.Reference]] }
    | join( raw_reads )
    | map { row -> [row[1], row[2]]}
    | set { samples }

    samples
    | branch {
        FindRef: it[0].reference == ''
            it[0].reference = "${it[0].sample}_ref"
            return [it[0].reference, it[1]]
        RefProvided: true
            def refFile = file(it[0].reference, checkIfExists: true)
            it[0].reference = refFile.simpleName
            return [ refFile.simpleName, refFile ]
    }
    | set { references }

    BWAIndex(references.RefProvided)

    FastQCRaw(samples, 'raw')
    Cutadapt(samples)
    FastQCClean(Cutadapt.out, 'clean')

    Cutadapt.out
    | map { row -> [row[0].reference, row[0], row[1]] }
    | combine( BWAIndex.out, by: 0)
    | map { row -> row.tail() }
    | BWAmem

    GenomeCov(BWAmem.out)
    | CoveragePlotter
}
