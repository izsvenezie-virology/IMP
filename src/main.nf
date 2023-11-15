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

    Channel.fromPath('references/*.fa')
    | map{ row -> [row.simpleName, row] }
    | set { references }

    BWAIndex(references)

    FastQCRaw(samples, 'raw')
    Cutadapt(samples)
    FastQCClean(Cutadapt.out, 'clean')

    Cutadapt.out
    | map { row -> [row[0].reference, row[0], row[1]] }
    | join( BWAIndex.out )
    | map { row -> row.tail() }
    | BWAmem

    GenomeCov(BWAmem.out)
}
