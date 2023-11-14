#! /usr/bin/env nextflow

include {
    FastQC as FastQCRaw;
    FastQC as FastQCClean;
} from './modules/fastqc.nf'
include {
    Cutadapt;
} from './modules/cutadapt.nf'

workflow {

    Channel.fromFilePairs("raw_reads/*_R{1,2}*fastq.gz")
    | map{row -> [row[0].split("_S")[0], row[1]]}  
    | set { raw_reads }

    Channel.fromPath('samplesheet.csv')
    | splitCsv( header:true, sep:'\t' )
    | map { row -> [row.Sample, [sample:row.Sample, name:row.Name]] }
    | join( raw_reads )
    | map { row -> [row[1], row[2]]}
    | set { samples }

    FastQCRaw(samples, 'raw')
    Cutadapt(samples)
    FastQCClean(Cutadapt.out, 'clean')
}
