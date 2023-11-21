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
include{
    FastqToFasta
} from './modules/seqtk.nf'
include{
    BlastN
} from './modules/blast.nf'
include{
    GetReferenceNames;
} from './modules/python.nf'
include{
    GetReference;
} from './modules/bash.nf'

workflow {
    // BLAST DB channels
    Channel.fromPath('ref_db/gisaid_epiflu_all_unique.fa.*')
    | toList
    | set { ref_db }
    Channel.fromPath('ref_db/gisaid_epiflu_all_unique.fa')
    | toList
    | set { ref_fasta }

    // Samples channels creation
    Channel.fromFilePairs("raw_reads/*_R{1,2}*fastq.gz")
    | map { row -> [row[0].split("_S")[0], row[1]] }  
    | set { raw_reads }

    Channel.fromPath('samplesheet.csv')
    | splitCsv( header:true, sep:'\t' )
    | map { row -> [row.Sample, [sample:row.Sample, name:row.Name, reference_file:row.Reference]] }
    | join( raw_reads )
    | map { row -> [row[1], row[2]] }
    | set { samples }

    // Clean reads
    FastQCRaw( samples, 'raw' )
    Cutadapt( samples )
    FastQCClean( Cutadapt.out, 'clean' )

    // References collection channel creation
    samples
    | branch {
        FindRef: it[0].reference_file == ''
            it[0].reference = "${it[0].sample}_ref"
            return [it[0].reference, it[1]]
        RefProvided: it[0].reference_file != ''
            def refFile = file(it[0].reference_file, checkIfExists: true)
            it[0].reference = refFile.simpleName
            return [it[0].reference, refFile]
    }
    | set { ref_collect }

    // Find missing references
    FastqToFasta( ref_collect.FindRef )
    BlastN( FastqToFasta.out, ref_db )
    GetReferenceNames( BlastN.out )
    GetReference( GetReferenceNames.out, ref_fasta )

    // References channel creation
    GetReference.out
    | mix ( ref_collect.RefProvided )
    | unique
    | set { references }

    // Reference index processes
    BWAIndex( references )

    // Reads alignment
    Cutadapt.out
    | map { row -> [row[0].reference, row[0], row[1]] }
    | combine( BWAIndex.out, by: 0 )
    | map { row -> row.tail() }
    | BWAmem

    GenomeCov( BWAmem.out )
    | CoveragePlotter
}
