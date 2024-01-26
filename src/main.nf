#! /usr/bin/env nextflow

include{
    GetReference;
} from './modules/bash.nf'
include{
    GenomeCov;
} from './modules/bedtools.nf'
include{
    BlastN
} from './modules/blast.nf'
include {
    BWAmem;
    BWAIndex;
} from './modules/bwa.nf'
include{
    Consensus as ConsensusDegenerated;
    Consensus as ConsensusNotDegenerated;
} from './modules/consenser.nf'
include{
    CoveragePlotter;
} from './modules/coverplotter.nf'
include {
    Cutadapt;
} from './modules/cutadapt.nf'
include {
    FastQC as FastQCRaw;
    FastQC as FastQCClean;
} from './modules/fastqc.nf'
include{
    DictIndex;
    FixBam;
    CleanBam;
    MarkDuplicates;
    BaseRecalibrator;
    IndexFeatureFile;
    ApplyBQSR;
} from './modules/gatk.nf'
include{
    Viterbi;
    Call as FakeVariantCall;
    Call as VariantCall;
} from './modules/lofreq.nf'
include{
    GetReferenceNames;
} from './modules/python.nf'
include{
    FaidxIndex;
    Sort as MDSort;
    BamIndex as MDBamIndex;
    BamIndex as BamIndex;
} from './modules/samtools.nf'
include{
    FastqToFasta
} from './modules/seqtk.nf'

workflow {
    // BLAST DB channels
    Channel.fromPath(params.references_database)
    | toList
    | set { references_db }
    Channel.fromPath([params.references_database, "${params.references_database}.*"])
    | toList
    | set { references_db_blast }

    // Samples channels creation
    Channel.fromFilePairs("${params.raw_reads_folder}/*_R{1,2}*fastq.gz")
    | map { row -> [row[0].split("_S")[0], row[1]] }  
    | set { raw_reads }

    Channel.fromPath(params.samples_metadata)
    | splitCsv( header:true, sep:'\t' )
    | map { row -> 
            reference = row.Reference ? file(row.Reference).simpleName : "${row.Sample}_ref"
            [row.Sample, [sample:row.Sample, name:row.Name, reference:reference, reference_file:row.Reference]] }
    | join( raw_reads )
    | map { row -> [row[1], row[2]] }
    | set { samples }

    // Clean reads
    FastQCRaw( samples, 'raw' )
    Cutadapt( samples )
    FastQCClean( Cutadapt.out, 'clean' )

    // References collection channel creation
    Cutadapt.out
    | branch {
        FindRef: it[0].reference_file == ''
            return [it[0].reference, it[1]]
        RefProvided: it[0].reference_file != ''
            return [it[0].reference, file(it[0].reference_file, checkIfExists: true)]
    }
    | set { ref_collect }

    // Find missing references
    FastqToFasta( ref_collect.FindRef )
    BlastN( FastqToFasta.out, references_db_blast )
    GetReferenceNames( BlastN.out )
    GetReference( GetReferenceNames.out, references_db )

    // References channel creation
    GetReference.out
    | filter { !it[1].isEmpty() }
    | mix ( ref_collect.RefProvided )
    | unique
    | set { References }

    // Reference index processes
    BWAIndex( References )
    FaidxIndex( References )
    DictIndex( References )

    // Reads alignment
    Cutadapt.out
    | map { row -> [row[0].reference] + row }
    | combine( References, by: 0)
    | combine( BWAIndex.out, by: 0 )
    | map { row -> row.tail() }
    | BWAmem

    BamIndex(BWAmem.out)

    GenomeCov( BWAmem.out )
    | CoveragePlotter

    // GATK best practices
    FixBam( BWAmem.out )

    CleanBam( FixBam.out )

    CleanBam.out
    | map { row -> [row[0].reference] + row }
    | combine( References, by: 0)
    | combine( BWAIndex.out, by: 0 )
    | map { row -> row.tail() }
    | Viterbi

    MDSort(Viterbi.out)

    MarkDuplicates( MDSort.out )

    MDBamIndex( MarkDuplicates.out )

    MarkDuplicates.out
    | combine( MDBamIndex.out, by: 0 )
    | map { row -> [row[0].reference] + row }
    | combine( References, by: 0 )
    | combine( FaidxIndex.out, by: 0 )
    | map { row -> row.tail() }
    | set { md_reference }
    FakeVariantCall(md_reference, false)

    IndexFeatureFile(FakeVariantCall.out)

    MarkDuplicates.out
    | map { row -> [row[0].reference] + row }
    | combine( References, by: 0)
    | combine( FaidxIndex.out, by: 0 )
    | combine( DictIndex.out, by: 0)
    | map { row -> row.tail() }
    | combine( FakeVariantCall.out, by: 0 )
    | combine( IndexFeatureFile.out, by: 0 )
    | BaseRecalibrator

    MarkDuplicates.out
    | combine( BaseRecalibrator.out, by: 0 )
    | ApplyBQSR

    ApplyBQSR.out
    | map { row -> [row[0].reference] + row }
    | combine( References, by: 0 )
    | combine( FaidxIndex.out, by: 0 )
    | map { row -> row.tail() }
    | set { bqsr_reference }
    VariantCall( bqsr_reference, true )

    VariantCall.out
    | map { row -> [row[0].reference] + row }
    | combine( References, by: 0 )
    | map { row -> row.tail() }
    | combine( GenomeCov.out, by: 0 )
    | set { vcf_reference_coverage } 
    ConsensusDegenerated( vcf_reference_coverage, true )
    ConsensusNotDegenerated( vcf_reference_coverage, false )

    ConsensusDegenerated.out.segments
    | map { row-> row[1] }
    | flatten
    | collectFile( storeDir: 'results' )
}
