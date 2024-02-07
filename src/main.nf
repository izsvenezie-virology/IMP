#! /usr/bin/env nextflow

include{
    GetReference;
} from './modules/bash.nf'
include{
    GenomeCov;
} from './modules/bedtools.nf'
include{
    BlastN
    MakeBlastDb
} from './modules/blast.nf'
include {
    BWAmem;
    BWAIndex;
} from './modules/bwa.nf'
include{
    DegeneratedConsensus;
    NonDegeneratedConsensus;
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
    CreateCutadaptPrimers;
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
    // BLAST DB channel
    references_db = file(params.references_database, checkIfExists: true)

    // Adapters channel
    adapters = file(params.adapters, checkIfExists: true)

    // Samples channels creation
    Channel.fromFilePairs("${params.raw_reads_folder}/*_R{1,2}*fastq.gz")
    | map { [it[0].split("_S")[0], it[1]] }  
    | set { raw_reads }

    Channel.fromPath(params.samples_metadata)
    | splitCsv( header:true, sep:'\t' )
    | map { 
            reference = it.Reference ? file(it.Reference).simpleName : "${it.Sample}_ref"
            primers = it.Primers ? file(it.Primers).simpleName : file(params.null_file).simpleName
            [it.Sample, [
                id:[sample: it.Sample, reference: reference],
                sample:it.Sample,
                name:it.Name,
                primers:primers,
                primers_file:it.Primers,
                reference:reference,
                reference_file:it.Reference,
                subset:it.Subset
            ]] }
    | set { metadata_ch }

    metadata_ch
    | map { 
        if (it[1].primers_file)
            [it[1].primers, file(it[1].primers_file, checkIfExists: true)] 
        else
            [it[1].primers, file(params.null_file, checkIfExists: true)]}
    | unique
    | CreateCutadaptPrimers

    // Clean reads
    raw_reads
    | join      ( metadata_ch )
    | map       { [it[2].primers, it[2].sample, [it[2].minQual, it[2].minLen], it[1]] }
    | combine   ( CreateCutadaptPrimers.out, by: 0 )
    | map       { it.tail() }
    | set       { to_cutadapt_ch }
    Cutadapt    ( to_cutadapt_ch, adapters )

    // Assess reads quality
    FastQCRaw   ( raw_reads, 'raw' )
    FastQCClean ( Cutadapt.out, 'clean' )

    // References collection channel creation
    metadata_ch
    | combine ( Cutadapt.out, by: 0 )
    | branch {
        FindRef: it[1].reference_file == ''
            return [it[1].reference, it[1].subset ?: 1, it[1]]
        RefProvided: true
            return [it[1].reference, file(it[1].reference_file, checkIfExists: true)]
    }
    | set { ref_collect }

    // Find missing references
    MakeBlastDb(references_db, ref_collect.FindRef.first().ifEmpty(false))
    FastqToFasta( ref_collect.FindRef )
    BlastN( FastqToFasta.out, MakeBlastDb.out )
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
    metadata_ch
    | map { [ it[1].sample, it[1].reference, it[1].id]}
    | combine (Cutadapt.out, by: 0)
    | map { it.tail() }
    | combine( References, by: 0)
    | combine( BWAIndex.out, by: 0 )
    | map { it.tail() }
    | BWAmem

    BamIndex(BWAmem.out)

    GenomeCov( BWAmem.out )
    | CoveragePlotter

    // GATK best practices
    FixBam( BWAmem.out )

    CleanBam( FixBam.out )

    CleanBam.out
    | map { [it[0].reference] + it }
    | combine( References, by: 0)
    | combine( BWAIndex.out, by: 0 )
    | map { it.tail() }
    | Viterbi

    MDSort(Viterbi.out)

    MarkDuplicates( MDSort.out )

    MDBamIndex( MarkDuplicates.out )

    MarkDuplicates.out
    | combine( MDBamIndex.out, by: 0 )
    | map { [it[0].reference] + it }
    | combine( References, by: 0 )
    | combine( FaidxIndex.out, by: 0 )
    | map { it.tail() }
    | set { md_reference }
    FakeVariantCall(md_reference, false)

    IndexFeatureFile(FakeVariantCall.out)

    MarkDuplicates.out
    | map { [it[0].reference] + it }
    | combine( References, by: 0)
    | combine( FaidxIndex.out, by: 0 )
    | combine( DictIndex.out, by: 0)
    | map { it.tail() }
    | combine( FakeVariantCall.out, by: 0 )
    | combine( IndexFeatureFile.out, by: 0 )
    | BaseRecalibrator

    MarkDuplicates.out
    | combine( BaseRecalibrator.out, by: 0 )
    | ApplyBQSR

    ApplyBQSR.out
    | map { [it[0].reference] + it }
    | combine( References, by: 0 )
    | combine( FaidxIndex.out, by: 0 )
    | map { it.tail() }
    | set { bqsr_reference }
    VariantCall( bqsr_reference, true )

    metadata_ch
    | map { [it[1].id, it[1].reference, it[1].id, [name: it[1].name]] }
    | combine (VariantCall.out, by: 0)
    | combine ( GenomeCov.out, by: 0 )
    | map { it.tail() }
    | combine ( References, by: 0 )
    | map { it.tail() }
    | set { vcf_coverage_reference } 
    DegeneratedConsensus( vcf_coverage_reference )
    NonDegeneratedConsensus( vcf_coverage_reference )

    DegeneratedConsensus.out.segments
    | map { it[1] }
    | flatten
    | collectFile( storeDir: 'results' )
}
