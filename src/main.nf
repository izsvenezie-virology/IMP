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
    RemoveDegenerations;
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
    raw_reads
    | join      ( metadata_ch )                                                     // join metadata to keep only samples to analyze
    | map       { [it[0], it[1]] }                                                  // keep only ID and raw reads files
    | set       { to_fastqcraw_ch }                                                 // set channel

    FastQCRaw   ( to_fastqcraw_ch, 'raw' )
    FastQCClean ( Cutadapt.out, 'clean' )

    // References collection channel creation
    metadata_ch
    | combine   ( Cutadapt.out, by: 0 )
    | branch    {
                  FindRef: it[1].reference_file == ''
                      return [it[1].reference, it[1].subset ?: 1, it[2]]
                  RefProvided: true
                      return [it[1].reference, file(it[1].reference_file, checkIfExists: true)]
    }
    | set       { ref_collect }

    // Find missing references
    MakeBlastDb (references_db, ref_collect.FindRef.first().ifEmpty(false))         // create blastDB only if at least one reference must be found
    FastqToFasta( ref_collect.FindRef )
    BlastN      ( FastqToFasta.out, MakeBlastDb.out )
    GetReferenceNames( BlastN.out )
    GetReference( GetReferenceNames.out, references_db )

    // References channel creation
    GetReference.out
    | filter    { !it[1].isEmpty() }                                                // remove empty references: in this case the sample using this reference is not processed
    | mix       ( ref_collect.RefProvided )                                         // merge provided references and found refereces
    | unique                                                                        // keeps only one copy for each reference
    | RemoveDegenerations
    
    RemoveDegenerations.out
    | set       { References }                                                      // set channel

    // Reference index processes
    BWAIndex    ( References )
    FaidxIndex  ( References )
    DictIndex   ( References )

    // Reads alignment
    metadata_ch
    | map       { [ it[1].sample, it[1].reference, it[1].id]}                       // .sample and .reference for channel merge, id is a parameter
    | combine   (Cutadapt.out, by: 0)                                               // combine cleaned reads
    | map       { it.tail() }                                                       // remove .sample
    | combine   ( References, by: 0)                                                // combine reference
    | combine   ( BWAIndex.out, by: 0 )                                             // combine reference index
    | map       { it.tail() }                                                       // remove .reference
    | BWAmem

    BamIndex    (BWAmem.out)

    GenomeCov   ( BWAmem.out )
    | CoveragePlotter

    // GATK best practices
    FixBam      ( BWAmem.out )

    CleanBam    ( FixBam.out )

    CleanBam.out
    | map       { [it[0].reference] + it }                                          // .refrence for channel merge
    | combine   ( References, by: 0)                                                // combine reference
    | combine   ( BWAIndex.out, by: 0 )                                             // combine reference index
    | map       { it.tail() }                                                       // remover .refrence
    | Viterbi

    MDSort      (Viterbi.out)

    MarkDuplicates( MDSort.out )

    MDBamIndex  ( MarkDuplicates.out )

    MarkDuplicates.out
    | combine   ( MDBamIndex.out, by: 0 )                                           // combine bam index
    | map       { [it[0].reference] + it }                                          // .reference for channel merge
    | combine   ( References, by: 0 )                                               // combine reference
    | combine   ( FaidxIndex.out, by: 0 )                                           // combine reference index
    | map       { it.tail() }                                                       // remove .reference
    | set       { md_reference }                                                    // set channel
    FakeVariantCall(md_reference, false)

    IndexFeatureFile(FakeVariantCall.out)

    MarkDuplicates.out
    | map       { [it[0].reference] + it }                                          // .reference for channel merge
    | combine   ( References, by: 0)                                                // combine reference
    | combine   ( FaidxIndex.out, by: 0 )                                           // combine reference index
    | combine   ( DictIndex.out, by: 0)                                             // combine reference dictionary
    | map       { it.tail() }                                                       // remove .reference
    | combine   ( FakeVariantCall.out, by: 0 )                                      // combine vcf
    | combine   ( IndexFeatureFile.out, by: 0 )                                     // combine vcf index
    | BaseRecalibrator

    MarkDuplicates.out
    | combine   ( BaseRecalibrator.out, by: 0 )                                     // combine recalibration table
    | ApplyBQSR

    ApplyBQSR.out
    | map       { [it[0].reference] + it }                                          // .reference for channel merge
    | combine   ( References, by: 0 )                                               // combine reference
    | combine   ( FaidxIndex.out, by: 0 )                                           // combine reference index
    | map       { it.tail() }                                                       // remove .reference
    | set       { bqsr_reference }                                                  // set channel
    VariantCall ( bqsr_reference, true )

    // Create consensuses
    metadata_ch
    | map       { [it[1].id, it[1].reference, it[1].id, [name: it[1].name]] }       // .id and .reference for channel merge, .id and .name are parameters
    | combine   ( VariantCall.out, by: 0 )                                          // combine vcf file
    | combine   ( GenomeCov.out, by: 0 )                                            // combine coverage file
    | map       { it.tail() }                                                       // remove .id
    | combine   ( References, by: 0 )                                               // combine reference file
    | map       { it.tail() }                                                       // remove .reference
    | set       { vcf_coverage_reference }                                          // set channel
    DegeneratedConsensus( vcf_coverage_reference )
    NonDegeneratedConsensus( vcf_coverage_reference )

    // Creates files with all sequences splitted by segment
    DegeneratedConsensus.out.segments
    | map       { it[1] }                                                           // Get only files
    | flatten                                                                       // Transform to list
    | collectFile( storeDir: 'results' )                                            // Merge content of files by name
}
