#! /usr/bin/env nextflow

nextflow.preview.output = true


include {
    GetReference ;
    ConcatenateConensus
} from './modules/bash.nf'
include {
    GenomeCov
} from './modules/bedtools.nf'
include {
    BlastN ;
    MakeBlastDb
} from './modules/blast.nf'
include {
    BWAMem ;
    BWAIndex
} from './modules/bwa.nf'
include {
    Consenser
} from './modules/consenser.nf'
include {
    Cutadapt
} from './modules/cutadapt.nf'
include {
    FastQC as FastQCRaw ;
    FastQC as FastQCClean
} from './modules/fastqc.nf'
include {
    DictIndex ;
    FixBam ;
    CleanBam ;
    MarkDuplicates ;
    BaseRecalibrator ;
    IndexFeatureFile ;
    ApplyBQSR
} from './modules/gatk.nf'
include {
    Viterbi ;
    Call as FakeVariantCall ;
    Call as VariantCall ;
    IndelQual
} from './modules/lofreq.nf'
include {
    CreateCutadaptPrimers ;
    GetReferenceNames ;
    PrepareReference
} from './modules/python.nf'
include {
    FaidxIndex ;
    Sort as MDSort ;
    BamIndex as IQBamIndex ;
    BamIndex as MDBamIndex ;
    BamIndex as BamIndex
} from './modules/samtools.nf'
include {
    FastqToFasta
} from './modules/seqtk.nf'
include {
    Tacos
} from './modules/tacos.nf'
include {
    ReadsStats as ReadsStatsRaw ;
    ReadsStats as ReadsStatsClean ;
    CoverageStats ;
    AlignmentStats ;
    VariantsStats
} from './modules/statistics/python.nf'

include {
    UpdateFluMut ;
    FluMut
} from './modules/aiv/flumut.nf'
include {
    UpdateGenin2 ;
    Genin2
} from './modules/aiv/genin2.nf'
include {
    AIVGetSubtype
} from './modules/aiv/python.nf'

workflow {
    main:
    log.info(
        """ 
                   .*%%%%%%%%%#.
                   =***********-
                   -***********=:.
                +%%%***********#%%%:
                =#%%%%%%%%%%%%%%%#=
                 *%%#####**#####%%:
                -%%%%%%%%%%%%%%%%%*.
               *%%%%%%%%%%%%%%%%%%%%+
             +%%%%%%%%%%%%%%%%%%%@@@@%+
           +%%%%%%%%%%%%%%%%%%%%%@@@@@%#-
          *%%%%%%%#%%%%%%%%%%%%%%%%@@@@%%-
         =%%%%%%%%  %%%%%%%%*:%%%%%%%@@@%*.
         =%%%%%%%%   *%%%%%%* .%%%%%%%@@%%=
         +%%%%%%%*             .#%%%%%@@%%=
         =%%%%%%+                 -%%%%@%%=
         =%%%%#                     +%%@@%=
         =%%%%                       =%%%%=
         =%%%-         .-+++-.        %%%%=
         =%%%.      %@@@@@@@#.=@+     #%%%=
         =%%%     :@@@@@@@@:    @@    *%%%=
         =%%%     #@@@@@@@@-    @@:   *%%%=
         =%%%     .@@@@@@@@@%-*@@#    +%%%=
         =%%%       +@@@@@@@@@@@      +%%%=
         =%%%           .::.          +%%%=
         =%%%                         +%%%=
         =%%%                         +%%%=
         =%%%                         +%%%=
         =%%%                         +%%%=
         =%%%                         +%%%=
         =%%%    ___   __  __   ___   +%%%=
         =%%%   |_ _| |  \\/  | | _ \\  +%%%=
         =%%%    | |  | |\\/| | |  _/  +%%%=
         =%%%   |___| |_|  |_| |_|    +%%%=
         =%%%                         +%%%=
         =%%%                         +%%%=
         +%%%                         +%%%=
         =%%%%.                      *%%%%=
          *%%%%%+                 .#%%%%%*
           :#%%%%%%#+:       .-+#%%%%%#=
              =#%%%%%%%%%%%%%%%%%%%%+.
                  =**%%%%%%%%%%*+.
    """
    )

    // Adapters channel
    adapters = file(params.adapters, checkIfExists: true)

    // Samples channels creation
    Channel.fromFilePairs("${params.raw_reads_folder}/${params.raw_reads_pattern}")
        | map { it -> [it[0].split("_S")[0], it.tail().flatten()] }
        | set { raw_reads }

    Channel.fromPath(params.samples_metadata)
        | splitCsv(header: true, sep: '\t')
        | multiMap { it ->
            // If a reference is not provided set the file to null, create a custom id and set the subset > 0.
            // If a reference is provided and the file exists use it and set the.baseName as ID. The subset must be null.
            // Else if the file does not exists set the file to null and use the.baseName as ID. The subset must be null.
            def reference_file = file(it.Reference ?: 'non_existing_file').exists() ? file(it.Reference) : null
            def reference_id = it.Reference ? file(it.Reference).baseName : "${it.Sample}_ref"

            def subset = it.Reference ? null : it.Subset ?: 1

            def primers_file = file(it.Primers ?: params.null_file, checkIfExists: true)
            def primers_id = primers_file.baseName
            metadata: [
                id: "${it.Sample}__${reference_id}",
                sample: it.Sample,
                name: it.Name,
                primers: primers_id,
                reference: reference_id,
                subset: subset,
                group: it.Group ?: params.run_name,
                minimum_coverage: it.MinimumCoverage ?: 20,
                no_filters: it.NoConsensusFilter ? true : false,
            ]
            primers: [primers_id, primers_file]
            references: [reference_id, reference_file]
        }
        | set { samples_config_ch }

    // Metadata channel creation
    samples_config_ch.metadata
        | unique
        | filter { it -> !it.sample.startsWith("#") }
        | set { metadata_ch }

    samples_config_ch.primers
        | unique
        | CreateCutadaptPrimers

    // Check if all samples has fastq files
    metadata_ch
        | map { meta -> [meta.sample] }
        | unique
        | join(raw_reads, remainder: true)
        | filter { it -> it[1] == null }
        | map { error("<${it[0]}> is missing fastq files") }

    // Clean reads
    metadata_ch
        | map { meta -> [meta.sample, meta.primers, meta.sample, [phred_threshold: meta.minQual, min_len: meta.minLen]] }
        | combine(raw_reads, by: 0)
        | map { it -> it.tail() }
        | combine(CreateCutadaptPrimers.out, by: 0)
        | map { it -> it.tail() }
        | unique
        | set { to_cutadapt_ch }
    Cutadapt(to_cutadapt_ch, adapters)

    metadata_ch
        | map { meta -> [meta.sample] }
        | combine(raw_reads, by: 0)
        | unique
        | set { to_fastqc_raw_ch }

    // Reads quality assesment
    FastQCRaw(to_fastqc_raw_ch)
    FastQCClean(Cutadapt.out)
    ReadsStatsRaw(to_fastqc_raw_ch, 'raw')
    ReadsStatsClean(Cutadapt.out, 'clean')


    // References collection channel creation
    samples_config_ch.references
        | filter { it -> it[1] }
        | set { provided_references }

    // Find missing references
    metadata_ch
        | filter { it -> it.subset }
        | first
        | map { file(params.references_database, checkIfExists: true) }
        | set { references_db }

    // create blastDB only if at least one reference must be found
    MakeBlastDb(references_db)
    metadata_ch
        | filter { it -> it.subset }
        | map { meta -> [meta.sample, meta.reference, meta.subset] }
        | combine(Cutadapt.out, by: 0)
        | map { it -> it.tail() }
        | FastqToFasta
    BlastN(FastqToFasta.out, MakeBlastDb.out)
    GetReferenceNames(BlastN.out)
    GetReference(GetReferenceNames.out, references_db)

    // Alert if some references are empty
    GetReference.out
        | filter { it -> it[1].isEmpty() }
        | view { it -> "<${it[0]}> is an empty reference" }
        | collectFile(
            storeDir: "${workflow.outputDir}/warnings"
        ) { it -> ['empty_references.txt', "${it[0]}\n"] }

    // Check if all the samples has a reference (even if it's empty)
    metadata_ch
        | map { meta -> [meta.reference, meta.id] }
        | unique
        | join(GetReference.out.mix(provided_references), remainder: true)
        | map { it -> it.tail() }
        | filter { it -> it[1] == null }
        | map { error("<${it[0]}> reference file does not exist") }

    // References channel creation
    GetReference.out
        | filter { it -> !it[1].isEmpty() }
        | mix(provided_references)
        | unique
        | PrepareReference
        | set { References }

    // Reference index processes
    BWAIndex(References)
    FaidxIndex(References)
    DictIndex(References)

    // Reads alignment
    metadata_ch
        | map { meta -> [meta.sample, meta.reference, meta.id] }
        | combine(Cutadapt.out, by: 0)
        | map { it -> it.tail() }
        | combine(References, by: 0)
        | combine(BWAIndex.out, by: 0)
        | map { it -> it.tail() }
        | BWAMem
        | BamIndex

    GenomeCov(BWAMem.out)
        | combine(metadata_ch.map { meta -> [meta.id, meta.minimum_coverage] }, by: 0)
        | Tacos

    BWAMem.out
        | AlignmentStats

    AlignmentStats.out
        | filter { it -> it[1].text =~ '\tTotal\tMapped reads\t0' }
        | view { it -> "<${it[0]}> has 0 mapped reads" }
        | collectFile(storeDir: "${workflow.outputDir}/warnings") { it -> ['no_mapped_reads.txt', "${it[0]}\n"] }

    metadata_ch
        | map { meta -> [meta.id, meta.minimum_coverage] }
        | combine(GenomeCov.out, by: 0)
        | CoverageStats

    // GATK best practices
    AlignmentStats.out
        | filter { it -> !(it[1].text =~ '\tTotal\tMapped reads\t0') }
        | map { it -> [it[0]] }
        | join(BWAMem.out)
        | FixBam
        | CleanBam

    metadata_ch
        | map { meta -> [meta.id, meta.reference, meta.id] }
        | combine(CleanBam.out, by: 0)
        | map { it -> it.tail() }
        | combine(References, by: 0)
        | combine(BWAIndex.out, by: 0)
        | map { it -> it.tail() }
        | Viterbi
        | MDSort
        | MarkDuplicates
        | MDBamIndex

    metadata_ch
        | map { meta -> [meta.id, meta.reference, meta.id, []] }
        | combine(MarkDuplicates.out, by: 0)
        | combine(MDBamIndex.out, by: 0)
        | map { it -> it.tail() }
        | combine(References, by: 0)
        | combine(FaidxIndex.out, by: 0)
        | map { it -> it.tail() }
        | FakeVariantCall
        | IndexFeatureFile

    metadata_ch
        | map { meta -> [meta.id, meta.reference, meta.id] }
        | combine(MarkDuplicates.out, by: 0)
        | map { it -> it.tail() }
        | combine(References, by: 0)
        | combine(FaidxIndex.out, by: 0)
        | combine(DictIndex.out, by: 0)
        | map { it -> it.tail() }
        | combine(FakeVariantCall.out, by: 0)
        | combine(IndexFeatureFile.out, by: 0)
        | BaseRecalibrator

    MarkDuplicates.out
        | combine(BaseRecalibrator.out, by: 0)
        | ApplyBQSR

    metadata_ch
        | map { meta -> [meta.id, meta.reference, meta.id] }
        | combine(ApplyBQSR.out, by: 0)
        | map { it -> it.tail() }
        | combine(References, by: 0)
        | map { it -> it.tail() }
        | IndelQual
        | IQBamIndex

    metadata_ch
        | map { meta -> [meta.id, meta.reference, meta.id, [no_filters: meta.no_filters]] }
        | combine(IndelQual.out, by: 0)
        | combine(IQBamIndex.out, by: 0)
        | map { it -> it.tail() }
        | combine(References, by: 0)
        | combine(FaidxIndex.out, by: 0)
        | map { it -> it.tail() }
        | VariantCall

    metadata_ch
        | map { meta -> [meta.id, meta.minimum_coverage] }
        | combine(VariantCall.out, by: 0)
        | VariantsStats

    VariantsStats.out
        | filter { it -> !(it[1].text =~ '\tTotal\tFrameshifts\t0') }
        | view { it -> "<${it[0]}> has frameshifts" }
        | collectFile(storeDir: "${workflow.outputDir}/warnings") { it -> ['frameshifts.txt', "${it[0]}\n"] }

    // Create consensuses
    metadata_ch
        | map { meta -> [meta.id, meta.reference, meta.id, [name: meta.name, minimum_coverage: meta.minimum_coverage]] }
        | combine(VariantCall.out, by: 0)
        | combine(GenomeCov.out, by: 0)
        | map { it -> it.tail() }
        | combine(References, by: 0)
        | map { it -> it.tail() }
        | Consenser

    metadata_ch
        | map { meta -> [meta.id, meta.group] }
        | combine(Consenser.out.degenerated, by: 0)
        | map { it -> it.tail() }
        | groupTuple(by: 0)
        | ConcatenateConensus

    if (params.virus == 'AIV') {
        AIVGetSubtype(BlastN.out)
            | collectFile(storeDir: "${workflow.outputDir}/results", sort: true) { it -> ['subtypes.tsv', "${it[0].replaceFirst(/_ref/, '')}\t${it[1].text}"] }

        UpdateFluMut()
        FluMut(ConcatenateConensus.out, UpdateFluMut.out)
        UpdateGenin2()
        Genin2(ConcatenateConensus.out, UpdateGenin2.out)
    }

    if (params.virus == 'SIV') {
        AIVGetSubtype(BlastN.out)
            | collectFile(storeDir: "${workflow.outputDir}/results", sort: true) { it -> ['subtypes.tsv', "${it[0].replaceFirst(/_ref/, '')}\t${it[1].text}"] }
    }

    Channel.topic('statistics')
        | map { it -> [it[1]] }
        | flatten
        | collectFile(name: 'statistics.tsv', storeDir: "${workflow.outputDir}/results")

    publish:
    fastqc_raw = FastQCRaw.out
    fastqc_clean = FastQCClean.out
    cutadapt = Cutadapt.out
    reference = PrepareReference.out
    reference_composition = GetReferenceNames.out
    bwa = BWAMem.out
    bwa_index = BamIndex.out
    coverage = GenomeCov.out
    tacos = Tacos.out
    vcfs = VariantCall.out
    consensus_deg = Consenser.out.degenerated
    consensus_undeg = Consenser.out.non_degenerated
    results = Channel.topic("results")
}

output {
    fastqc_raw {
        path { sample ->
            sample[1][0] >> "reads_quality/raw/${sample[0]}_R1.html"
            sample[1][1] >> "reads_quality/raw/${sample[0]}_R2.html"
        }
    }
    fastqc_clean {
        path { sample ->
            sample[1][0] >> "reads_quality/clean/${sample[0]}_R1.html"
            sample[1][1] >> "reads_quality/clean/${sample[0]}_R2.html"
        }
    }
    cutadapt {
        path { sample ->
            sample[1][0] >> "cleaned_reads/${sample[0]}_R1_cleaned.fastq.gz"
            sample[1][1] >> "cleaned_reads/${sample[0]}_R2_cleaned.fastq.gz"
        }
        mode 'symlink'
    }
    reference {
        path { sample -> sample[1] >> "alignments/references/${sample[0]}.fa" }
    }
    reference_composition {
        path { sample -> sample[1] >> "alignments/references/headers/${sample[0]}_reference_headers.txt" }
    }
    bwa {
        path { sample -> sample[1] >> "alignments/${sample[0]}.bam" }
    }
    bwa_index {
        path { sample -> sample[1] >> "alignments/${sample[0]}.bam.bai" }
    }
    coverage {
        path { sample -> sample[1] >> "coverage/raw/${sample[0]}_coverage.tsv" }
    }
    tacos {
        path { sample -> sample[1] >> "coverage/${sample[0]}_coverage.pdf" }
    }
    vcfs {
        path { sample -> sample[1] >> "vcfs/${sample[0]}.vcf" }
    }
    consensus_deg {
        path { sample -> sample[1] >> "consensus/${sample[0]}_consensus.fa" }
    }
    consensus_undeg {
        path { sample -> sample[1] >> "consensus/unambiguous/${sample[0]}_consensus.fa" }
    }
    results {
        path "results/"
    }
}
