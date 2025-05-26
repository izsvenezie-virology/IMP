#! /usr/bin/env nextflow

nextflow.preview.output = true

include {
    GetReference ;
    ConcatenateConensus as CC_group ;
    ConcatenateConensus as CC_run
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
    Consenser as DegeneratedConsensus ;
    Consenser as NonDegeneratedConsensus
} from './modules/consenser.nf'
include {
    Cutadapt
} from './modules/cutadapt.nf'
include {
    FastQC as FastQCRaw ;
    FastQC as FastQCClean
} from './modules/fastqc.nf'
include {
    UpdateFluMut ;
    FluMut
} from './modules/flumut.nf'
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
    UpdateGenin2 ;
    Genin2
} from './modules/genin2.nf'
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
    AIVSubtype
} from './workflows/aiv_subtype.nf'

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

    // BLAST DB channel
    references_db = file(params.references_database, checkIfExists: true)

    // Adapters channel
    adapters = file(params.adapters, checkIfExists: true)

    // Samples channels creation
    Channel.fromFilePairs("${params.raw_reads_folder}/*_R{1,2}*fastq.gz")
        | map { it -> [it[0].split("_S")[0], it.tail().flatten()] }
        | set { raw_reads }

    Channel.fromPath(params.samples_metadata)
        | splitCsv(header: true, sep: '\t')
        | multiMap { it ->
            // If a reference is not provided set the file to null, create a custom id and set the subset > 0.
            // If a reference is provided and the file exists use it and set the simpleName as ID. The subset must be null.
            // Else if the file does not exists set the file to null and use the simpleName as ID. The subset must be null.
            def reference_file = file(it.Reference ?: 'non_existing_file').exists() ? file(it.Reference) : null
            def reference_id = it.Reference ? file(it.Reference).simpleName : "${it.Sample}_ref"

            def subset = it.Reference ? null : it.Subset ?: 0.2

            def primers_file = file(it.Primers ?: params.null_file, checkIfExists: true)
            def primers_id = primers_file.simpleName
            metadata: [
                id: "${it.Sample}__${reference_id}",
                sample: it.Sample,
                name: it.Name,
                primers: primers_id,
                reference: reference_id,
                subset: subset,
                group: it.Group,
                minimum_coverage: it.MinimumCoverage ?: 20,
            ]
            primers: [primers_id, primers_file]
            references: [reference_id, reference_file]
        }
        | set { samples_config_ch }

    metadata_ch = samples_config_ch.metadata
    samples_config_ch.primers
        | unique
        | CreateCutadaptPrimers

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
    FastQCRaw(to_fastqc_raw_ch, 'raw')
    FastQCClean(Cutadapt.out, 'clean')

    // References collection channel creation
    samples_config_ch.references
        | filter { it -> it[1] }
        | set { provided_references }

    // Find missing references
    metadata_ch
        | filter { it -> it.subset }
        | first
        | map { references_db }
        | MakeBlastDb
    // create blastDB only if at least one reference must be found
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
        | collectFile(storeDir: 'warnings') { it -> ['empty_references.txt', "${it[0]}\t${it[1]}\n"] }

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
        | Tacos

    // GATK best practices
    FixBam(BWAMem.out)

    CleanBam(FixBam.out)

    metadata_ch
        | map { meta -> [meta.id, meta.reference, meta.id] }
        | combine(CleanBam.out, by: 0)
        | map { it -> it.tail() }
        | combine(References, by: 0)
        | combine(BWAIndex.out, by: 0)
        | map { it -> it.tail() }
        | Viterbi

    MDSort(Viterbi.out)

    MarkDuplicates(MDSort.out)
        | MDBamIndex

    metadata_ch
        | map { meta -> [meta.id, meta.reference, meta.id] }
        | combine(MarkDuplicates.out, by: 0)
        | combine(MDBamIndex.out, by: 0)
        | map { it -> it.tail() }
        | combine(References, by: 0)
        | combine(FaidxIndex.out, by: 0)
        | map { it -> it.tail() }
        | FakeVariantCall

    IndexFeatureFile(FakeVariantCall.out)

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
        | map { meta -> [meta.id, meta.reference, meta.id] }
        | combine(IndelQual.out, by: 0)
        | combine(IQBamIndex.out, by: 0)
        | map { it -> it.tail() }
        | combine(References, by: 0)
        | combine(FaidxIndex.out, by: 0)
        | map { it -> it.tail() }
        | VariantCall

    // Create consensuses
    metadata_ch
        | map { meta -> [meta.id, meta.reference, meta.id, [name: meta.name, minimum_coverage: meta.minimum_coverage]] }
        | combine(VariantCall.out, by: 0)
        | combine(GenomeCov.out, by: 0)
        | map { it -> it.tail() }
        | combine(References, by: 0)
        | map { it -> it.tail() }
        | set { vcf_coverage_reference }
    DegeneratedConsensus(vcf_coverage_reference, true)
    NonDegeneratedConsensus(vcf_coverage_reference, false)

    metadata_ch
        | map { meta -> [meta.id, meta.group] }
        | combine(DegeneratedConsensus.out.consensus, by: 0)
        | map { it -> it.tail() }
        | groupTuple(by: 0)
        | CC_group

    if (params.virus == 'AIV') {
        AIVSubtype(Cutadapt.out)
            | map { it -> [it[1]] }
            | flatten
            | collectFile(name: 'subtypes.tsv', sort: { file -> file.text }, storeDir: 'results')

        CC_group.out
            | map { it -> it[1] }
            | toList
            | map { it -> [params.run, it] }
            | CC_run

        UpdateFluMut()
        FluMut(CC_group.out, UpdateFluMut.out)
        UpdateGenin2()
        Genin2(CC_run.out, UpdateGenin2.out)
    }

    publish:
    fastqc = Channel.topic('reads_quality')
    cutadapt = Cutadapt.out
    reference = GetReference.out
    bwa = BWAMem.out
    bwa_index = BamIndex.out
    coverage = GenomeCov.out
    tacos = Tacos.out
    vcfs = VariantCall.out
    consensus_deg = DegeneratedConsensus.out.consensus
    consensus_undeg = NonDegeneratedConsensus.out.consensus
    conensus_segments_deg = DegeneratedConsensus.out.segments
    consensus_segments_undeg = NonDegeneratedConsensus.out.segments
    results = Channel.topic("results")
}

output {
    fastqc {
        path { sample ->
            sample[2][0] >> "reads_quality/${sample[1]}/${sample[0]}_R1.html"
            sample[2][1] >> "reads_quality/${sample[1]}/${sample[0]}_R2.html"
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
        path { sample ->
            sample[1] >> "refs/${sample[0]}.fa"
        }
    }
    bwa {
        path { sample ->
            sample[1] >> "alignments/${sample[0]}.bam"
        }
    }
    bwa_index {
        path { sample ->
            sample[1] >> "alignments/${sample[0]}.bam.bai"
        }
    }
    coverage {
        path { sample ->
            sample[1] >> "coverage/raw/${sample[0]}_coverage.tsv"
        }
    }
    tacos {
        path { sample ->
            sample[1] >> "coverage/${sample[0]}_coverage.pdf"
        }
    }
    vcfs {
        path { sample ->
            sample[1] >> "vcfs/${sample[0]}.vcf"
        }
    }
    consensus_deg {
        path "consensus/"
    }
    consensus_undeg {
        path "consensus/unambiguous/"
    }
    conensus_segments_deg {
        path "consensus/segments/"
    }
    consensus_segments_undeg {
        path "consensus/unambiguous/segments/"
    }
    results {
        path "results/"
    }
}
