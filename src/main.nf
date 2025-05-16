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
    DegeneratedConsensus ;
    NonDegeneratedConsensus
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
        | map { [it[0].split("_S")[0], it[1]] }
        | set { raw_reads }

    Channel.fromPath(params.samples_metadata)
        | splitCsv(header: true, sep: '\t')
        | map {
            def reference = it.Reference ? file(it.Reference).simpleName : "${it.Sample}_ref"
            def primers = it.Primers ? file(it.Primers).simpleName : file(params.null_file).simpleName
            [
                it.Sample,
                [
                    id: [sample: it.Sample, reference: reference],
                    sample: it.Sample,
                    name: it.Name,
                    primers: primers,
                    primers_file: it.Primers,
                    reference: reference,
                    reference_file: it.Reference,
                    subset: it.Subset,
                    group: it.Group,
                ],
            ]
        }
        | set { metadata_ch }

    metadata_ch
        | map {
            if (it[1].primers_file) {
                [it[1].primers, file(it[1].primers_file, checkIfExists: true)]
            }
            else {
                [it[1].primers, file(params.null_file, checkIfExists: true)]
            }
        }
        | unique
        | CreateCutadaptPrimers

    // Clean reads
    raw_reads
        | join(metadata_ch)
        | map { [it[2].primers, it[2].sample, [it[2].minQual, it[2].minLen], it[1]] }
        | combine(CreateCutadaptPrimers.out, by: 0)
        | map { it.tail() }
        | set { to_cutadapt_ch }
    Cutadapt(to_cutadapt_ch, adapters)

    // Assess reads quality
    raw_reads
        | join(metadata_ch)
        | map { [it[0], it[1]] }
        | set { to_fastqcraw_ch }

    FastQCRaw(to_fastqcraw_ch, 'raw')
    FastQCClean(Cutadapt.out, 'clean')

    // References collection channel creation
    metadata_ch
        | combine(Cutadapt.out, by: 0)
        | branch {
            FindRef: it[1].reference_file == ''
            return [it[1].reference, it[1].subset ?: 1, it[2]]
            RefProvided: true
            return [it[1].reference, file(it[1].reference_file, checkIfExists: true)]
        }
        | set { ref_collect }

    // Find missing references
    ref_collect.FindRef
        | first
        | ifEmpty(false)
        | map { it ? true : false }
        | set { build_blast_db }
    MakeBlastDb(references_db, build_blast_db)
    // create blastDB only if at least one reference must be found
    FastqToFasta(ref_collect.FindRef)
    BlastN(FastqToFasta.out, MakeBlastDb.out)
    GetReferenceNames(BlastN.out)
    GetReference(GetReferenceNames.out, references_db)

    // References channel creation
    GetReference.out
        | filter { !it[1].isEmpty() }
        | mix(ref_collect.RefProvided)
        | unique
        | PrepareReference

    PrepareReference.out
        | set { References }
    // set channel

    // Reference index processes
    BWAIndex(References)
    FaidxIndex(References)
    DictIndex(References)

    // Reads alignment
    metadata_ch
        | map { [it[1].sample, it[1].reference, it[1].id] }
        | combine(Cutadapt.out, by: 0)
        | map { it.tail() }
        | combine(References, by: 0)
        | combine(BWAIndex.out, by: 0)
        | map { it.tail() }
        | BWAMem

    BamIndex(BWAMem.out, true)

    GenomeCov(BWAMem.out)
        | Tacos

    // GATK best practices
    FixBam(BWAMem.out)

    CleanBam(FixBam.out)

    CleanBam.out
        | map { [it[0].reference] + it }
        | combine(References, by: 0)
        | combine(BWAIndex.out, by: 0)
        | map { it.tail() }
        | Viterbi

    MDSort(Viterbi.out)

    MarkDuplicates(MDSort.out)

    MDBamIndex(MarkDuplicates.out, false)

    MarkDuplicates.out
        | combine(MDBamIndex.out, by: 0)
        | map { [it[0].reference] + it }
        | combine(References, by: 0)
        | combine(FaidxIndex.out, by: 0)
        | map { it.tail() }
        | set { md_reference }
    // set channel
    FakeVariantCall(md_reference, false)

    IndexFeatureFile(FakeVariantCall.out)

    MarkDuplicates.out
        | map { [it[0].reference] + it }
        | combine(References, by: 0)
        | combine(FaidxIndex.out, by: 0)
        | combine(DictIndex.out, by: 0)
        | map { it.tail() }
        | combine(FakeVariantCall.out, by: 0)
        | combine(IndexFeatureFile.out, by: 0)
        | BaseRecalibrator

    MarkDuplicates.out
        | combine(BaseRecalibrator.out, by: 0)
        | ApplyBQSR

    ApplyBQSR.out.bam
        | map { [it[0].reference] + it }
        | combine(References, by: 0)
        | map { it.tail() }
        | IndelQual

    IQBamIndex(IndelQual.out, false)

    IndelQual.out
        | combine(IQBamIndex.out, by: 0)
        | map { [it[0].reference] + it }
        | combine(References, by: 0)
        | combine(FaidxIndex.out, by: 0)
        | map { it.tail() }
        | set { to_variantcall_ch }
    // set channel
    VariantCall(to_variantcall_ch, true)

    // Create consensuses
    metadata_ch
        | map { [it[1].id, it[1].reference, it[1].id, [name: it[1].name]] }
        | combine(VariantCall.out, by: 0)
        | combine(GenomeCov.out, by: 0)
        | map { it.tail() }
        | combine(References, by: 0)
        | map { it.tail() }
        | set { vcf_coverage_reference }
    // set channel
    DegeneratedConsensus(vcf_coverage_reference)
    NonDegeneratedConsensus(vcf_coverage_reference)

    metadata_ch
        | map { [it[1].id, it[1].group] }
        | combine(DegeneratedConsensus.out.consensus, by: 0)
        | map { it.tail() }
        | groupTuple(by: 0)
        | CC_group

    if (params.virus == 'AIV') {
        AIVSubtype(Cutadapt.out)
            | map { [it[1]] }
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
            sample[1] >> "alignments/${sample[0].sample}__${sample[0].reference}.bam"
        }
    }
    bwa_index {
        path { sample ->
            sample[1] >> "alignments/${sample[0].sample}__${sample[0].reference}.bam.bai"
        }
    }
    coverage {
        path { sample ->
            sample[1] >> "coverage/raw/${sample[0].sample}__${sample[0].reference}_coverage.tsv"
        }
    }
    tacos {
        path { sample ->
            sample[1] >> "coverage/${sample[0].sample}__${sample[0].reference}_coverage.pdf"
        }
    }
    vcfs {
        path { sample ->
            sample[1] >> "vcfs/${sample[0].sample}__${sample[0].reference}.vcf"
        }
    }
    consensus_deg {
        path { sample ->
            sample[1] >> "consensus/${sample[0].sample}__${sample[0].reference}.fa"
        }
    }
    consensus_undeg {
        path { sample ->
            sample[1] >> "consensus/undegenerated/${sample[0].sample}__${sample[0].reference}.fa"
        }
    }
    results {
        path "results/"
    }
}
