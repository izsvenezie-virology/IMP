process DictIndex {
    tag "${id}"

    cpus 4
    memory '1 GB'
    time '1h'

    input:
    tuple val(id), path(reference)

    output:
    tuple val(id), path('*')

    script:
    """
    gatk --java-options "-XX:ConcGCThreads=1" CreateSequenceDictionary -R ${reference}
    """
}

process IndexFeatureFile {
    tag "${id}"

    cpus 4
    memory '1 GB'
    time '1h'

    input:
    tuple val(id), path(feature_file)

    output:
    tuple val(id), path('*')

    script:
    """
    gatk --java-options "-XX:ConcGCThreads=1" IndexFeatureFile -I ${feature_file}
    """
}

process FixBam {
    tag "${id}"

    cpus 4
    memory '20 GB'
    time '1h'

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("fixed.bam")

    script:
    """
    gatk --java-options "-XX:ConcGCThreads=1" FixMateInformation -I ${bam} -O fixed.bam
    """
}

process CleanBam {
    tag "${id}"

    cpus 4
    memory '5 GB'
    time '1h'

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("cleaned.bam")

    script:
    """
    gatk --java-options "-XX:ConcGCThreads=1" CleanSam -I ${bam} -O cleaned.bam
    """
}

process MarkDuplicates {
    tag "${id}"

    cpus 4
    memory '25 GB'
    time '1h'

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path("marked_duplicates.bam")

    script:
    """
    gatk --java-options "-XX:ConcGCThreads=1" MarkDuplicates -I ${bam} -O marked_duplicates.bam -M marked_duplicates.metrics
    """
}

process BaseRecalibrator {
    tag "${id}"

    cpus 4
    memory '5 GB'
    time '1h'

    input:
    tuple val(id), path(bam), path(reference), path(ref_idx), path(dict_idx), path(known_sites), path(feature_idx)

    output:
    tuple val(id), path("recalibration_report.txt")

    script:
    """
    gatk --java-options "-XX:ConcGCThreads=1" BaseRecalibrator -R ${reference} -I ${bam} --known-sites ${known_sites} -O recalibration_report.txt
    """
}

process ApplyBQSR {
    tag "${id}"

    cpus 4
    memory '5 GB'
    time '1h'

    input:
    tuple val(id), path(bam), path(recalibration)

    output:
    tuple val(id), path("bqsr.bam")

    script:
    """
    gatk --java-options "-XX:ConcGCThreads=1" ApplyBQSR -I ${bam} -O bqsr.bam -bqsr ${recalibration}
    """
}
