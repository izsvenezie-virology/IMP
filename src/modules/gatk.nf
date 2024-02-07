process DictIndex{
    tag "$id"

    memory '1 GB'
    time '30s'

    input:
        tuple val(id), path(reference)
    output:
        tuple val(id), path('*')

    """
    gatk --java-options "-XX:ConcGCThreads=1" CreateSequenceDictionary -R $reference
    """
}

process IndexFeatureFile{
    tag "$id.sample"

    memory '1 GB'
    time '30s'

    input:
        tuple val(id), path(feature_file)
    output:
        tuple val(id), path('*')

    """
    gatk --java-options "-XX:ConcGCThreads=1" IndexFeatureFile -I $feature_file
    """
}

process FixBam{
    tag "$id.sample"

    memory '20 GB'
    time '5m'

    input:
        tuple val(id), path(bam)
    output:
        tuple val(id), path("*")

    """
    gatk --java-options "-XX:ConcGCThreads=1" FixMateInformation -I ${bam} -O fixed.bam
    """
}

process CleanBam{
    tag "$id.sample"

    memory '5 GB'
    time '5m'

    input:
        tuple val(id), path(bam)
    output:
        tuple val(id), path("*")

    """
    gatk --java-options "-XX:ConcGCThreads=1" CleanSam -I ${bam} -O cleaned.bam
    """
}

process MarkDuplicates{
    tag "$id.sample"

    memory '25 GB'
    time '5m'

    input:
        tuple val(id), path(bam)
    output:
        tuple val(id), path("*.bam")

    """
    gatk --java-options "-XX:ConcGCThreads=1" MarkDuplicates -I ${bam} -O marked_duplicates.bam -M marked_duplicates.metrics
    """
}

process BaseRecalibrator{
    tag "$id.sample"

    memory '5 GB'
    time '5m'

    input:
        tuple val(id), path(bam), path(reference), path(ref_idx), path(dict_idx), path(known_sites), path(feature_idx)
    output:
        tuple val(id), path("*")

    """
    gatk --java-options "-XX:ConcGCThreads=1" BaseRecalibrator -R $reference -I $bam --known-sites $known_sites -O recalibration_report.txt
    """
}

process ApplyBQSR{
    tag "$id.sample"
    
    memory '5 GB'
    time '5m'

    input:
        tuple val(id), path(bam), path(recalibration)
    output:
        tuple val(id), path("*.bam"), path("*.bai")

    """
    gatk --java-options "-XX:ConcGCThreads=1" ApplyBQSR -I $bam -O bqsr.bam -bqsr $recalibration
    """
}
