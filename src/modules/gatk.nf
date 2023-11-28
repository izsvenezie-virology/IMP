process DictIndex{
    input:
        tuple val(meta), path(reference)
    output:
        tuple val(meta), path('*')

    """
    gatk CreateSequenceDictionary -R $reference
    """
}

process FixBam{
    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path("*")

    """
    gatk FixMateInformation -I ${bam} -O fixed.bam
    """
}

process CleanBam{
    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path("*")

    """
    gatk CleanSam -I ${bam} -O cleaned.bam
    """
}

process MarkDuplicates{
    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path("*.bam")

    script:
    """
    gatk MarkDuplicates -I ${bam} -O marked_duplicates.bam -M marked_duplicates.metrics
    """
}
