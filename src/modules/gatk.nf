process DictIndex{
    input:
        tuple val(meta), path(reference)
    output:
        tuple val(meta), path(reference, includeInputs: true), path('*')

    """
    gatk CreateSequenceDictionary -R $reference
    """
}

process FixBam{
    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*")

    script:
    """
    gatk FixMateInformation -I ${bam} -O fixed.bam
    """
}
