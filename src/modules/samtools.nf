process FaidxIndex{
    input:
        tuple val(meta), path(reference)
    output:
        tuple val(meta), path(reference, includeInputs: true), path('*')

    """
    samtools faidx $reference
    """
}

process Sort{
    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path('*')
    """
    samtools sort -O bam -o sorted.bam $bam
    """
}
