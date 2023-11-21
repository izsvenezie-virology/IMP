process FaidxIndex{
    input:
        tuple val(meta), path(reference)
    output:
        tuple val(meta), path(reference, includeInputs: true), path('*')

    """
    samtools faidx $reference
    """
}
