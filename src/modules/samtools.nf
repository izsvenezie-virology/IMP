process FaidxIndex{
    input:
        tuple val(meta), path('reference.fa')
    output:
        tuple val(meta), path('reference.fa', includeInputs: true), path('*')

    """
    samtools faidx reference.fa
    """
}
