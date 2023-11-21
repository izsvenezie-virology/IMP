process DictIndex{
    input:
        tuple val(meta), path(reference)
    output:
        tuple val(meta), path(reference, includeInputs: true), path('*')

    """
    gatk CreateSequenceDictionary -R $reference
    """
}
