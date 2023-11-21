process DictIndex{
    input:
        tuple val(meta), path('reference.fa')
    output:
        tuple val(meta), path('reference.fa', includeInputs: true), path('*')

    """
    gatk CreateSequenceDictionary -R reference.fa
    """
}
