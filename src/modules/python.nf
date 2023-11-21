process GetReferenceNames{
    input:
        tuple val(meta), path(best_hits)
    output:
        tuple val(meta), path('*')

    """
    extract_reference.py $best_hits >ref_names.txt
    """
}
