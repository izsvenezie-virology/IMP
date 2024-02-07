process GetReferenceNames{
    tag "$meta"

    memory '500 MB'
    time '5m'

    input:
        tuple val(meta), path(best_hits)
    output:
        tuple val(meta), path('*')

    """
    extract_reference.py $best_hits >ref_names.txt
    """
}

process CreateCutadaptPrimers{
    tag "$id"

    memory '500 MB'
    time '30s'

    input:
        tuple val(id), path(primers_tsv)
    output:
        tuple val(id), path('*')

    """
    create_cutadapt_primers.py $primers_tsv
    """
}
