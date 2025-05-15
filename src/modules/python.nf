process GetReferenceNames {
    tag "${id}"

    memory '500 MB'
    time '5m'

    input:
    tuple val(id), path(best_hits)

    output:
    tuple val(id), path('*')

    script:
    """
    extract_reference.py ${best_hits} >ref_names.txt
    """
}

process CreateCutadaptPrimers {
    tag "${id}"

    memory '500 MB'
    time '30s'

    input:
    tuple val(id), path(primers_tsv)

    output:
    tuple val(id), path('*')

    script:
    """
    create_cutadapt_primers.py ${primers_tsv}
    """
}

process PrepareReference {
    tag "${id}"

    memory '500 MB'
    time '30s'

    input:
    tuple val(id), path('raw_reference.fa')

    output:
    tuple val(id), path('*')

    script:
    """
    prepare_reference.py raw_reference.fa >reference.fa
    """
}
