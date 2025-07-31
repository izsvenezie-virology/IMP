process GetReferenceNames {
    tag "${id}"

    memory '500 MB'
    time '1h'

    input:
    tuple val(id), path(best_hits)

    output:
    tuple val(id), path("ref_names.txt")

    script:
    """
    extract_reference.py ${best_hits} >ref_names.txt
    """
}

process CreateCutadaptPrimers {
    tag "${id}"

    memory '500 MB'
    time '1h'

    input:
    tuple val(id), path(primers_tsv)

    output:
    tuple val(id), path("primers_??.fa")

    script:
    """
    create_cutadapt_primers.py ${primers_tsv}
    """
}

process PrepareReference {
    tag "${id}"

    memory '500 MB'
    time '1h'

    input:
    tuple val(id), path(reference)

    output:
    tuple val(id), path("cleaned_reference.fa")

    script:
    """
    prepare_reference.py ${reference} >cleaned_reference.fa
    """
}

process FindSubtypes {
    tag "${id}"

    memory '500 MB'
    time '1h'

    input:
    tuple val(id), path(coverage)

    output:
    tuple val(id), path("alignment_stats.tsv")

    script:
    """
    find_subtypes.py ${id} ${coverage} >alignment_stats.tsv
    """
}
