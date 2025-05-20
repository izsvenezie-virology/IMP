process GetReference {
    tag "${id}"

    memory '500 MB'
    time '30s'

    input:
    tuple val(id), path(ref_names)
    path db_fasta

    output:
    tuple val(id), path('reference.fa')

    script:
    """
    grep --no-group-separator -A 1 -f ${ref_names} ${db_fasta} | sed -E 's/^>.+\\|/>/g' >reference.fa
    """
}

process ConcatenateConensus {
    tag "${id}"

    memory '500MB'
    time '30s'

    input:
    tuple val(id), path(consensuses)

    output:
    tuple val(id), path("${id}_consensus.fa"), topic: 'results'

    script:
    """
    cat * > ${id}_consensus.fa
    """
}
