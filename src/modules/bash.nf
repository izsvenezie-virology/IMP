process GetReference{
    tag "$id"
    publishDir "refs", saveAs: { "${id}.fa" }, mode: 'copy'

    memory '500 MB'
    time '30s'

    input:
        tuple val(id), path(ref_names)
        path(db_fasta)
    output:
        tuple val(id), path('*')

    """
    grep --no-group-separator -A 1 -f ${ref_names} ${db_fasta} | sed -E 's/^>.+\\|/>/g' >reference.fa
    """    
}

process ConcatenateConensus{
    tag "$id"
    publishDir "results", mode: 'copy'

    memory '500MB'
    time '30s'

    input:
        tuple val(id), path(consensuses)
    output:
        tuple val(id), path("${id}_consensus.fa")

    """
    cat * > ${id}_consensus.fa
    """
}
