process UpdateGenin2 {
    tag "Genin2"

    memory '500 MB'
    time '1h'

    output:
    val 'true'

    script:
    """
    pip install --upgrade genin2
    """
}


process Genin2 {
    tag "${id}"

    memory '2 GB'
    time '1h'

    input:
    tuple val(id), path(fasta)
    val updated

    output:
    tuple val(id), path("${id}_genin2.tsv"), topic: 'results'

    script:
    """
    genin2 -o ${id}_genin2.tsv ${fasta}
    """
}
