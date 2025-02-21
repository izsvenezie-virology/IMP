process UpdateGenin2{
    tag "Genin2"

    memory '500 MB'
    time '5m'

    output:
        val 'true'

    script:
    """
    pip install --upgrade genin2
    """
}


process Genin2{
    tag "$id"
    publishDir "results", mode: 'copy'

    memory '2 GB'
    time '30m'

    input:
        tuple val(id), path(fasta)
        val(updated)
    output:
        tuple val(id), path('*')
    
    """
    genin2 -o ${id}_genin2.tsv $fasta
    """
}
