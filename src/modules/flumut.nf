process UpdateFluMut {
    tag "FluMutDB"

    memory '500 MB'
    time '5m'

    output:
    val 'true'

    script:
    """
    flumut --update
    """
}


process FluMut {
    tag "${id}"

    memory '2 GB'
    time '30m'

    input:
    tuple val(id), path(fasta)
    val updated

    output:
    tuple val(id), path('*'), topic: 'results'

    script:
    """
    flumut -x ${id}_flumut.xlsm ${fasta}
    """
}
