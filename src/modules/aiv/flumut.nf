process UpdateFluMut {
    tag "FluMutDB"

    memory '500 MB'
    time '1h'

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
    time '1h'

    input:
    tuple val(id), path(fasta)
    val updated

    output:
    tuple val(id), path("${id}_flumut.xlsm"), topic: 'results'

    script:
    """
    flumut -x ${id}_flumut.xlsm ${fasta}
    """
}
