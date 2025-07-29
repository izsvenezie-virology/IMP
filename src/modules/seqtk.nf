process FastqToFasta {
    tag "${id}"

    memory '500 MB'
    time '1h'

    input:
    tuple val(id), val(subset), path(reads)

    output:
    tuple val(id), path("reads.fa")

    script:
    """
    cat ${reads} | seqtk seq -A -f ${subset} - >reads.fa
    """
}
