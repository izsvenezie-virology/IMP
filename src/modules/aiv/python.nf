process AIVGetSubtype {
    tag "${id}"

    memory '500 MB'
    time '1h'

    input:
    tuple val(id), path(best_hits)

    output:
    tuple val(id), path("subtypes.txt")

    script:
    """
    aiv_get_subtypes.py ${best_hits} >subtypes.txt
    """
}
