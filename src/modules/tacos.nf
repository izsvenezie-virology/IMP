process Tacos {
    tag "${id}"
    errorStrategy 'ignore'

    memory '1 GB'
    time '1h'

    input:
    tuple val(id), path(coverage), val(min_coverage)

    output:
    tuple val(id), path("coverage.pdf")

    script:
    sample_id = id.split("__")[0]
    """
    tacos -m ${min_coverage} -s ${sample_id} ${coverage} coverage.pdf
    """
}
