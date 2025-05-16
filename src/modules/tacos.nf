process Tacos {
    tag "${id.sample}"
    errorStrategy 'ignore'

    memory '1 GB'
    time '5m'

    input:
    tuple val(id), path(coverage)

    output:
    tuple val(id), path('*')

    script:
    """
    tacos -o coverage.pdf ${coverage}
    """
}
