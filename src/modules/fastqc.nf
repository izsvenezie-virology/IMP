process FastQC {
    tag "${id}"

    cpus 2
    memory '2 GB'
    time '5m'

    input:
    tuple val(id), path(reads)
    val type

    output:
    tuple val(id), val(type), path('*.html'), topic: reads_quality

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}
