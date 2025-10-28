process FastQC {
    tag "${id}"

    cpus 2
    memory '2 GB'
    time '1h'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path('*.html')

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}
