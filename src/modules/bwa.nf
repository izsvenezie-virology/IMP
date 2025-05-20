process BWAMem {
    tag "${id}"
    label 'multiThread'

    memory '10 GB'
    time '5m'

    input:
    tuple val(id), path(reads), path(reference), path(index)

    output:
    tuple val(id), path('sorted.bam')

    script:
    """
    bwa mem -t ${task.cpus} -R '@RG\\tID:${id}\\tSM:${id}\\tPL:ILLUMINA' -M ${reference} ${reads} |
        samtools sort -@ ${task.cpus} -O bam -o sorted.bam
    """
}

process BWAIndex {
    tag "${id}"

    memory '50 MB'
    time '30s'

    input:
    tuple val(id), path(reference)

    output:
    tuple val(id), path("${reference}*")

    script:
    """
    bwa index ${reference}
    """
}
