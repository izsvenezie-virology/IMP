process FastQC {
    tag "${id}"

    publishDir "reads_quality/${type}", saveAs: { "${id}_R1_${type}.html" }, mode: 'copy', pattern: '*_R1_*'
    publishDir "reads_quality/${type}", saveAs: { "${id}_R2_${type}.html" }, mode: 'copy', pattern: '*_R2_*'

    cpus 2
    memory '2 GB'
    time '5m'

    input:
    tuple val(id), path(reads)
    val type

    output:
    tuple val(id), path('*html')

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}
