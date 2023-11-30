process FastQC {
    tag "$meta.sample"
    cpus 2

    publishDir "reads_quality/$type", saveAs: { "${meta.sample}_R1_${type}.html" }, mode: 'copy', pattern: '*_R1_*'
    publishDir "reads_quality/$type", saveAs: { "${meta.sample}_R2_${type}.html" }, mode: 'copy', pattern: '*_R2_*'

    input:
        tuple val(meta), path(reads)
        val(type) // "raw"/"clean" before/after reads clean up
        
    output:
        tuple val(meta), path('*html')

    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}
