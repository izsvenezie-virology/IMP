process FastQC {
    cpus 2

    publishDir "reads_quality/$type", saveAs: { fileName -> fileName.contains("_R1_") ? "${meta.sample}_R1_${type}.html" : "${meta.sample}_R2_${type}.html" } 

    input:
        tuple val(meta), path(reads)
        val(type) // "raw"/"clean" before/after reads clean up
        
    output:
        tuple val(meta), path('*html')

    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}
