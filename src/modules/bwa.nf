process BWAmem {
    publishDir "alignments", saveAs: { "${meta.sample}__${meta.reference}.bam" }, mode: 'symlink'

    input:
        tuple val(meta), path(reads), path(reference), path(index)
    output:
        tuple val(meta), path('*')
    
    """
    bwa mem -t $task.cpus -R '@RG\\tID:${meta.sample}\\tSM:${meta.sample}\\tPL:ILLUMINA' -M $reference $reads |
        samtools sort -@ $task.cpus -O bam -o sorted.bam
    """
}

process BWAIndex{
    input:
        tuple val(meta), path(reference)
    output:
        tuple val(meta), path(reference, includeInputs: true), path('*')
    
    """
    bwa index $reference
    """
}
