process BWAmem {
    input:
        tuple val(meta), path(reads), path(reference), path(index)
    output:
        tuple val(meta), path('sorted.bam')
    
    """
    bwa mem -t $task.cpus -R '@RG\\tID:${meta.sample}\\tSM:${meta.sample}\\tPL:ILLUMINA' -M $reference $reads |
        samtools sort -@ $task.cpus -O bam -o sorted.bam
    """
}

process BWAIndex{
    debug true

    input:
        tuple val(meta), path('reference.fa')
    output:
        tuple val(meta), path('reference.fa', includeInputs: true), path('*')
    
    """
    bwa index reference.fa
    """
}
