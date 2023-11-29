process FaidxIndex{
    input:
        tuple val(meta), path(reference)
    output:
        tuple val(meta), path('*')

    """
    samtools faidx $reference
    """
}

process Sort{
    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path('*')
    """
    samtools sort -O bam -o sorted.bam $bam
    """
}

process BamIndex{
    publishDir "alignments", saveAs: { "${meta.sample}__${meta.reference}.bai" }

    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path('*')
    """
    samtools index $bam
    """
}
