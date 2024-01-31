process FaidxIndex{
    tag "$meta"

    memory '500 MB'
    time '30s'

    input:
        tuple val(meta), path(reference)
    output:
        tuple val(meta), path('*')

    """
    samtools faidx $reference
    """
}

process Sort{
    tag "$meta.sample"

    memory '5 GB'
    time '5m'

    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path('*')
    """
    samtools sort -O bam -o sorted.bam $bam
    """
}

process BamIndex{
    tag "$meta.sample"
    publishDir "alignments", saveAs: { "${meta.sample}__${meta.reference}.bai" }, mode: 'copy'

    memory '500 MB'
    time '30s'

    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path('*')
    """
    samtools index $bam
    """
}
