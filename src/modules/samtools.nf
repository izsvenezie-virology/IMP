process FaidxIndex{
    tag "$id"

    memory '500 MB'
    time '30s'

    input:
        tuple val(id), path(reference)
    output:
        tuple val(id), path('*')

    """
    samtools faidx $reference
    """
}

process Sort{
    tag "$id.sample"

    memory '5 GB'
    time '5m'

    input:
        tuple val(id), path(bam)
    output:
        tuple val(id), path('*')
    """
    samtools sort -O bam -o sorted.bam $bam
    """
}

process BamIndex{
    tag "$id.sample"
    publishDir "alignments", saveAs: { "${id.sample}__${id.reference}.bai" }, mode: 'copy'

    memory '500 MB'
    time '30s'

    input:
        tuple val(id), path(bam)
    output:
        tuple val(id), path('*')
    """
    samtools index $bam
    """
}
