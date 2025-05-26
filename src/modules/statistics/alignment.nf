process MappedCount {
    tag "${id}"

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path('mapped_reads_count.txt'), topic: 'alignment_stats'

    script:
    """
    echo "${id}\t-\tMapped reads\t\$(samtools view -c -F 260 ${bam})" >mapped_reads_count.txt
    """
}

process UnmappedCount {
    tag "${id}"

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path('unmapped_reads_count.txt'), topic: 'alignment_stats'

    script:
    """
    echo "${id}\t-\tUnmapped reads\t\$(samtools view -c -f 4 ${bam})" >unmapped_reads_count.txt
    """
}
