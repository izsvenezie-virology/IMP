process CoverageStats {
    tag "${id}"

    input:
    tuple val(id), val(minimum_coverage), path(coverage)

    output:
    tuple val(id), path('coverage_stats.txt'), topic: 'statistics'

    script:
    """
    statistics_coverage.py ${id} ${minimum_coverage} ${coverage} >coverage_stats.txt
    """
}

process AlignmentStats {
    tag "${id}"

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path('alignment_stats.txt'), topic: 'statistics'

    script:
    """
    statistics_alignment.py ${id} ${bam} >alignment_stats.txt
    """
}

process ReadsStats {
    tag "${id}"

    input:
    tuple val(id), path(reads)
    val type

    output:
    tuple val(id), path('alignment_stats.txt'), topic: 'statistics'

    script:
    """
    statistics_reads.py ${id} ${type} ${reads} >alignment_stats.txt
    """
}
