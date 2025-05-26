process CoverageStats {
    tag "${id}"

    input:
    tuple val(id), val(minimum_coverage), path(coverage)

    output:
    tuple val(id), path('coverage_stats.txt'), topic: 'alignment_stats'

    script:
    """
    coverage_stats.py ${id} ${minimum_coverage} ${coverage} >coverage_stats.txt
    """
}
