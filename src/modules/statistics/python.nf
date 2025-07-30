process CoverageStats {
    tag "${id}"

    memory '2 GB'
    time '1h'

    input:
    tuple val(id), val(minimum_coverage), path(coverage)

    output:
    tuple val(id), path('statistics_coverage.txt'), topic: 'statistics'

    script:
    """
    statistics_coverage.py ${id} ${minimum_coverage} ${coverage} >statistics_coverage.txt
    """
}

process AlignmentStats {
    tag "${id}"

    memory '2 GB'
    time '1h'

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path('statistics_alignment.txt'), topic: 'statistics'

    script:
    """
    statistics_alignment.py ${id} ${bam} >statistics_alignment.txt
    """
}

process ReadsStats {
    tag "${id}"

    memory '2 GB'
    time '1h'

    input:
    tuple val(id), path(reads)
    val type

    output:
    tuple val(id), path('statistics_reads.txt'), topic: 'statistics'

    script:
    """
    statistics_reads.py ${id} ${type} ${reads} >statistics_reads.txt
    """
}

process VariantsStats {
    tag "${id}"

    memory '2 GB'
    time '1h'

    input:
    tuple val(id), val(minimum_coverage), path(vcf)

    output:
    tuple val(id), path('statistics_variants.txt'), topic: 'statistics'

    script:
    """
    statistics_variants.py ${id} ${minimum_coverage} ${vcf} >statistics_variants.txt
    """
}
