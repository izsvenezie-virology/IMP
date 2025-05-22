process GenomeCov {
    tag "${id}"

    memory '500 MB'
    time '30s'

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path('coverage.tsv')

    script:
    """
    bedtools genomecov -d -ibam ${bam} >coverage.tsv
    """
}
