process GenomeCov {
    tag "${id.sample}"

    memory '500 MB'
    time '30s'

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path('*')

    script:
    """
    bedtools genomecov -d -ibam ${bam} >coverage.tsv
    """
}
