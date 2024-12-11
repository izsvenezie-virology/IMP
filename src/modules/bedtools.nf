process GenomeCov {
    tag "$id.sample"
    publishDir "coverage/raw", saveAs: { "${id.sample}__${id.reference}.tsv" }, mode: 'copy'

    memory '500 MB'
    time '30s'

    input:
        tuple val(id), path(bam)
    output:
        tuple val(id), path('*')
    
    """
    bedtools genomecov -d -ibam $bam >coverage.tsv
    """
}
