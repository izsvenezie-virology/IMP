process GenomeCov {
    tag "$meta.sample"
    publishDir "coverage/raw", saveAs: { "${meta.sample}__${meta.reference}.tsv" }, mode: 'copy'

    memory '500 MB'
    time '30s'

    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path('*')
    
    """
    bedtools genomecov -d -ibam $bam >coverage.tsv
    """
}
