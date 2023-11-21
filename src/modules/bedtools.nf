process GenomeCov {
    publishDir "coverage/raw", saveAs: { "${meta.sample}__${meta.reference}.tsv" }, mode: 'symlink'

    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path('*')
    
    """
    bedtools genomecov -d -ibam $bam >coverage.tsv
    """
}
