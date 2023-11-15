process GenomeCov {
    publishDir "coverage", saveAs: { fileName -> "${meta.sample}__${meta.reference}.tsv" }, mode: 'symlink'

    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path('coverage.tsv')
    
    """
    bedtools genomecov -d -ibam $bam >coverage.tsv1
    """
}
