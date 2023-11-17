process CoveragePlotter {
    errorStrategy 'ignore'
    publishDir "coverage", saveAs: { fileName -> "${meta.sample}__${meta.reference}.pdf" }, mode: 'copy'

    input:
        tuple val(meta), path(coverage)
    output:
        tuple val(meta), path('coverage.pdf')
    
    """
    coverplotter -o coverage.pdf $coverage
    """
}
