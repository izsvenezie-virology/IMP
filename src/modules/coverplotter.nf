process CoveragePlotter {
    tag "$meta.sample"
    errorStrategy 'ignore'
    
    publishDir "coverage", saveAs: { "${meta.sample}__${meta.reference}.pdf" }, mode: 'copy'

    memory '1 GB'
    time '30s'

    input:
        tuple val(meta), path(coverage)
    output:
        tuple val(meta), path('*')
    
    """
    coverplotter -o coverage.pdf $coverage
    """
}
