process CoveragePlotter {
    tag "$id.sample"
    errorStrategy 'ignore'
    
    publishDir "coverage", saveAs: { "${id.sample}__${id.reference}.pdf" }, mode: 'copy'

    memory '1 GB'
    time '30s'

    input:
        tuple val(id), path(coverage)
    output:
        tuple val(id), path('*')
    
    """
    coverplotter -o coverage.pdf $coverage
    """
}
