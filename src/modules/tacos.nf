process Tacos {
    tag "$id.sample"
    errorStrategy 'ignore'
    
    publishDir "coverage", saveAs: { "${id.sample}__${id.reference}.pdf" }, mode: 'copy'

    memory '1 GB'
    time '5m'

    input:
        tuple val(id), path(coverage)
    output:
        tuple val(id), path('*')
    
    """
    tacos -o coverage.pdf $coverage
    """
}
