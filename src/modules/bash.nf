process GetReference{
    publishDir "refs", saveAs: { "${meta}.fa" }, mode: 'copy'

    input:
        tuple val(meta), path(ref_names)
        path(db_fasta)
    output:
        tuple val(meta), path('*')

    """
    grep --no-group-separator -A 1 -f ${ref_names} ${db_fasta} | sed -E 's/^>.+\\|/>/g' >reference.fa
    """    
}

process ConcatFiles{
    publishDir "results", saveAs: { "${name}.fa" }, mode: 'copy'

    input:
        tuple val(name), path(files)
    output:
        tuple val(name), path("*")
    """
    cat $files >all.fa
    """
}
