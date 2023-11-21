process GetReference{
    input:
        tuple val(meta), path(ref_names)
        path(db_fasta)
    output:
        tuple val(meta), path('*')

    """
    grep --no-group-separator -A 1 -f $ref_names $db_fasta | sed -E 's/.+\|/>/g' >reference.fa
    """
    
}
