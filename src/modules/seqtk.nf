process FastqToFasta {
    tag "$meta"
    
    input:
        tuple val(meta), val (subset), path(reads)
    output:
        tuple val(meta), path('*')
    
    """
    cat $reads | seqtk seq -A -f $subset - >reads.fa
    """
}
