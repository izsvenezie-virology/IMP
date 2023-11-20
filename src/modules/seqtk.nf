process FastqToFasta {
    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path('reads.fa')
    
    """
    cat $reads | seqtk seq -A -f 0.3 - >reads.fa
    """
}
