process FastqToFasta {
    tag "$meta"
    
    memory '500 MB'
    time '1m'

    input:
        tuple val(meta), val (subset), path(reads)
    output:
        tuple val(meta), path('*')
    
    """
    cat $reads | seqtk seq -A -f $subset - >reads.fa
    """
}
