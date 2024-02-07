process FastqToFasta {
    tag "$id"
    
    memory '500 MB'
    time '1m'

    input:
        tuple val(id), val (subset), path(reads)
    output:
        tuple val(id), path('*')
    
    """
    cat $reads | seqtk seq -A -f $subset - >reads.fa
    """
}
