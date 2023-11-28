process Viterbi{
    input:
        tuple val(meta), path(bam), path(reference), path(reference_index)
    output:
        tuple val(meta), path("viterbi.bam")
    
    """
    lofreq viterbi -f $reference -o viterbi.bam $bam 
    """
}
