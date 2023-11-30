process Viterbi{
    tag "$meta.sample"
    
    input:
        tuple val(meta), path(bam), path(reference), path(reference_index)
    output:
        tuple val(meta), path("viterbi.bam")
    
    """
    lofreq viterbi -f $reference -o viterbi.bam $bam 
    """
}

process Call{
    tag "$meta.sample"
    label 'multiThread'
    publishDir 'vcfs', saveAs: { "${meta.sample}__${meta.reference}.vcf" }, mode: 'copy', enabled: "$call_indels"

    input:
        tuple val(meta), path(bam), path(reference)
        val(call_indels)
    output:
        tuple val(meta), path("variants.vcf")

    script:
    def call_indels_opt = call_indels ? '--call-indels' : '' // If call indels is true the call indels option is set
    """
    lofreq call-parallel --pp-threads $task.cpus -f $reference -o variants.vcf $bam $call_indels_opt
    """
}
