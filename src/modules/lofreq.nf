process Viterbi{
    tag "$id.sample"
    
    memory '500 MB'
    time '5m'

    input:
        tuple val(id), path(bam), path(reference), path(reference_index)
    output:
        tuple val(id), path("viterbi.bam")
    
    """
    lofreq viterbi -f $reference -o viterbi.bam $bam 
    """
}

process Call{
    tag "$id.sample"
    label 'multiThread'
    publishDir 'vcfs', saveAs: { "${id.sample}__${id.reference}.vcf" }, mode: 'copy', enabled: "$call_indels"

    memory '5 GB'
    time '15m'

    input:
        tuple val(id), path(bam), path(bam_index), path(reference), path(reference_index)
        val(call_indels)
    output:
        tuple val(id), path("variants.vcf")

    script:
    def call_indels_opt = call_indels ? '--call-indels' : '' // If call indels is true the call indels option is set
    """
    lofreq call-parallel --pp-threads $task.cpus -f $reference -o variants.vcf $bam $call_indels_opt
    """
}

process IndelQual{
    tag "$id.sample"
    
    memory '500 MB'
    time '5m'

    input:
        tuple val(id), path(bam), path(reference)
    output:
        tuple val(id), path("*.bam")
    
    """
    lofreq indelqual --dindel -f $reference -o indelqual.bam $bam
    """
}
