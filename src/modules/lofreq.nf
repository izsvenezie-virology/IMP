process Viterbi {
    tag "${id}"

    memory '500 MB'
    time '1h'

    input:
    tuple val(id), path(bam), path(reference), path(reference_index)

    output:
    tuple val(id), path("viterbi.bam")

    script:
    """
    lofreq viterbi -f ${reference} -o viterbi.bam ${bam} 
    """
}

process Call {
    tag "${id}"
    label 'multiThread'

    memory '5 GB'
    time '1h'

    input:
    tuple val(id), path(bam), path(bam_index), path(reference), path(reference_index)

    output:
    tuple val(id), path("variants.vcf")

    script:
    """
    lofreq call-parallel --call-indels --pp-threads ${task.cpus} -f ${reference} -o variants.vcf ${bam}
    """
}

process IndelQual {
    tag "${id}"

    memory '500 MB'
    time '1h'

    input:
    tuple val(id), path(bam), path(reference)

    output:
    tuple val(id), path("indelqual.bam")

    script:
    """
    lofreq indelqual --dindel -f ${reference} -o indelqual.bam ${bam}
    """
}
