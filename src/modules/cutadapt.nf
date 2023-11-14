process Cutadapt {
    publishDir "clean_reads", saveAs:  { fileName -> fileName.contains("_R1_") ? "${meta.sample}_R1_clean.fastq.gz" : "${meta.sample}_R2_clean.fastq.gz" }

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path('*_clean.fastq.gz')

    """
    cutadapt --interleaved --action=trim --pair-filter=any -q 20 -m 80 -o reads_R1_clean.fastq.gz -p reads_R2_clean.fastq.gz $reads
    """
}
