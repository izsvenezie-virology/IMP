process Cutadapt {
    publishDir "clean_reads", saveAs:  { fileName -> fileName.startsWith("R1_") ? "${meta.sample}_R1_clean.fastq.gz" : "${meta.sample}_R2_clean.fastq.gz" }

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path('R?_clean.fastq.gz')

    """
    cutadapt --interleaved --action=trim --pair-filter=any -q 20 -m 80 -o R1_clean.fastq.gz -p R2_clean.fastq.gz $reads
    """
}
