process Cutadapt {
    tag "$meta.sample"
    label 'multiThread'
    publishDir "clean_reads", saveAs: { "${meta.sample}_R1_clean.fastq.gz" }, mode: 'symlink', pattern: '*_R1_*'
    publishDir "clean_reads", saveAs: { "${meta.sample}_R2_clean.fastq.gz" }, mode: 'symlink', pattern: '*_R2_*'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path('*_clean.fastq.gz')

    """
    cutadapt --interleaved -j $task.cpus --action=trim --pair-filter=any -q 20 -m 80 -o reads_R1_clean.fastq.gz -p reads_R2_clean.fastq.gz $reads
    """
}
