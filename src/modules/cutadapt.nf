process Cutadapt {
    tag "$meta.sample"
    label 'multiThread'
    publishDir "clean_reads", saveAs: { "${meta.sample}_R1_clean.fastq.gz" }, mode: 'symlink', pattern: '*_R1_*'
    publishDir "clean_reads", saveAs: { "${meta.sample}_R2_clean.fastq.gz" }, mode: 'symlink', pattern: '*_R2_*'

    memory '4 GB'
    time '1m'

    input:
        tuple val(meta), path(reads), path(primers)
        path(adapters)
    output:
        tuple val(meta), path('*_clean.fastq.gz')

    script:
    def phred_threshold = meta.phred_threshold ?: 20
    def min_len = meta.min_len ?: 80
    def remove_primers_opt = ( primers.name == file(params.null_file).name ) ? '' : """-a file:primers_3a.fa -A file:primers_3a.fa \
                                                                                       -g file:primers_5g.fa -G file:primers_5g.fa \
                                                                                       --times=3"""
    """
    cutadapt --interleaved -j $task.cpus \
        --action=trim --pair-filter=any \
        -a file:$adapters -A file:$adapters \
        $remove_primers_opt \
        -q $phred_threshold -m $min_len \
        -o reads_R1_clean.fastq.gz -p reads_R2_clean.fastq.gz \
        $reads
    """
}
