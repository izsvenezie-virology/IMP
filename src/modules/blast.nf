process BlastN{
    tag "$meta"
    
    input:
        tuple val(meta), path(reads_fasta)
        path(db)
    output:
        tuple val(meta), path('*')
    
    """
    blastn -db ${db[0].baseName} -query $reads_fasta -out best_hits.tsv -num_threads $task.cpus \
        -outfmt '6 sseqid qseqid evalue bitscore score ' \
        -evalue 1e-50 \
        -max_target_seqs 100 
    """
}
