process BlastN{
    tag "$meta"
    label 'multiThread'

    input:
        tuple val(meta), path(reads_fasta)
        tuple val(db_name), path(db)
    output:
        tuple val(meta), path('*')
    
    """
    blastn -db ${db_name} -query $reads_fasta \
        -out best_hits.tsv -outfmt '6 sseqid qseqid evalue bitscore score' \
        -task megablast -num_threads $task.cpus \
        -evalue 1e-50 -max_target_seqs 100
    """
}

process MakeBlastDb{
    tag "$references_fasta.baseName"

    input:
        path(references_fasta)
        val (enabled)
    output:
        tuple val('blastdb'), path('*')
    when:
        enabled != false

    """
    makeblastdb -in $references_fasta -dbtype nucl -out blastdb
    """
}
