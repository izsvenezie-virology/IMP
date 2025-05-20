process BlastN {
    tag "${id}"
    label 'multiThread'

    memory '10 GB'
    time '48h'

    input:
    tuple val(id), path(reads_fasta)
    tuple val(db_name), path(db)

    output:
    tuple val(id), path('best_hits.tsv')

    script:
    """
    blastn -db ${db_name} -query ${reads_fasta} \
        -out best_hits.tsv -outfmt '6 sseqid qseqid evalue bitscore score' \
        -task megablast -num_threads ${task.cpus} \
        -evalue 1e-50 -max_target_seqs 100
    """
}

process MakeBlastDb {
    tag "${references_fasta.baseName}"

    memory '1 GB'
    time '5m'

    input:
    path references_fasta
    val enabled

    output:
    tuple val('blastdb'), path('blastdb*')

    when:
    enabled != false

    script:
    """
    makeblastdb -in ${references_fasta} -dbtype nucl -out blastdb
    """
}
