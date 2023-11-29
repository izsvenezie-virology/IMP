process Consensus{
    input:
        tuple val(meta), path(vcf), path(reference), path(coverage)
        val(degenerated)
    output:
        tuple val(meta), path("consensus.fa"), emit: consensus
        tuple val(meta), path("*.fasta"), emit: segments

    script:
    def degenerated_opt = degenerated ? '-d' : '' // If degenerated is true the degenerated consensus is produced
    """
    consenser --force $degenerated_opt -s CHROMNAME.fasta -a ${meta.name}_CHROMNAME -o consensus.fa -c $coverage $reference $vcf
    """
}
