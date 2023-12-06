process Consensus{
    tag "$meta.sample"

    publishDir "consensus", saveAs: { "${meta.sample}__${meta.reference}.fa" }, mode: 'copy', pattern: 'consensus.fasta'
    publishDir "consensus/chroms", saveAs: { "${meta.sample}__${meta.reference}_${it}" }, mode: 'copy', pattern: '*.fa'
    
    input:
        tuple val(meta), path(vcf), path(reference), path(coverage)
        val(degenerated)
    output:
        tuple val(meta), path("consensus.fasta"), emit: consensus
        tuple val(meta), path("*.fa"), emit: segments

    script:
    def degenerated_opt = degenerated ? '-d' : '' // If degenerated is true the degenerated consensus is produced
    """
    consenser --force $degenerated_opt -s CHROMNAME.fa -a ${meta.name}_CHROMNAME -o consensus.fasta -c $coverage $reference $vcf
    """
}
