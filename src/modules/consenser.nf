process DegeneratedConsensus{
    tag "$meta.sample"

    publishDir "consensus", saveAs: { "${meta.sample}__${meta.reference}.fa" }, mode: 'copy', pattern: 'consensus.fasta'
    publishDir "consensus/chroms", saveAs: { "${meta.sample}__${meta.reference}_${it}" }, mode: 'copy', pattern: '*.fa'
    
    input:
        tuple val(meta), path(vcf), path(reference), path(coverage)
    output:
        tuple val(meta), path("consensus.fasta"), emit: consensus
        tuple val(meta), path("*.fa"), emit: segments

    """
    consenser --force -d -s CHROMNAME.fa -a ${meta.name}_CHROMNAME -o consensus.fasta -c $coverage $reference $vcf
    """
}

process NonDegeneratedConsensus{
    tag "$meta.sample"

    publishDir "consensus_no_degenerations", saveAs: { "${meta.sample}__${meta.reference}.fa" }, mode: 'copy', pattern: 'consensus.fasta'
    publishDir "consensus_no_degenerations/chroms", saveAs: { "${meta.sample}__${meta.reference}_${it}" }, mode: 'copy', pattern: '*.fa'
    
    input:
        tuple val(meta), path(vcf), path(reference), path(coverage)
    output:
        tuple val(meta), path("consensus.fasta"), emit: consensus
        tuple val(meta), path("*.fa"), emit: segments

    """
    consenser --force -s CHROMNAME.fa -a ${meta.name}_CHROMNAME -o consensus.fasta -c $coverage $reference $vcf
    """
}
