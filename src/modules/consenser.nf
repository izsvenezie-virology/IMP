process DegeneratedConsensus{
    tag "$id.sample"

    memory '500 MB'
    time '30s'

    publishDir "consensus", saveAs: { "${id.sample}__${id.reference}.fa" }, mode: 'copy', pattern: 'consensus.fasta'
    publishDir "consensus/chroms", saveAs: { "${id.sample}__${id.reference}_${it}" }, mode: 'copy', pattern: '*.fa'
    
    input:
        tuple val(id), val(parameters), path(vcf), path(coverage), path(reference)
    output:
        tuple val(id), path("consensus.fasta"), emit: consensus
        tuple val(id), path("*.fa"), emit: segments

    """
    consenser --force -d -s CHROMNAME.fa -a '${parameters.name}' -o consensus.fasta -c $coverage $reference $vcf
    """
}

process NonDegeneratedConsensus{
    tag "$id.sample"

    publishDir "consensus_no_degenerations", saveAs: { "${id.sample}__${id.reference}.fa" }, mode: 'copy', pattern: 'consensus.fasta'
    publishDir "consensus_no_degenerations/chroms", saveAs: { "${id.sample}__${id.reference}_${it}" }, mode: 'copy', pattern: '*.fa'
    
    memory '500 MB'
    time '30s'

    input:
        tuple val(id), val(parameters), path(vcf), path(coverage), path(reference)
    output:
        tuple val(id), path("consensus.fasta"), emit: consensus
        tuple val(id), path("*.fa"), emit: segments

    """
    consenser --force -s CHROMNAME.fa -a '${parameters.name}' -o consensus.fasta -c $coverage $reference $vcf
    """
}
