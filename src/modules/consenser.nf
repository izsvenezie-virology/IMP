process DegeneratedConsensus {
    tag "${id.sample}"

    memory '500 MB'
    time '30s'

    publishDir "consensus", mode: 'copy', pattern: "${id.sample}__${id.reference}.fa"
    publishDir "consensus/chroms", mode: 'copy', pattern: "${id.sample}__${id.reference}_*.fa"

    input:
    tuple val(id), val(parameters), path(vcf), path(coverage), path(reference)

    output:
    tuple val(id), path("${id.sample}__${id.reference}.fa"), emit: consensus
    tuple val(id), path("${id.sample}__${id.reference}_*.fa"), emit: segments

    script:
    """
    consenser --force -d -s ${id.sample}__${id.reference}_CHROMNAME.fa -a '${parameters.name}' -o ${id.sample}__${id.reference}.fa -c ${coverage} ${reference} ${vcf}
    """
}

process NonDegeneratedConsensus {
    tag "${id.sample}"

    publishDir "consensus_no_degenerations", mode: 'copy', pattern: "${id.sample}__${id.reference}.fa"
    publishDir "consensus_no_degenerations/chroms", mode: 'copy', pattern: "${id.sample}__${id.reference}_*.fa"

    memory '500 MB'
    time '30s'

    input:
    tuple val(id), val(parameters), path(vcf), path(coverage), path(reference)

    output:
    tuple val(id), path("${id.sample}__${id.reference}.fa"), emit: consensus
    tuple val(id), path("${id.sample}__${id.reference}_*.fa"), emit: segments

    script:
    """
    consenser --force -s ${id.sample}__${id.reference}_CHROMNAME.fa -a '${parameters.name}' -o ${id.sample}__${id.reference}.fa -c ${coverage} ${reference} ${vcf}
    """
}
