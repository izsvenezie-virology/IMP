process Consenser {
    tag "${id}"

    memory '500 MB'
    time '1h'

    input:
    tuple val(id), val(parameters), path(vcf), path(coverage), path(reference)

    output:
    tuple val(id), path("${id}_consensus.fa"), emit: degenerated
    tuple val(id), path("${id}_consensus_non_degenerated.fa"), emit: non_degenerated

    script:
    """
    consenser -H '${parameters.name}' -m ${parameters.minimum_coverage} -o ${id}_consensus.fa -c ${coverage} ${reference} ${vcf}
    consenser -H '${parameters.name}' -m ${parameters.minimum_coverage} -s 1 -o ${id}_consensus_non_degenerated.fa -c ${coverage} ${reference} ${vcf}
    """
}
