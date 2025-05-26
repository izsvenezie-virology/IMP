process Consenser {
    tag "${id}"

    memory '500 MB'
    time '30s'

    input:
    tuple val(id), val(parameters), path(vcf), path(coverage), path(reference)
    val degenerated

    output:
    tuple val(id), path("${id}.fa"), emit: consensus
    tuple val(id), path("${id}_*.fa"), emit: segments

    script:
    def deg_option = degenerated ? '-d' : ''
    """
    consenser --force ${deg_option} -s ${id}_CHROMNAME.fa \
    -a '${parameters.name}' --min-cov ${parameters.minimum_coverage} -o ${id}.fa -c ${coverage} ${reference} ${vcf}
    """
}
