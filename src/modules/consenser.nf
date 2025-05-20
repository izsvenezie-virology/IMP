process Consenser {
    tag "${id.sample}"

    memory '500 MB'
    time '30s'

    input:
    tuple val(id), val(parameters), path(vcf), path(coverage), path(reference)
    val degenerated

    output:
    tuple val(id), path("${id.sample}__${id.reference}.fa"), emit: consensus
    tuple val(id), path("${id.sample}__${id.reference}_*.fa"), emit: segments

    script:
    def deg_option = degenerated ? '-d' : ''
    """
    consenser --force ${deg_option} -s ${id.sample}__${id.reference}_CHROMNAME.fa \
    -a '${parameters.name}' -o ${id.sample}__${id.reference}.fa -c ${coverage} ${reference} ${vcf}
    """
}
