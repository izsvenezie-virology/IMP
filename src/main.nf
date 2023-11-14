#! /usr/bin/env nextflow


process FastQC {
    cpus 2

    publishDir "reads_quality/$type", saveAs: { fileName -> fileName.contains("_R1_") ? "${meta.sample}_R1_${type}.html" : "${meta.sample}_R2_${type}.html" } 

    input:
        tuple val(meta), path(reads)
        val(type) // raw/clean before/after reads clean up
        
    output:
        tuple val(meta), path('*html')

    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}


workflow {
    Channel.fromFilePairs("raw_reads/*_R{1,2}*fastq.gz")
    | map{row -> [row[0].split("_S")[0], row[1]]}  
    | set { raw_reads }

    Channel.fromPath('samplesheet.csv')
    | splitCsv( header:true, sep:'\t' )
    | map { row -> [row.Sample, [sample:row.Sample, name:row.Name]] }
    | join( raw_reads )
    | map { row -> [row[1], row[2]]}
    | set { samples }

    FastQC(samples, 'raw')

}
