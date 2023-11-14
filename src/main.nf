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

process Cutadapt {
    publishDir "clean_reads", saveAs:  { fileName -> fileName.startsWith("R1_") ? "${meta.sample}_R1_clean.fastq.gz" : "${meta.sample}_R2_clean.fastq.gz" }

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path('R?_clean.fastq.gz')

    """
    cutadapt --interleaved --action=trim --pair-filter=any -q 20 -m 80 -o R1_clean.fastq.gz -p R2_clean.fastq.gz $reads
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
    Cutadapt(samples)
}
