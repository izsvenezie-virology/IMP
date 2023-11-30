process DictIndex{
    tag "$meta"

    input:
        tuple val(meta), path(reference)
    output:
        tuple val(meta), path('*')

    """
    gatk CreateSequenceDictionary -R $reference
    """
}

process IndexFeatureFile{
    tag "$meta.sample"

    input:
        tuple val(meta), path(feature_file)
    output:
        tuple val(meta), path('*')

    """
    gatk IndexFeatureFile -I $feature_file
    """
}

process FixBam{
    tag "$meta.sample"

    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path("*")

    """
    gatk FixMateInformation -I ${bam} -O fixed.bam
    """
}

process CleanBam{
    tag "$meta.sample"

    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path("*")

    """
    gatk CleanSam -I ${bam} -O cleaned.bam
    """
}

process MarkDuplicates{
    tag "$meta.sample"

    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path("*.bam")

    """
    gatk MarkDuplicates -I ${bam} -O marked_duplicates.bam -M marked_duplicates.metrics
    """
}

process BaseRecalibrator{
    tag "$meta.sample"

    input:
        tuple val(meta), path(bam), path(reference), path(ref_idx), path(dict_idx), path(known_sites), path(feature_idx)
    output:
        tuple val(meta), path("*")

    """
    gatk BaseRecalibrator -R $reference -I $bam --known-sites $known_sites -O recalibration_report.txt
    """
}

process ApplyBQSR{
    tag "$meta.sample"
    
    input:
        tuple val(meta), path(bam), path(recalibration)
    output:
        tuple val(meta), path("*.bam"), path("*.bai")

    """
    gatk ApplyBQSR -I $bam -O bqsr.bam -bqsr $recalibration
    """
}
