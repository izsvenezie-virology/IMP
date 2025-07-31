#! /usr/bin/env python

import sys

import pysam

sample_id = sys.argv[1]
reads_type = sys.argv[2]
fastq_r1_file = sys.argv[3]
fastq_r2_file = sys.argv[4]


def parse_fastq(fastq: str, r: str):
    quality = 0
    quality_count = 0
    reads_count = 0
    with pysam.FastxFile(fastq) as fq:
        for read in fq:
            quality += sum(read.get_quality_array())
            quality_count += len(read.sequence)
            reads_count += 1
    print(f"{sample_id}\t{r}\tReads {reads_type} count\t{reads_count}")
    print(
        f"{sample_id}\t{r}\tReads {reads_type} mean quality\t{quality / quality_count:.2f}"
    )


parse_fastq(fastq_r1_file, "R1")
parse_fastq(fastq_r2_file, "R2")
