#! /usr/bin/env python

import sys
from collections import defaultdict

import pysam

sample_id = sys.argv[1]
bam_file = sys.argv[2]


with pysam.AlignmentFile(bam_file, "rb") as bam:
    mapped_reads = defaultdict(set)
    unmapped_reads = set()
    for read in bam.fetch(until_eof=True):
        if read.is_secondary:
            continue
        if read.is_unmapped:
            unmapped_reads.add(read)
            continue
        mapped_reads[read.reference_name].add((read.query_name, read.is_read1))

reads_count = defaultdict(int)
for reference in mapped_reads:
    reads_count[reference] = len(mapped_reads[reference])
reads_count["Total"] = sum(reads_count.values())

for chrom, count in reads_count.items():
    print(f"{sample_id}\t{chrom}\tMapped reads\t{count}")
print(f"{sample_id}\t-\tUnmapped reads\t{len(unmapped_reads)}")
