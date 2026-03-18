#! /usr/bin/env python

import sys
from collections import defaultdict

import pysam

sample_id = sys.argv[1]
minimum_coverage = int(sys.argv[2])
vcf_file = sys.argv[3]

with pysam.VariantFile(vcf_file, "r") as vcf:
    indel_count = 0
    frameshifts_count = 0
    muts_consensus = 0
    muts_degeneration = defaultdict(int)

    for variant in vcf.fetch():
        info = variant.info
        if info["DP"] < minimum_coverage:
            continue

        if info["INDEL"]:
            if info["AF"] > 0.5:
                indel_count += 1
                if (len(variant.alleles[0]) + len(variant.alleles[1]) - 2) % 3:  # type: ignore
                    frameshifts_count += 1
            continue

        if info["AF"] > 0.75:
            muts_consensus += 1
            continue
        if info["AF"] > 0.25:
            muts_degeneration[variant.chrom] += 1
            continue

print(f"{sample_id}\tTotal\tIndels\t{indel_count}")
print(f"{sample_id}\tTotal\tFrameshifts\t{frameshifts_count}")

print(f"{sample_id}\tTotal\tConsensus mutations\t{muts_consensus}")
for chrom, count in muts_degeneration.items():
    print(f"{sample_id}\t{chrom}\tDegenerated mutations\t{count}")
