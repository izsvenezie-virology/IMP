#! /usr/bin/env python

import sys
from collections import defaultdict

sample = sys.argv[1]
coverage_file = sys.argv[2]

coverages = defaultdict(int)

with open(coverage_file, "r") as f:
    for line in f:
        chrom, _, coverage = line.split("\t")
        coverages[chrom] += int(coverage)

main_chrom = max(coverages.values())

for chrom, coverage in coverages.items():
    if coverage / main_chrom >= 0.05:
        print(f"{sample}\t{chrom}\t{coverage}", end="\n")
