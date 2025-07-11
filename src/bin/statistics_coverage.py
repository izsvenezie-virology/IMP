#! /usr/bin/env python

from collections import defaultdict
from statistics import mean
import sys

sample_id = sys.argv[1]
horizontal_threshold = int(sys.argv[2])
coverage_file = sys.argv[3]

with open(coverage_file, 'r') as f:
    coverages = defaultdict(list)
    for line in f:
        chrom, _, coverage = line.strip().split('\t')
        coverages['Total'].append(int(coverage))
        coverages[chrom].append(int(coverage))

for chrom, covs in coverages.items():
    mean_cov = mean(covs)
    print(f'{sample_id}\t{chrom}\tMean coverage\t{mean_cov:.2f}')
    
    uniformity_threshold = mean_cov * 0.2
    uniformity = sum(1 for cov in covs if cov >= uniformity_threshold)
    print(f'{sample_id}\t{chrom}\tCoverage uniformity\t{uniformity/len(covs):.2f}')

    horizontal_coverage = sum(1 for cov in covs if cov >= horizontal_threshold)    
    print(f'{sample_id}\t{chrom}\tHorizontal coverage\t{horizontal_coverage/len(covs):.2f}')
