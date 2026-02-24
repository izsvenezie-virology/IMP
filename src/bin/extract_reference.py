#! /usr/bin/env python

# Reads Blast output
# Counts best hit for each sequence
# Prints the name of sequences to be used for each segment

import sys
from collections import defaultdict

best_hits_file = sys.argv[1]

sequences = defaultdict(lambda: defaultdict(int))

with open(best_hits_file, "r") as f:
    for line in f:
        ref = line.split("\t", 1)[0]
        try:
            seg = ref.rsplit("|", 1)[1]
        except IndexError:
            continue
        sequences[seg][ref] += 1

for seg in sequences.keys():
    print(f">{max(sequences[seg], key=sequences[seg].get)}", end="\n")  # type: ignore
