#! /usr/bin/env python

# Reads Blast output
# Counts best hit for each sequence
# Prints the name of sequences to be used for each segment

import sys
from collections import defaultdict


def select_subtypes(subtypes, prefix):
    max_count = max(sum(refs.values()) for refs in subtypes.values()) if subtypes else 0

    detected_subtypes = []
    total_hits = []
    best_refs = []
    for subtype, refs in sorted(subtypes.items()):
        if sum(refs.values()) >= max_count * 0.1:
            subtype_name = f"{prefix}{str(subtype)}"
            detected_subtypes.append(subtype_name)
            total_hits.append(str(sum(refs.values())))
            best_refs.append(
                f"{subtype_name}: {max(refs, key=refs.get)} ({str(max(refs.values()))})"
            )
    return detected_subtypes, total_hits, best_refs


best_hits_file = sys.argv[1]

sequences_ha = defaultdict(lambda: defaultdict(int))
sequences_na = defaultdict(lambda: defaultdict(int))

with open(best_hits_file, "r") as f:
    for line in f:
        ref = line.split("\t", 1)[0]
        id, subtype, seg = ref.split("|")
        if not subtype:
            continue
        subtype = subtype
        if seg == "HA":
            sequences_ha[subtype][ref] += 1
        elif seg == "NA":
            sequences_na[subtype][ref] += 1

ha_subtypes, ha_hits, ha_refs = select_subtypes(sequences_ha, "H")
na_subtypes, na_hits, na_refs = select_subtypes(sequences_na, "N")

subtypes = "".join(ha_subtypes + na_subtypes)
hits = ",".join(ha_hits + na_hits)
refs = ", ".join(ha_refs + na_refs)

result = f"{subtypes}\t{hits}\t{refs}"

print(result, end="\n")
