#! /usr/bin/env python

# This script apply following changes to reference:
# - Remove degenerations
# - Remove gaps
# - Remove spaces and forbidden characters from FASTA headers
# - Writes sequences on one line
# - Set new line to LF
# - Sequences to uppercase
# - Raise error if non-standard IUPAC nucleotides are found in sequence

import sys

# fmt: off
deg_dict = {
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G']
}
# fmt: on


def remove_degenerations(sequence):
    seq = list(sequence)
    for pos, nucl in enumerate(seq):
        if nucl not in "ACGTRYSWKMBDHVN-":
            raise ValueError(f"Found invalid character in sequence: {nucl}")
        if nucl == "-":
            seq[pos] = ""
        if nucl in deg_dict:
            near = seq[max(pos - 5, 0) : min(pos + 5, len(seq) - 1)]
            near_count = dict.fromkeys(deg_dict[nucl], 0)
            for key in near_count.keys():
                near_count[key] = near.count(key)
            seq[pos] = min(near_count, key=near_count.get)  # type: ignore
    return "".join(seq)


def format_header(header):
    h = header[1:].strip()
    h = h.replace(" ", "_")
    h = h.replace("/", "_")
    h = h.replace("\\", "_")
    h = h.replace("<", "_")
    h = h.replace(">", "_")
    h = h.replace(".", "_")
    h = h.replace(",", "_")
    h = h.replace(":", "_")
    h = h.replace(";", "_")
    h = h.replace('"', "_")
    h = h.replace("'", "_")
    h = h.replace("|", "_")
    h = h.replace("?", "_")
    h = h.replace("*", "_")
    return f">{h}"


input_fasta = sys.argv[1]

sequence = ""

with open(input_fasta, "r") as f:
    for line in f:
        if line.startswith(">"):
            if sequence:
                print(remove_degenerations(sequence), end="\n")
                sequence = ""
            print(format_header(line), end="\n")
            continue
        sequence += line.strip().upper()
print(remove_degenerations(sequence), end="\n")
