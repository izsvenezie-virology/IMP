#! /usr/bin/env python

# Reads a tsv file. Columns: primer name, primer sequence
# Writes a fasta file for 5' primers (X[primer sequence])
# Writes a fasta file for 3' primers ([primer sequence reverse complement]X)

import sys

primers_tsv = sys.argv[1]

complements = {
    'A': 'T', 'T': 'A', 'U': 'A',
    'C': 'G', 'G': 'C',
    'Y': 'R', 'R': 'Y',
    'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K',
    'B': 'V', 'V': 'B',
    'D': 'H', 'H': 'D',
    'N': 'N'
}

with open(primers_tsv, 'r') as f_in:
    with open('primers_5g', 'w') as f_5first:
        with open('primers_3a', 'w') as f_3first:
            for line in f_in:
                id, sequence = line.strip().split('\t')
                f_5first.write(f'>{id}\nX{sequence}\n')
                rc_sequence = ''.join([complements[c] for c in reversed(sequence.upper())])
                f_3first.write(f'>{id}_rc\n{rc_sequence}X\n')
