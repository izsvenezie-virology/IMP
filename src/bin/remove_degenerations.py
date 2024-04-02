#! /usr/bin/env python

import sys

deg_dict = {
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G']
}

input_fasta = sys.argv[1]

with open(input_fasta, 'r') as f:
    for line in f:
        if line.startswith('>'):
            print(line.strip(), end='\n')
            continue
        undeg_seq = list(line.strip().upper())
        for pos, nucl in enumerate(undeg_seq):
            if nucl not in 'ACGTRYSWKMBDHV':
                raise ValueError(f'Found invalid character in sequence: {nucl}')
            if nucl in deg_dict:
                near = undeg_seq[max(pos - 5, 0):min(pos + 5, len(undeg_seq)-1)]
                near_count = dict.fromkeys(deg_dict[nucl],0)
                for key in near_count.keys():
                    near_count[key] = near.count(key)
                undeg_seq[pos] = min(near_count, key=near_count.get)
        print(''.join(undeg_seq), end='\n')
