#!/usr/bin/env python3
"""
Build a per-nucleotide G4 meta-profile around a region anchor (e.g. TSS) from merged
G4Hunter output.

For every predicted window, its score is added to each nucleotide the
window covers, using the window Start/End columns (which are relative
to the extracted sequence). Positions are then re-centred on the anchor point
using --upstream, so x = 0 is the anchor (e.g. the TSS).

This replaces the earlier approach of aggregating scores by genomic
window centre, which mixed absolute genomic coordinates across genes
and did not produce a meaningful profile.

Output (tab-separated, with header):
    position           relative to anchor (negative = upstream)
    sense_score_sum    summed positive scores per position
    antisense_score_sum summed |negative| scores per position
    n_regions          number of sequences contributing windows there

Usage:
    python3 05_g4_metaprofile.py <g4hunter_merged.txt> <output.tab> --upstream 100

Author: Nejc Haberman
"""

import argparse
import sys
from collections import defaultdict


def metaprofile(input_file, output_file, upstream):
    sense = defaultdict(float)
    antisense = defaultdict(float)
    n_regions = defaultdict(set)
    region_id = -1

    with open(input_file) as fin:
        for line in fin:
            if line.startswith(">"):
                region_id += 1
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 5:
                continue
            try:
                w_start = int(columns[0])
                w_end = int(columns[1])
                score = float(columns[4])
            except ValueError:
                continue  # column header line
            for pos in range(w_start, w_end):
                rel = pos - upstream
                if score > 0:
                    sense[rel] += score
                else:
                    antisense[rel] += abs(score)
                n_regions[rel].add(region_id)

    if not sense and not antisense:
        sys.stderr.write("No G4Hunter windows found - check input format\n")
        sys.exit(1)

    positions = sorted(set(sense) | set(antisense))
    with open(output_file, "w") as fout:
        fout.write("position\tsense_score_sum\tantisense_score_sum\tn_regions\n")
        for pos in range(positions[0], positions[-1] + 1):
            fout.write(f"{pos}\t{sense.get(pos, 0.0):.3f}\t"
                       f"{antisense.get(pos, 0.0):.3f}\t"
                       f"{len(n_regions.get(pos, ()))}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Per-nucleotide G4 meta-profile around the region anchor")
    parser.add_argument("input", help="merged G4Hunter output (.txt)")
    parser.add_argument("output", help="output profile table (.tab)")
    parser.add_argument("--upstream", type=int, default=100,
                        help="bp upstream of the anchor in the extracted "
                             "sequences (default: 100)")
    args = parser.parse_args()
    metaprofile(args.input, args.output, args.upstream)


if __name__ == "__main__":
    main()
