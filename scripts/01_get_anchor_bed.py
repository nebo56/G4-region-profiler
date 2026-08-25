#!/usr/bin/env python3
"""
Extract single-nucleotide start (TSS) and end (TES) positions from a BED6+ file.

For each interval, writes two BED files:
  <prefix>-Start.bed : the strand-aware 5' position (e.g. TSS for genes/transcripts)
  <prefix>-End.bed   : the strand-aware 3' position (TES)

Coordinates follow BED convention (0-based, half-open). For minus-strand
intervals the 5'-most base is [end-1, end) and the 3'-most base is
[start, start+1).

Usage:
    python3 01_get_anchor_bed.py input.bed output_prefix

Author: Nejc Haberman
"""

import sys


def get_start_end(fin_fname, fout_prefix):
    n_skipped = 0
    with open(fin_fname) as fin, \
            open(fout_prefix + "-Start.bed", "w") as fout_start, \
            open(fout_prefix + "-End.bed", "w") as fout_end:
        for line in fin:
            col = line.rstrip("\n").split("\t")
            if len(col) < 6:
                n_skipped += 1
                continue
            chrom, start, end = col[0], int(col[1]), int(col[2])
            name, score, strand = col[3], col[4], col[5]

            if strand == "+":
                five_prime = start          # TSS base: [start, start+1)
                three_prime = end - 1       # TES base: [end-1, end)
            elif strand == "-":
                five_prime = end - 1        # TSS base: [end-1, end)
                three_prime = start         # TES base: [start, start+1)
            else:
                n_skipped += 1
                continue

            fout_start.write("\t".join(map(str, (
                chrom, five_prime, five_prime + 1, name, score, strand))) + "\n")
            fout_end.write("\t".join(map(str, (
                chrom, three_prime, three_prime + 1, name, score, strand))) + "\n")

    if n_skipped:
        sys.stderr.write(f"Skipped {n_skipped} lines (fewer than 6 columns "
                         f"or unknown strand)\n")


def main():
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage: python3 01_get_anchor_bed.py input.bed output_prefix\n")
        sys.exit(1)
    get_start_end(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
