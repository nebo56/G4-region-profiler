#!/usr/bin/env python3
"""
Flank BED intervals in a strand-aware way, keeping all extra columns.

"Upstream" and "downstream" are defined relative to the strand:
on '+' the interval is extended `upstream` bp to the left and
`downstream` bp to the right; on '-' the reverse.

The original position is appended as a final column (chrom:start:end)
so regions can be traced back after flanking.

Intervals whose flanked start would become negative are skipped
(a warning with the count is printed to stderr).

Usage:
    python3 02_flank_bed.py input.bed output.bed <upstream_bp> <downstream_bp>

Example (TSS -100/+100, matching the example FASTA in this repo):
    python3 02_flank_bed.py genes-Start.bed genes-Start-flanked100.bed 100 100

Author: Nejc Haberman
"""

import sys


def flank_positions(fin_fname, fout_fname, upstream, downstream):
    n_skipped = 0
    with open(fin_fname) as fin, open(fout_fname, "w") as fout:
        for line in fin:
            col = line.rstrip("\n").split("\t")
            if len(col) < 6:
                n_skipped += 1
                continue
            chrom, start, end, strand = col[0], int(col[1]), int(col[2]), col[5]

            if strand == "+":
                new_start = start - upstream
                new_end = end + downstream
            elif strand == "-":
                new_start = start - downstream
                new_end = end + upstream
            else:
                n_skipped += 1
                continue

            if new_start < 0 or new_start >= new_end:
                n_skipped += 1
                continue

            original = f"{chrom}:{start}:{end}"
            extra = "\t".join(col[3:])
            fout.write(f"{chrom}\t{new_start}\t{new_end}\t{extra}\t{original}\n")

    if n_skipped:
        sys.stderr.write(f"Skipped {n_skipped} lines (short lines, unknown "
                         f"strand, or negative coordinates after flanking)\n")


def main():
    if len(sys.argv) != 5:
        sys.stderr.write("Usage: python3 02_flank_bed.py input.bed output.bed "
                         "<upstream_bp> <downstream_bp>\n")
        sys.exit(1)
    flank_positions(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))


if __name__ == "__main__":
    main()
