#!/usr/bin/env python3
"""
Convert merged G4Hunter output into a per-region score table (BED6+3).

Expects FASTA headers produced by `bedtools getfasta -s -name`, i.e.:
    >NAME::chrom:start-end(strand)
followed by G4Hunter window lines where the score is in column 5.

Because `getfasta -s` already reverse-complements minus-strand regions,
positive G4Hunter scores correspond to G4 propensity on the SENSE
(transcribed) strand and negative scores to the ANTISENSE strand.

Output columns (BED6 + 3 extra):
    1  chrom
    2  start
    3  end
    4  name (e.g. transcript ID)
    5  score  = sum of positive window scores (sense G4 propensity)
    6  strand
    7  max positive window score        (0 if none)
    8  sum of |negative| window scores  (antisense G4 propensity)
    9  max |negative| window score      (0 if none)

Regions with no windows above threshold are written with zero scores,
so the output has one line per input sequence and can be joined
directly with gene annotation.

Usage:
    python3 04_g4hunter_to_bed.py <g4hunter_merged.txt> <output.bed>

Author: Nejc Haberman
"""

import re
import sys

HEADER_RE = re.compile(
    r"^>(?P<name>.+?)::(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)"
    r"\((?P<strand>[+-.])\)"
)


class Region:
    __slots__ = ("chrom", "start", "end", "name", "strand",
                 "pos_scores", "neg_scores")

    def __init__(self, chrom, start, end, name, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.strand = strand
        self.pos_scores = []
        self.neg_scores = []


def write_region(fout, r):
    sum_pos = sum(r.pos_scores)
    max_pos = max(r.pos_scores) if r.pos_scores else 0.0
    sum_neg = sum(abs(s) for s in r.neg_scores)
    max_neg = max((abs(s) for s in r.neg_scores), default=0.0)
    fout.write("\t".join((
        r.chrom, str(r.start), str(r.end), r.name,
        f"{sum_pos:.3f}", r.strand,
        f"{max_pos:.3f}", f"{sum_neg:.3f}", f"{max_neg:.3f}")) + "\n")


def g4hunter_to_bed(input_file, output_file):
    n_regions = 0
    n_windows = 0
    region = None
    with open(input_file) as fin, open(output_file, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                if region is not None:
                    write_region(fout, region)
                m = HEADER_RE.match(line.rstrip("\n"))
                if m is None:
                    sys.stderr.write(
                        "ERROR: cannot parse header (expected "
                        "'>NAME::chrom:start-end(strand)' from "
                        "bedtools getfasta -s -name):\n" + line)
                    sys.exit(1)
                region = Region(m["chrom"], int(m["start"]), int(m["end"]),
                                m["name"], m["strand"])
                n_regions += 1
            elif region is not None:
                columns = line.rstrip("\n").split("\t")
                if len(columns) < 5:
                    continue
                try:
                    score = float(columns[4])
                except ValueError:
                    continue  # column header or malformed line
                n_windows += 1
                if score > 0:
                    region.pos_scores.append(score)
                elif score < 0:
                    region.neg_scores.append(score)
        if region is not None:
            write_region(fout, region)

    sys.stderr.write(f"Processed {n_regions} regions, "
                     f"{n_windows} G4Hunter windows\n")


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python3 04_g4hunter_to_bed.py "
                         "<g4hunter_merged.txt> <output.bed>\n")
        sys.exit(1)
    g4hunter_to_bed(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
