# G4-region-profiler

**Profile G-quadruplex (G4) propensity across any set of genomic regions:
per-region scores and a per-nucleotide meta-profile around an anchor point
(TSS, 3′UTR start, protein binding sites, ...), built on
[G4Hunter](https://doi.org/10.1093/nar/gkw006).**

Take a BED file of regions, extract strand-aware anchor points and flanks,
run G4Hunter (Bedrat, Lacroix & Mergny, 2016) at one or more
window/threshold settings, and summarise the predictions two ways:

- **per-region scores** (BED6+3): summed and maximum G4 score on the
  sense and antisense strand of every region — ready to join with gene
  annotation or overlap with CLIP/ChIP peaks;
- **meta-profile**: summed G4 score at each nucleotide position relative
  to the anchor, aggregated across all regions.

Developed on mouse (mm10, GENCODE) TSS regions, but nothing is
genome- or feature-specific: any BED6 annotation + genome FASTA works.

```
region BED6
   │  01_get_anchor_bed.py               → single-nt 5'/3' anchor positions
   ▼
anchor BED
   │  02_flank_bed.py (e.g. 100 100)     → strand-aware flanks
   ▼
flanked BED
   │  bedtools getfasta -s -name         → sequences (see below)
   ▼
FASTA
   │  03_run_g4hunter.sh                 → G4Hunter, several W/S settings
   ▼
G4Hunter output (merged .txt)
   ├─ 04_g4hunter_to_bed.py              → per-region scores (BED6+3)
   ├─ 05_g4_metaprofile.py               → per-nucleotide anchor profile
   └─ 06_annotate_and_plot.R             → gene symbols + figures
```

Steps are independent: if you already have a FASTA of regions (or your
regions are intervals rather than anchor flanks), start directly at
step 3.

## Dependencies

- Python ≥ 3.8 (pipeline scripts; standard library only)
- Python 2.7 (required by G4Hunter itself)
- [G4Hunter](https://github.com/AnneCambon/G4-hunter)
- [bedtools](https://bedtools.readthedocs.io)
- R with `ggplot2` and `cowplot`

## Usage (example: mouse TSS ±100 bp)

### 1. Anchor positions and flanking regions

```bash
python3 scripts/01_get_anchor_bed.py mm10.genes.bed mm10.genes
# writes mm10.genes-Start.bed (5' ends / TSS) and mm10.genes-End.bed (3' ends)

python3 scripts/02_flank_bed.py mm10.genes-Start.bed \
    mm10.genes-Start-flanked100.bed 100 100
```

Both scripts are strand-aware: for minus-strand features the 5' end is
taken from the BED `end` coordinate and "upstream" extends to the
right. Any single-nucleotide anchor works the same way — 3′UTR starts,
polyA sites, CLIP crosslink positions.

### 2. Extract sequences

```bash
bedtools getfasta -s -name \
    -fi mm10.fa \
    -bed mm10.genes-Start-flanked100.bed \
    -fo mm10.genes-Start-flanked100.fa
```

**Both `-s` and `-name` are required.** `-s` reverse-complements
minus-strand regions (so positive G4Hunter scores always refer to the
sense strand), and `-name` produces headers of the form
`>NAME::chrom:start-end(strand)` that the downstream parsers rely on.

### 3. Run G4Hunter

```bash
./scripts/03_run_g4hunter.sh mm10.genes-Start-flanked100.fa \
    ~/Programs/G4Hunter/G4Hunter.py
```

Runs the window/threshold settings recommended in the G4Hunter paper
(FDR 10%: W30/S1.4 or W20/S1.7; FDR 5%: W30/S1.4 or W25/S1.6) plus one
exploratory low-stringency run (W25/S0.8). Edit the `SETTINGS` array
to change them. Concatenate G4Hunter's per-sequence output into one
merged `.txt` before the next step if your G4Hunter version splits it.

### 4. Per-region scores

```bash
python3 scripts/04_g4hunter_to_bed.py \
    G4Hunter.RESULTS_S1.6_W25/merged.txt \
    mm10.TSS_flank100.W25_S1.6.scores.bed
```

Output is BED6+3, one line per input sequence (regions without any
window above threshold get zeros):

| col | content |
|-----|---------|
| 1–3 | chrom, start, end (genomic, from the FASTA header) |
| 4   | name (transcript/gene ID) |
| 5   | **sum of positive window scores** (sense-strand G4 propensity) |
| 6   | strand |
| 7   | max positive window score |
| 8   | sum of \|negative\| window scores (antisense G4 propensity) |
| 9   | max \|negative\| window score |

### 5. Anchor meta-profile

```bash
python3 scripts/05_g4_metaprofile.py \
    G4Hunter.RESULTS_S1.6_W25/merged.txt \
    mm10.TSS_flank100.W25_S1.6.profile.tab \
    --upstream 100
```

Each window's score is added to every nucleotide it covers, using the
window coordinates *relative to the extracted sequence*, then shifted
by `--upstream` so that position 0 = the anchor (here, the TSS).

### 6. Annotate and plot

```bash
Rscript scripts/06_annotate_and_plot.R \
    mm10.TSS_flank100.W25_S1.6.scores.bed \
    mm10.TSS_flank100.W25_S1.6.profile.tab \
    mm10.genes.tab \
    mm10.TSS_flank100.W25_S1.6
```

`mm10.genes.tab` is a headerless table with columns
`gene.id chr strand start end gene.id2 gene.symbol`.
Writes the annotated score table plus three PDFs (anchor meta-profile,
sum- and max-score distributions).

## Testing

A small mock G4Hunter output is provided:

```bash
python3 scripts/04_g4hunter_to_bed.py example/example_g4hunter_output.txt /tmp/out.bed
python3 scripts/05_g4_metaprofile.py example/example_g4hunter_output.txt /tmp/profile.tab --upstream 100
```

## Citation

If you use this pipeline, please cite G4Hunter:

Bedrat A, Lacroix L, Mergny J-L. *Re-evaluation of G-quadruplex
propensity with G4Hunter.* Nucleic Acids Research 44(4):1746–1759
(2016). doi:10.1093/nar/gkw006



