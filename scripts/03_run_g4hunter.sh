#!/bin/bash
# ------------------------------------------------------------------
# Run G4Hunter (Bedrat et al. 2016, NAR; PMC4770238) on a FASTA of
# extracted regions (e.g. TSS flanks), for one or more window/threshold settings.
#
# Recommended settings from the G4Hunter paper:
#   FDR 10%: W30/S1.4 or W20/S1.7
#   FDR  5%: W30/S1.4 or W25/S1.6
# An exploratory low-stringency run (W25/S0.8) can be added to
# inspect the full score distribution.
#
# G4Hunter itself requires Python 2:
#   https://github.com/AnneCambon/G4-hunter (or the original release)
#
# Usage:
#   ./03_run_g4hunter.sh <input.fa> [path/to/G4Hunter.py]
# ------------------------------------------------------------------
set -euo pipefail

FASTA="${1:?Usage: $0 <input.fa> [path/to/G4Hunter.py]}"
G4HUNTER="${2:-$HOME/Programs/G4Hunter/G4Hunter.py}"

# window_size:threshold pairs to run
SETTINGS=("30:1.4" "25:1.6" "20:1.7" "25:0.8")

for setting in "${SETTINGS[@]}"; do
    W="${setting%%:*}"
    S="${setting##*:}"
    OUTDIR="G4Hunter.RESULTS_S${S}_W${W}"
    mkdir -p "${OUTDIR}"
    echo ">>> G4Hunter: window=${W} threshold=${S} -> ${OUTDIR}"
    python2 "${G4HUNTER}" -i "${FASTA}" -o "${OUTDIR}" -w "${W}" -s "${S}"
done
