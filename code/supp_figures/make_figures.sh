#!/usr/bin/env bash
# Regenerate Supplemental Figs S3 and S4, Supplemental Tables S2-S5 (TSV) and numbers.json
# from the summary tables in results/.   Usage (from the repository root): bash code/supp_figures/make_figures.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
python3 -I "$ROOT/code/supp_figures/make_supp_items.py"
python3 -I "$ROOT/code/supp_figures/make_supp_figS4.py"
cp "$ROOT"/results/supp_items/Supplemental_Fig_S{3,4}.{pdf,png} "$ROOT/figures/supplemental/"
echo "done: figures/supplemental/Supplemental_Fig_S{3,4}.{pdf,png}, results/supp_items/"
