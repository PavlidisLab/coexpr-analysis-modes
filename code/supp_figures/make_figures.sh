#!/usr/bin/env bash
# Regenerate Supplementary Figures 5 and 6, Supplementary Tables 2-5 (TSV) and numbers.json
# from the summary tables in results/.   Usage (from the repository root): bash code/supp_figures/make_figures.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
python3 -I "$ROOT/code/supp_figures/make_supp_items.py"
python3 -I "$ROOT/code/supp_figures/make_supp_fig6.py"
cp "$ROOT"/results/supp_items/Supplementary_Figure_{5,6}.{pdf,png} "$ROOT/figures/supplementary/"
echo "done: figures/supplementary/Supplementary_Figure_{5,6}.{pdf,png}, results/supp_items/"
