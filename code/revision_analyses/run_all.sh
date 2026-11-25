#!/bin/bash
# Order in which the analyses were run (FRX, 2026-10-15/16). Heavy steps ran as background pools; see eric/logs. Python = homl4 env with -I.
set -e; PY=/home/nelazzabi/anaconda3/envs/homl4/bin/python; cd ~/xcell-xbulk-revision
Rscript --vanilla eric/R/e0_export_saved_tables.R                  # Eric's saved fisher tables -> TSV
for l in sc sbj bk; do Rscript --vanilla eric/R/e1_export_edges.R $l 10; done   # edge sets from his saved matrices
$PY -I eric/baseline/e2_baseline.py; $PY -I eric/baseline/e4_bulk_baseline.py; $PY -I eric/make_A5_A6_A9.py
Rscript --vanilla eric/R/e3_export_counts.R; Rscript --vanilla eric/R/e5_export_bulk_expr.R   # raw counts / bulk expression for reruns
OPENBLAS_NUM_THREADS=4 $PY -I lib/eric_pass.py 10            # xCell variants (std/A7/A8/raw-r), validated vs saved matrices
OPENBLAS_NUM_THREADS=2 $PY -I lib/eric_xsbj.py 4              # xSubject variants (A4)
for a in A7 A8 A6raw A4; do $PY -I lib/eric_compare.py $a 10; done
$PY -I A3_expression_strata/a3_eric.py; $PY -I A3_expression_strata/a3_high_compare.py
$PY -I lib/eric_a10.py groups 6; $PY -I lib/eric_a10.py bulk; $PY -I lib/eric_a10.py compare 10
for s in S1 S2; do $PY -I lib/eric_a1.py $s 14 50; $PY -I lib/eric_a1_compare.py $s 16; done   # A1, seed 20261015, 50 draws
$PY -I lib/eric_a2.py xsubject 6 0 50; $PY -I lib/eric_a2.py xbulk 2 10 20; $PY -I lib/eric_a2.py full 14 1 10; $PY -I lib/eric_a2.py S1 6 10 20; $PY -I lib/eric_a2.py S2 10 3 10
$PY -I lib/eric_a2_filtered.py
