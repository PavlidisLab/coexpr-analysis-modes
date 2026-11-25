#!/usr/bin/env python3
"""Build every number, supplemental table and Supplemental Fig S3 for the
Genome Research revision (GENOME/2026/282552) from the FRX result files.

Input : ../../results/  (summary tables of the revision analyses; see code/revision_analyses)
Output: out/numbers.json, out/supp_table_{2,3,4,5}.tsv, out/Supplemental_Fig_S3.{pdf,png},
        out/manifest.json (sha256 + timestamps of every input; data_version)

No number is typed by hand: the manuscript, supplement and response-letter text
are filled from out/numbers.json by the other build scripts.
Run: python3 -I make_supp_items.py
"""
import csv, json, hashlib, math, os, statistics as st, sys, datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
FRX = REPO / "results"
OUT = REPO / "results" / "supp_items"  # numbers.json, Supplemental Tables S2-S5, Supplemental Fig S3, manifest.json
OUT.mkdir(exist_ok=True)
USED = {}  # input path -> sha256


def path(rel):
    p = FRX / rel
    if rel not in USED:
        USED[rel] = hashlib.sha256(p.read_bytes()).hexdigest()
    return p


def rows(rel):
    with open(path(rel), newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fnum(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def summ(ors):
    v = [o for o in ors if not math.isnan(o) and not math.isinf(o)]
    return {"n": len(ors), "n_defined": len(v), "mean": st.mean(v), "median": st.median(v),
            "min": min(v), "max": max(v), "n_le_1": sum(o <= 1 for o in v)}


def pick(rs, **kw):
    out = []
    for r in rs:
        if all(str(r.get(k)) == str(v) for k, v in kw.items()):
            out.append(r)
    return out


def frac_rows(rs, frac, **kw):
    return [r for r in pick(rs, **kw) if abs(fnum(r["frac"]) - frac) < 1e-9]


N = {}  # numbers.json

# ---------------------------------------------------------------- A6 thresholds (baseline + R2)
a6 = rows("A6_thresholds/results.tsv")
FRACS = [0.005, 0.01, 0.05, 0.1]
LV = {"sc": "xCell", "sbj": "xSubject"}
t2 = []
for an, ref, test, label in [("replicability", "sc", "sc", "xCell replicability"),
                             ("replicability", "sbj", "sbj", "xSubject replicability"),
                             ("preservation", "sc", "sbj", "xCell to xSubject preservation")]:
    for fr in FRACS:
        r = frac_rows(a6, fr, analysis=an, reference_level=ref, test_level=test)
        s = summ([fnum(x["odds_ratio"]) for x in r])
        s["mean_F1"] = st.mean(fnum(x["F1"]) for x in r)
        s["mean_jaccard"] = st.mean(fnum(x["jaccard"]) for x in r)
        key = f"a6_{an}_{ref}_{test}_{fr}"
        N[key] = s
        t2.append([label, f"{fr*100:g}%", s["n"], s["mean"], s["median"], s["min"], s["max"], s["n_le_1"],
                   s["mean_F1"], s["mean_jaccard"]])
N["baseline_repl_sc"] = N["a6_replicability_sc_sc_0.01"]
N["baseline_repl_sbj"] = N["a6_replicability_sbj_sbj_0.01"]
N["baseline_pres"] = N["a6_preservation_sc_sbj_0.01"]

# pairs for panel C (both directions of each unordered pair, 1%)
pairsC = {}
for lvl in ("sc", "sbj"):
    d = {}
    for r in frac_rows(a6, 0.01, analysis="replicability", reference_level=lvl, test_level=lvl):
        d[(r["cell_type"], r["reference_dataset"], r["test_dataset"])] = fnum(r["odds_ratio"])
    pts, logr = [], []
    for (ct, a, b), o in d.items():
        if a < b and (ct, b, a) in d:
            o2 = d[(ct, b, a)]
            pts.append((o, o2))
            logr.append(abs(math.log(o) - math.log(o2)))
    pairsC[lvl] = pts
    N[f"asym_{lvl}"] = {"n_pairs": len(pts), "median_abs_log_ratio": st.median(logr), "max_abs_log_ratio": max(logr)}
allr = [abs(math.log(a) - math.log(b)) for l in pairsC for a, b in pairsC[l]]
N["asym_all"] = {"n_ordered": 2 * len(allr), "median_abs_log_ratio": st.median(allr), "max_abs_log_ratio": max(allr)}

# ---------------------------------------------------------------- A1 matched sample (R1-M1)
for sch in ("S1", "S2"):
    for r in rows(f"A1_matched_sample/summary_{sch}.tsv"):
        N[f"a1_{sch}_{r['analysis']}_{float(r['frac'])}"] = {
            "draws": int(r["draws"]), "mean": fnum(r["mean_of_mean_OR"]), "sd": fnum(r["sd_of_mean_OR"]),
            "mean_n_le_1": fnum(r["mean_n_le_1"]), "baseline": fnum(r["baseline_mean_OR"])}
nsub = [int(fnum(r["n_subjects"])) for r in rows("A2_permutation_null/results_full.tsv")]
N["n_subj_range"] = [min(nsub), max(nsub)]

# ---------------------------------------------------------------- A5 precision/recall/F1/Jaccard (R2)
a5 = rows("A5_prf_jaccard/results.tsv")
t3 = []
for an, ref, test, label in [("replicability", "sc", "sc", "xCell replicability"),
                             ("replicability", "sbj", "sbj", "xSubject replicability"),
                             ("preservation", "sc", "sbj", "xCell to xSubject preservation")]:
    r = pick(a5, analysis=an, reference_level=ref, test_level=test, frac="0.01")
    m = {k: st.mean(fnum(x[k]) for x in r) for k in ("odds_ratio", "precision", "recall", "F1", "jaccard")}
    m["n"] = len(r)
    m["jaccard_min"] = min(fnum(x["jaccard"]) for x in r)
    m["jaccard_max"] = max(fnum(x["jaccard"]) for x in r)
    N[f"a5_{an}_{ref}_{test}"] = m
    t3.append([label, m["n"], m["odds_ratio"], m["precision"], m["recall"], m["F1"], m["jaccard"],
               m["jaccard_min"], m["jaccard_max"]])

# ---------------------------------------------------------------- sensitivity analyses at 1% (Supplemental Table S4)
def s_from(rel, an, spec_ref, frac=0.01, universe="test_group_genes", spec_key="ref_spec"):
    rs = rows(rel)
    r = [x for x in rs if x["analysis"] == an and x[spec_key] == spec_ref and x["universe"] == universe
         and (frac is None or abs(fnum(x["frac"]) - frac) < 1e-9)]
    if not r:
        raise SystemExit(f"no rows: {rel} {an} {spec_ref} {frac} {universe}")
    return summ([fnum(x["odds_ratio"]) for x in r])


def s4(rel, an, spec, frac=0.01):
    for r in rows(rel):
        if r["analysis"] == an and (r["ref_spec"] == spec or r["test_spec"] == spec) and abs(fnum(r["frac"]) - frac) < 1e-9:
            return {"n": int(r["n"]), "mean": fnum(r["mean_OR"]), "median": fnum(r["median_OR"]),
                    "min": fnum(r["min_OR"]), "n_le_1": int(r["n_OR_le_1"])}
    raise SystemExit(f"A4 missing {an} {spec}")


A4 = "A4_xsubject_pseudobulk/summary.tsv"
for spec in ("s:std", "s:A4_sum_counts_cpm", "s:A4_sum_cellcpm", "s:A4_mean_nolog_qn", "s:A4_log_noqn", "s:A4_nolog_noqn"):
    k = spec.split(":")[1]
    N[f"a4_{k}_repl"] = s4(A4, "replicability", spec)
    N[f"a4_{k}_pres"] = s4(A4, "preservation", spec)

A7 = "A7_no_imputation/results.tsv"
N["a7_none_repl"] = s_from(A7, "replicability", "x:A7")
N["a7_none_pres"] = s_from(A7, "preservation", "x:A7")
N["a7_std_shared_repl"] = s_from(A7, "replicability", "x:std", universe="shared_genes")
N["a7_none_shared_repl"] = s_from(A7, "replicability", "x:A7", universe="shared_genes")

A8 = "A8_weighted/results.tsv"
for spec in ("A8_wrank", "A8_r", "A8_z"):
    N[f"a8_{spec}_repl"] = s_from(A8, "replicability", "x:" + spec)
    N[f"a8_{spec}_pres"] = s_from(A8, "preservation", "x:" + spec)

A6R = "A6_raw/results.tsv"
N["a6raw_rawmean_repl"] = s_from(A6R, "replicability", "x:rawmean")
N["a6raw_rawmean_pres"] = s_from(A6R, "preservation", "x:rawmean")
N["a6raw_sign_repl_sc"] = s_from(A6R, "replicability", "x:rawsign", frac=None)
N["a6raw_sign_repl_sbj"] = s_from(A6R, "replicability", "s:std", frac=None)
N["a6raw_sign_pres"] = s_from(A6R, "preservation", "x:rawsign", frac=None)

a3s = rows("A3_expression_strata/summary_stratified_baseline_edges.tsv")
for r in a3s:
    N[f"a3_strat_{r['analysis']}_{r['reference_level']}_{r['stratum']}"] = {
        "n": int(r["n"]), "n_nan": int(r["n_nan"]), "mean": fnum(r["mean_OR"]), "median": fnum(r["median_OR"]),
        "min": fnum(r["min_OR"]), "n_le_1": int(r["n_le_1"]), "mean_genes": fnum(r["mean_genes"])}
for r in rows("A3_expression_strata/summary_high_recomputed.tsv"):
    if abs(fnum(r["frac"]) - 0.01) < 1e-9:
        N[f"a3_high_{r['analysis']}_{r['reference_level']}"] = {
            "n": int(r["n"]), "n_nan": int(r["n_nan"]), "mean": fnum(r["mean_OR"]), "median": fnum(r["median_OR"]),
            "min": fnum(r["min_OR"]), "n_le_1": int(r["n_le_1"]), "mean_genes": fnum(r["mean_genes"])}
g = [int(fnum(r["genes_mean_ge1"])) for r in rows("A3_expression_strata/group_strata_counts.tsv")]
N["a3_genes_ge1"] = {"median": st.median(g), "min": min(g), "max": max(g), "n_groups": len(g)}
for r in rows("A3_expression_strata/size_correlations.tsv"):
    N[f"a3_rho_{r['level']}_{r['covariate']}_{r['outcome']}"] = fnum(r["spearman"])

for r in rows("A10_common_geneset/summary_by_level.tsv"):
    if abs(fnum(r["frac"]) - 0.01) < 1e-9:
        N[f"a10_{r['analysis']}_{r['ref_level']}_{r['rule']}"] = {
            "n": int(r["n"]), "mean": fnum(r["mean_OR"]), "median": fnum(r["median_OR"]),
            "min": fnum(r["min_OR"]), "n_le_1": int(r["n_le_1"]), "mean_genes": fnum(r["mean_genes"])}

# ---------------------------------------------------------------- A9 OR <= 1 (R2)
for r in rows("A9_or_le1/min_OR_summary.tsv"):
    N[f"a9_{r['analysis']}_{r['reference_level']}_{float(r['frac'])}"] = {
        "n": int(r["n"]), "min": fnum(r["min_OR"]), "n_le_1": int(r["n_OR_le_1"])}

# ---------------------------------------------------------------- A2 permutation null (R1-M2)
t5 = []
def a2(rel, label, filt=None):
    rs = rows(rel)
    if filt:
        rs = [r for r in rs if filt(r)]
    def m(k):
        v = [fnum(r[k]) for r in rs if r.get(k) not in (None, "")]
        return (st.mean(v), max(v)) if v else (float("nan"), float("nan"))
    d = {"n_units": len(rs), "n_perm": sorted({int(fnum(r["n_perm"])) for r in rs}),
         "draws": len({r["draw"] for r in rs}),
         "null_max_mean": m("null_max_overall")[0],
         "pos_above_max": m("frac_pos_above_null_max"), "neg_below_min": m("frac_neg_below_null_min"),
         "pos_above_q999": m("frac_pos_above_null_q0.999")}
    t5.append([label, d["n_units"], d["draws"], "/".join(map(str, d["n_perm"])), d["null_max_mean"],
               d["pos_above_max"][0] * 100, d["pos_above_max"][1] * 100, d["neg_below_min"][0] * 100,
               d["pos_above_q999"][0] * 100])
    return d
_spc = [int(fnum(r["cells_per_corr"])) for r in rows("A2_permutation_null/results_xsubject.tsv")]
N["a2_xsubject_samples_range"] = [min(_spc), max(_spc)]
N["a2_xcell_full"] = a2("A2_permutation_null/results_full.tsv", "xCell, all cells (39 groups)")
N["a2_xsubject"] = a2("A2_permutation_null/results_xsubject.tsv", "xSubject (39 groups)")
N["a2_xcell_S2"] = a2("A2_permutation_null/results_S2.tsv", "xCell, cells per subject capped at N_subj (S2)")
N["a2_xcell_S1"] = a2("A2_permutation_null/results_S1.tsv", "xCell, N_subj pooled cells (S1)")
N["a2_rosmap_full"] = a2("A2_permutation_null/results_xbulk.tsv", "xBulk ROSMAP, 636 samples",
                         lambda r: r["dataset"] == "rosmap" and r["level"] == "xbulk_full")
N["a2_rosmap_matched"] = a2("A2_permutation_null/results_xbulk.tsv", "xBulk ROSMAP, 48 samples",
                            lambda r: r["dataset"] == "rosmap" and r["level"] == "xbulk_matched")
N["a2_velmeshev_full"] = a2("A2_permutation_null/results_xbulk.tsv", "xBulk Velmeshev, 41 samples",
                            lambda r: r["dataset"] == "velmeshev" and r["level"] == "xbulk_full")

# ---------------------------------------------------------------- bulk check (R1-m2)
bc = rows("A_bulk_check/results.tsv")
N["bulk_check"] = {f"{r['stage']}|{r['reading']}": fnum(r["value"]) for r in bc}

# ---------------------------------------------------------------- write tables
def wtsv(name, header, body, caption):
    with open(OUT / name, "w", newline="") as f:
        f.write("# " + caption.replace("\n", "\n# ") + "\n")
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for r in body:
            w.writerow([f"{x:.4g}" if isinstance(x, float) else x for x in r])

wtsv("supp_table_2_thresholds.tsv",
     ["Comparison", "Edge threshold", "N comparisons", "Mean OR", "Median OR", "Min OR", "Max OR", "N OR<=1", "Mean F1", "Mean Jaccard"],
     t2, "Supplemental Table S2. Edge-threshold sensitivity (source A6_thresholds/results.tsv).")
wtsv("supp_table_3_prf_jaccard.tsv",
     ["Comparison", "N comparisons", "Mean OR", "Mean precision", "Mean recall", "Mean F1", "Mean Jaccard", "Min Jaccard", "Max Jaccard"],
     t3, "Supplemental Table S3. Agreement measures at the 1% threshold (source A5_prf_jaccard/results.tsv).")

def cell(s, nan_note=False):
    if s is None:
        return "–"
    t = f"{s['mean']:.3g} ({s['n_le_1']})"
    if nan_note and s.get("n_nan"):
        t += f" [{s['n_nan']} undefined]"
    return t

def a1row(sch):
    r = N[f"a1_{sch}_replicability_0.01"]; p = N[f"a1_{sch}_preservation_xsubject_0.01"]
    return [f"{r['mean']:.3g} ({r['mean_n_le_1']:.2g})", "–", f"{p['mean']:.3g} ({p['mean_n_le_1']:.2g})"]

t4 = [
    ["Submitted analysis (baseline)", "–", cell(N["baseline_repl_sc"]), cell(N["baseline_repl_sbj"]), cell(N["baseline_pres"])],
    ["Matched sample size, cells per subject capped at N_subj (S2; 50 draws)", "R1 M1", *a1row("S2")],
    ["Matched sample size, N_subj cells pooled, no subject step (S1; 50 draws)", "R1 M1", *a1row("S1")],
    ["xSubject pseudobulk by summing raw counts", "R1 M4", "–", cell(N["a4_A4_sum_counts_cpm_repl"]), cell(N["a4_A4_sum_counts_cpm_pres"])],
    ["xSubject pseudobulk by summing cell CPM", "R1 M4", "–", cell(N["a4_A4_sum_cellcpm_repl"]), cell(N["a4_A4_sum_cellcpm_pres"])],
    ["xSubject, no log2 transform", "R1 M4", "–", cell(N["a4_A4_mean_nolog_qn_repl"]), cell(N["a4_A4_mean_nolog_qn_pres"])],
    ["xSubject, no quantile normalization", "R1 M4", "–", cell(N["a4_A4_log_noqn_repl"]), cell(N["a4_A4_log_noqn_pres"])],
    ["xSubject, no log2 transform and no quantile normalization", "R1 M4", "–", cell(N["a4_A4_nolog_noqn_repl"]), cell(N["a4_A4_nolog_noqn_pres"])],
    ["xCell, no imputation of undetected genes", "R1 m3; R2", cell(N["a7_none_repl"]), "–", cell(N["a7_none_pres"])],
    ["xCell, genes present in both groups only: median imputation / no imputation", "R1 m3; R2",
     f"{cell(N['a7_std_shared_repl'])} / {cell(N['a7_none_shared_repl'])}", "–", "–"],
    ["xCell, subjects weighted by cell number (ranks)", "R1 m3", cell(N["a8_A8_wrank_repl"]), "–", cell(N["a8_A8_wrank_pres"])],
    ["xCell, subjects weighted by cell number (raw r)", "R1 m3", cell(N["a8_A8_r_repl"]), "–", cell(N["a8_A8_r_pres"])],
    ["xCell, Fisher z, weighted by cell number", "R1 m3", cell(N["a8_A8_z_repl"]), "–", cell(N["a8_A8_z_pres"])],
    ["xCell, raw per-subject r averaged, ranked once", "R2", cell(N["a6raw_rawmean_repl"]), "–", cell(N["a6raw_rawmean_pres"])],
    ["Edges defined by sign (mean r > 0) instead of top 1%", "R2", cell(N["a6raw_sign_repl_sc"]), cell(N["a6raw_sign_repl_sbj"]), cell(N["a6raw_sign_pres"])],
    ["Genes with mean UMI/cell >= 1 only: submitted edges restricted (all-gene reference in this universe: "
     f"{N['a3_strat_replicability_sc_all']['mean']:.3g} / {N['a3_strat_replicability_sbj_all']['mean']:.3g})", "R1 M3",
     cell(N["a3_strat_replicability_sc_high"], True), cell(N["a3_strat_replicability_sbj_high"], True), cell(N["a3_strat_preservation_sc_high"], True)],
    ["Genes with mean UMI/cell < 1 only: submitted edges restricted", "R1 M3",
     cell(N["a3_strat_replicability_sc_low"], True), cell(N["a3_strat_replicability_sbj_low"], True), cell(N["a3_strat_preservation_sc_low"], True)],
    ["Genes with mean UMI/cell >= 1 only: networks recomputed", "R1 M3",
     cell(N["a3_high_replicability_sc"], True), cell(N["a3_high_replicability_sbj"], True), cell(N["a3_high_preservation_sc"], True)],
    ["Common gene set per cell type (genes in >= 70% of datasets), networks recomputed", "R1 M3",
     cell(N["a10_replicability_sc_p70"]), cell(N["a10_replicability_sbj_p70"]), cell(N["a10_preservation_xsubject_sc_p70"])],
]
wtsv("supp_table_4_sensitivity.tsv",
     ["Analysis", "Reviewer point", "xCell replicability, mean OR (N OR<=1)", "xSubject replicability, mean OR (N OR<=1)",
      "xCell to xSubject preservation, mean OR (N OR<=1)"], t4,
     "Supplemental Table S4. Sensitivity analyses at the 1% edge threshold. For the matched-sample rows, N OR<=1 is the mean over draws.")
wtsv("supp_table_5_permutation_null.tsv",
     ["Level and sample", "N groups x draws", "Draws", "Permutations", "Mean null maximum r",
      "Positive r above null max, mean %", "Positive r above null max, max %", "Negative r below null min, mean %",
      "Positive r above null 99.9% quantile, mean %"], t5,
     "Supplemental Table S5. Comparison with a permutation null (source A2_permutation_null).")

# ---------------------------------------------------------------- Supplemental Fig S3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import random

C = {"sc": "#2a78d6", "sbj": "#eb6834", "S2": "#1baf7a", "S1": "#eda100"}
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dd"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                     "font.size": 7, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.linewidth": 0.6, "pdf.fonttype": 42})
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={"wspace": 0.42})
rng = random.Random(20261113)

ax = axs[0]
xs = list(range(len(FRACS)))
for off, lvl in ((-0.18, "sc"), (0.18, "sbj")):
    for i, fr in enumerate(FRACS):
        ors = [fnum(x["odds_ratio"]) for x in frac_rows(a6, fr, analysis="replicability", reference_level=lvl, test_level=lvl)]
        ax.scatter([i + off + rng.uniform(-0.07, 0.07) for _ in ors], ors, s=3, color=C[lvl], alpha=0.35, linewidths=0)
        ax.scatter([i + off], [st.mean(ors)], s=26, marker="o" if lvl == "sc" else "s", color=C[lvl],
                   edgecolors="white", linewidths=0.8, zorder=3)
ax.set_yscale("log", base=2)
ax.axhline(1, color=MUTED, ls="--", lw=0.7)
ax.set_xticks(xs, [f"{f*100:g}%" for f in FRACS])
ax.set_xlabel("Edge threshold (top % of pairs)")
ax.set_ylabel("Replicability odds ratio")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=C["sc"], label="xCell"),
                   Line2D([], [], marker="s", ls="", color=C["sbj"], label="xSubject")],
          frameon=False, loc="upper right", handletextpad=0.2)

ax = axs[1]
series = [("xCell, all cells", [N[f"a6_replicability_sc_sc_{f}"]["mean"] for f in FRACS], C["sc"], "o", "-"),
          ("xCell, capped (S2)", [N[f"a1_S2_replicability_{f}"]["mean"] for f in FRACS], C["S2"], "D", "-"),
          ("xCell, pooled (S1)", [N[f"a1_S1_replicability_{f}"]["mean"] for f in FRACS], C["S1"], "^", "-"),
          ("xSubject, all cells", [N[f"a6_replicability_sbj_sbj_{f}"]["mean"] for f in FRACS], C["sbj"], "s", "--")]
for lab, ys, col, mk, ls in series:
    ax.plot(xs, ys, color=col, marker=mk, ms=4, lw=1.6, ls=ls, label=lab, markeredgecolor="white", markeredgewidth=0.6)
ax.set_yscale("log", base=2)
ax.axhline(1, color=MUTED, ls=":", lw=0.7)
ax.set_xticks(xs, [f"{f*100:g}%" for f in FRACS])
ax.set_xlabel("Edge threshold (top % of pairs)")
ax.set_ylabel("Mean replicability odds ratio")
ax.legend(frameon=False, loc="upper right", fontsize=6, handlelength=2.2)

ax = axs[2]
lo, hi = 0.5, 512
for lvl, mk in (("sc", "o"), ("sbj", "s")):
    p = pairsC[lvl]
    ax.scatter([a for a, b in p], [b for a, b in p], s=9, marker=mk, color=C[lvl], alpha=0.75,
               edgecolors="white", linewidths=0.4, label=LV[lvl])
ax.plot([lo, hi], [lo, hi], color=MUTED, lw=0.7)
ax.axhline(1, color=MUTED, ls="--", lw=0.6); ax.axvline(1, color=MUTED, ls="--", lw=0.6)
ax.set_xscale("log", base=2); ax.set_yscale("log", base=2)
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_xlabel("OR, dataset A as reference")
ax.set_ylabel("OR, dataset B as reference")
ax.legend(frameon=False, loc="upper left", handletextpad=0.2)

from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
TICKS = [0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
fmt = FuncFormatter(lambda v, _: f"{v:g}")
for a in axs:
    for axis in ([a.yaxis] + ([a.xaxis] if a is axs[2] else [])):
        axis.set_major_locator(FixedLocator(TICKS)); axis.set_major_formatter(fmt); axis.set_minor_locator(NullLocator())
for axis in (axs[2].xaxis, axs[2].yaxis):
    axis.set_major_locator(FixedLocator([1, 4, 16, 64, 256]))
axs[0].set_ylim(0.3, 900); axs[1].set_ylim(0.8, 128)
for a, l in zip(axs, "ABC"):
    a.text(-0.2, 1.04, l, transform=a.transAxes, fontsize=10, fontweight="bold", color=INK)
    a.spines[["top", "right"]].set_visible(False)
    a.grid(axis="y", color=GRID, lw=0.5); a.set_axisbelow(True)
fig.savefig(OUT / "Supplemental_Fig_S3.pdf", bbox_inches="tight")
fig.savefig(OUT / "Supplemental_Fig_S3.png", dpi=600, bbox_inches="tight")

# ---------------------------------------------------------------- manifest
metas = {}
for p in sorted(FRX.glob("*/results*_meta.json")):
    try:
        m = json.loads(p.read_text())
        metas[str(p.relative_to(FRX))] = {k: m.get(k) for k in ("timestamp", "script", "script_sha256", "seed") if k in m}
    except Exception:
        pass
manifest = {
    "built": datetime.datetime.now().isoformat(timespec="seconds"),
    "builder": str(Path(__file__).name),
    "builder_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "data_version": "FRX xcell-xbulk-revision, final 2026-10-20 (frink:/home/paul/xcell-xbulk-revision), Eric Chu pipeline 041_CH4_FINAL",
    "inputs_sha256": USED,
    "input_meta": metas,
    "outputs_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.glob("*")) if p.name != "manifest.json"},
}
(OUT / "numbers.json").write_text(json.dumps(N, indent=1, default=float))
manifest["outputs_sha256"]["numbers.json"] = hashlib.sha256((OUT / "numbers.json").read_bytes()).hexdigest()
(OUT / "manifest.json").write_text(json.dumps(manifest, indent=1))
print("ok:", len(N), "numbers;", len(USED), "inputs")
