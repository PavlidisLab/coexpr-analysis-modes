#!/usr/bin/env python3
"""Supplemental Fig S4: permutation-null comparison per data group (R1 M2).
Input : ../../results/A2_permutation_null/results_{full,S2,S1,xsubject,xbulk}.tsv
Output: out/Supplemental_Fig_S4.{pdf,png}; input sha256 appended to out/manifest_figS4.json
"""
import csv, json, hashlib, statistics as st, random
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
A2 = REPO / "results" / "A2_permutation_null"
OUT = REPO / "results" / "supp_items"
SHA = {}


def rows(name):
    p = A2 / name
    SHA[name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return list(csv.DictReader(open(p), delimiter="\t"))


def per_group(rs, key):
    """Average over draws within a data group (or bulk dataset)."""
    g = defaultdict(list)
    for r in rs:
        g[(r["dataset"], r["cell_type"])].append(r)
    return [{k: st.mean(float(x[k]) for x in v) for k in key} for v in g.values()]


K = ["frac_pos_above_null_max", "frac_pos_above_null_q0.999", "null_max_overall", "cells_per_corr"]
bulk = rows("results_xbulk.tsv")
LEVELS = [  # label, rows, colour, marker
    ("xCell, all cells", per_group(rows("results_full.tsv"), K), "#2a78d6", "o"),
    ("xCell, capped (S2)", per_group(rows("results_S2.tsv"), K), "#1baf7a", "D"),
    ("xCell, pooled (S1)", per_group(rows("results_S1.tsv"), K), "#eda100", "^"),
    ("xSubject", per_group(rows("results_xsubject.tsv"), K), "#eb6834", "s"),
    ("xBulk, ROSMAP (636)", per_group([r for r in bulk if r["dataset"] == "rosmap" and r["level"] == "xbulk_full"], K), "#e87ba4", "v"),
    ("xBulk, ROSMAP (48)", [{k: float(r[k]) for k in K} for r in bulk if r["dataset"] == "rosmap" and r["level"] == "xbulk_matched"], "#e87ba4", "v"),
    ("xBulk, Velmeshev (41)", per_group([r for r in bulk if r["dataset"] == "velmeshev" and r["level"] == "xbulk_full"], K), "#e87ba4", "v"),
]

INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dd"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                     "font.size": 7, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.linewidth": 0.6, "pdf.fonttype": 42})
fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.1), gridspec_kw={"wspace": 0.38, "width_ratios": [1.25, 1.25, 1]})
rng = random.Random(20261113)
FLOOR = 1e-4  # percent; zeros are drawn here and labelled "0"
fmt = FuncFormatter(lambda v, _: "0" if v <= FLOOR * 1.01 else f"{v:g}")

for ax, key, ylab in ((axs[0], "frac_pos_above_null_max", "% of positive r above null maximum"),
                      (axs[1], "frac_pos_above_null_q0.999", "% of positive r above null 99.9% quantile")):
    for i, (lab, rs, col, mk) in enumerate(LEVELS):
        ys = [max(r[key] * 100, FLOOR) for r in rs]
        xs = [i + rng.uniform(-0.15, 0.15) if len(ys) > 1 else i for _ in ys]
        ax.scatter(xs, ys, s=10, color=col, marker=mk, alpha=0.75, edgecolors="white", linewidths=0.3, zorder=3)
        if len(ys) > 1:
            m = max(st.mean(r[key] for r in rs) * 100, FLOOR)
            ax.plot([i - 0.28, i + 0.28], [m, m], color=INK, lw=1.2, zorder=4)
    ax.set_yscale("log")
    ax.set_ylim(FLOOR * 0.6, 100)
    ax.yaxis.set_major_locator(FixedLocator([FLOOR, 0.001, 0.01, 0.1, 1, 10, 100]))
    ax.yaxis.set_major_formatter(fmt); ax.yaxis.set_minor_locator(NullLocator())
    ax.axhline(FLOOR * 2.2, color=GRID, lw=0.6, ls=":")
    ax.set_xticks(range(len(LEVELS)), [l for l, *_ in LEVELS], fontsize=6, rotation=40, ha="right", rotation_mode="anchor")
    ax.set_ylabel(ylab)

ax = axs[2]
seen = set()
for lab, rs, col, mk in LEVELS:
    lab = "xBulk" if lab.startswith("xBulk") else lab
    ax.scatter([r["cells_per_corr"] for r in rs], [r["null_max_overall"] for r in rs], s=10, color=col, marker=mk,
               alpha=0.75, edgecolors="white", linewidths=0.3, label=None if lab in seen else lab)
    seen.add(lab)
ax.set_xscale("log")
ax.set_xlabel("Samples per correlation")
ax.set_ylabel("Null maximum r")
ax.set_ylim(0, 1.05)
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
ax.legend(frameon=False, fontsize=5.6, loc="center right", bbox_to_anchor=(1.0, 0.6), handletextpad=0.1)

for a, l in zip(axs, "ABC"):
    a.text(-0.22, 1.03, l, transform=a.transAxes, fontsize=10, fontweight="bold", color=INK)
    a.spines[["top", "right"]].set_visible(False)
    a.grid(axis="y", color=GRID, lw=0.5); a.set_axisbelow(True)
# guards for the caption's qualitative statements
lv = {l.split(",")[0] if not l.startswith("xBulk") else l: rs for l, rs, *_ in LEVELS}
mean = lambda rs, k: st.mean(r[k] for r in rs)
assert mean(LEVELS[0][1], "frac_pos_above_null_max") > mean(LEVELS[3][1], "frac_pos_above_null_max")
assert mean(LEVELS[3][1], "frac_pos_above_null_q0.999") > mean(LEVELS[0][1], "frac_pos_above_null_q0.999")
assert min(r["null_max_overall"] for r in LEVELS[2][1]) > 0.99 and mean(LEVELS[3][1], "null_max_overall") > 0.8
assert max(r["null_max_overall"] for r in LEVELS[0][1]) < 0.5
fig.savefig(OUT / "Supplemental_Fig_S4.pdf", bbox_inches="tight")
fig.savefig(OUT / "Supplemental_Fig_S4.png", dpi=600, bbox_inches="tight")
(OUT / "manifest_figS4.json").write_text(json.dumps({"inputs_sha256": SHA}, indent=1))
print("ok")
