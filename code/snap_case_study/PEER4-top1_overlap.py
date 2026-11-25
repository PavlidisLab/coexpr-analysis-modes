#!/usr/bin/env python3
"""
Compute overlap between gene-pair sets from two gene lists (SNAP vs top1)
per (mode x study x cell_type), from the correlation CSVs produced by PEER4 extraction.

Expected input filename pattern:
  {study}_PEER4_top1_gene_correlations_with_background{suffix}.csv

Where {suffix} is optional and typically encodes the mode, e.g.:
  _xCell, _xCellxSubject, _xSubject

Outputs:
- SNAP_vs_top1_overlap_summary_ALL_MODES.csv
- SNAP_vs_top1_overlapping_edges_ALL_MODES.csv
- optionally: per-mode CSVs (if --write_per_mode)

Author: Nairuz
Date: 2026-01-07
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np

MARKER = "_PEER4_top1_gene_correlations_with_background"


def normalize_pair(g1: str, g2: str) -> tuple:
    """Canonicalize an undirected gene pair (gene1, gene2)."""
    return (g1, g2) if g1 <= g2 else (g2, g1)


def build_pair_set(df: pd.DataFrame) -> set:
    """Build a set of normalized (gene1, gene2) pairs from a dataframe."""
    sub = df[["gene1", "gene2"]].dropna()
    g1 = sub["gene1"].astype(str).to_numpy()
    g2 = sub["gene2"].astype(str).to_numpy()
    return set(normalize_pair(a, b) for a, b in zip(g1, g2))


def overlap_metrics(a: set, b: set) -> dict:
    """Compute overlap metrics between two sets."""
    inter = a & b
    union = a | b
    return {
        "n_snap": len(a),
        "n_top1": len(b),
        "n_intersection": len(inter),
        "n_union": len(union),
        "jaccard": (len(inter) / len(union)) if union else np.nan,
        "overlap_snap_fraction": (len(inter) / len(a)) if a else np.nan,
        "overlap_top1_fraction": (len(inter) / len(b)) if b else np.nan,
    }


def parse_study_and_mode(filename: str) -> tuple[str, str]:
    """
    From:
      {study}_PEER4_top1_gene_correlations_with_background{suffix}.csv
    return (study, mode)

    mode is derived from {suffix}. If no suffix, mode="no_suffix".
    """
    if MARKER not in filename:
        raise ValueError(f"Filename does not contain expected marker {MARKER}: {filename}")

    study = filename.split(MARKER)[0]
    after = filename.split(MARKER, 1)[1]  # e.g. "_xCell.csv" or ".csv"

    if after.lower().endswith(".csv"):
        after = after[:-4]

    mode = after if after else "no_suffix"
    return study, mode


def extract_overlapping_edges(
    snap_df: pd.DataFrame,
    top1_df: pd.DataFrame,
    mode: str,
    study: str,
    cell_type: str
) -> pd.DataFrame:
    """Return overlapping SNAP ∩ top1 edges with metadata + correlations from both sources."""
    if snap_df.empty or top1_df.empty:
        return pd.DataFrame()

    snap_df = snap_df.copy()
    top1_df = top1_df.copy()

    snap_df["pair"] = snap_df.apply(lambda x: tuple(sorted((x["gene1"], x["gene2"]))), axis=1)
    top1_df["pair"] = top1_df.apply(lambda x: tuple(sorted((x["gene1"], x["gene2"]))), axis=1)

    merged = snap_df.merge(
        top1_df[["pair", "correlation"]],
        on="pair",
        how="inner",
        suffixes=("_snap", "_top1")
    )

    if merged.empty:
        return pd.DataFrame()

    merged[["gene1", "gene2"]] = pd.DataFrame(merged["pair"].tolist(), index=merged.index)

    return merged[[
        "gene1",
        "gene2",
        "correlation_snap",
        "correlation_top1"
    ]].assign(
        mode=mode,
        study=study,
        cell_type=cell_type
    )


def main(input_dir: str, output_dir: str, pattern: str, write_per_mode: bool):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(input_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files found in {input_dir} matching pattern: {pattern}")

    summary_rows = []
    overlap_edge_rows = []

    for fp in files:
        study_from_name, mode = parse_study_and_mode(fp.name)
        print(f"\nReading: {fp.name}  | inferred mode: {mode}")

        df = pd.read_csv(fp)

        required_cols = {"gene1", "gene2", "cell_type", "gene_list", "study", "correlation"}
        missing = required_cols - set(df.columns)
        if missing:
            raise ValueError(f"{fp.name} is missing required columns: {sorted(missing)}")

        # group by actual study column + cell_type
        for (study, cell_type), g in df.groupby(["study", "cell_type"], dropna=False):
            snap_df = g[g["gene_list"] == "SNAP"]
            top1_df = g[g["gene_list"] == "top1"]

            snap_pairs = build_pair_set(snap_df)
            top1_pairs = build_pair_set(top1_df)

            m = overlap_metrics(snap_pairs, top1_pairs)
            summary_rows.append({
                "mode": mode,
                "study": study,
                "cell_type": cell_type,
                **m,
                "source_file": fp.name,
                "study_from_filename": study_from_name,
            })

            overlap_edges = extract_overlapping_edges(
                snap_df=snap_df,
                top1_df=top1_df,
                mode=mode,
                study=study,
                cell_type=cell_type
            )

            if not overlap_edges.empty:
                overlap_edge_rows.append(overlap_edges)

            print(
                f"  {mode} | {study} | {cell_type}: "
                f"SNAP={m['n_snap']}  top1={m['n_top1']}  "
                f"∩={m['n_intersection']}  J={m['jaccard']:.6f}"
            )

    summary_df = pd.DataFrame(summary_rows).sort_values(["mode", "study", "cell_type"]).reset_index(drop=True)
    summary_path = output_dir / "SNAP_vs_top1_overlap_summary_ALL_MODES.csv"
    summary_df.to_csv(summary_path, index=False)

    edges_path = None
    if overlap_edge_rows:
        overlap_edges_df = pd.concat(overlap_edge_rows, ignore_index=True)
        edges_path = output_dir / "SNAP_vs_top1_overlapping_edges_ALL_MODES.csv"
        overlap_edges_df.to_csv(edges_path, index=False)

    # optional per-mode outputs
    if write_per_mode:
        for mode_val, sub in summary_df.groupby("mode"):
            safe_mode = str(mode_val).replace("/", "_")
            sub.to_csv(output_dir / f"SNAP_vs_top1_overlap_summary_{safe_mode}.csv", index=False)

        if overlap_edge_rows:
            overlap_edges_df = pd.concat(overlap_edge_rows, ignore_index=True)
            for mode_val, sub in overlap_edges_df.groupby("mode"):
                safe_mode = str(mode_val).replace("/", "_")
                sub.to_csv(output_dir / f"SNAP_vs_top1_overlapping_edges_{safe_mode}.csv", index=False)

    print(f"\nSaved:\n- {summary_path}")
    if edges_path:
        print(f"- {edges_path}")
    print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute SNAP vs top1 gene-pair overlap per mode, study, and cell type."
    )
    parser.add_argument("--input_dir", required=True, help="Directory containing per-study correlation CSVs for all modes")
    parser.add_argument("--output_dir", required=True, help="Directory to write overlap summaries")
    parser.add_argument(
        "--pattern",
        default=f"*{MARKER}*.csv",
        help="Glob pattern to match input files"
    )
    parser.add_argument(
        "--write_per_mode",
        action="store_true",
        help="If set, also write one CSV per mode"
    )

    args = parser.parse_args()
    main(args.input_dir, args.output_dir, args.pattern, args.write_per_mode)

# python PEER4-top1_overlap.py \
#  --input_dir /cosmos/data/project-data/NW-rewiring/coExpr/ling-PEER4 \
#  --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/ling-PEER4/PEER4-top1-overlaps \
#  --write_per_mode {optional} 
