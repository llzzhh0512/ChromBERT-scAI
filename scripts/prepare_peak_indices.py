#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Map scATAC peaks to ChromBERT genomic bin indices.

Requires bedtools. Writes/updates anndata.varm['indices'].

Example:
  # hg38
  python scripts/prepare_peak_indices.py \
    --input peaks.h5ad \
    --bins annotations/hg38_6k_1kb_region.bed \
    --output peaks_indexed.h5ad

  # mm10
  python scripts/prepare_peak_indices.py \
    --input peaks_mm10.h5ad \
    --bins annotations/mm10_5k_1kb_region.bed \
    --output peaks_mm10_indexed.h5ad
"""

import argparse
import io
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def peaks_to_centers(var_names):
    """Parse peak names like chr1:100-200 into center BED."""
    rows = []
    for name in var_names.astype(str):
        chrom, rest = name.split(":", 1)
        start, end = rest.replace("-", " ").split()[:2]
        start_i, end_i = int(start), int(end)
        center = (start_i + end_i) // 2
        rows.append((chrom, center, center + 1, name))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "peak"])


def map_with_bedtools(centers, bins_bed):
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        peak_bed = td / "peaks_center.bed"
        centers[["chrom", "start", "end", "peak"]].to_csv(
            peak_bed, sep="\t", header=False, index=False
        )
        cmd = [
            "bedtools", "intersect",
            "-a", str(peak_bed),
            "-b", str(bins_bed),
            "-wa", "-wb",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if not result.stdout.strip():
            raise RuntimeError("bedtools intersect returned no overlaps")
        cols = [
            "chrom", "start", "end", "peak",
            "bin_chrom", "bin_start", "bin_end", "bin_index",
        ]
        mapped = pd.read_csv(
            io.StringIO(result.stdout),
            sep="\t", header=None, names=cols,
        )
        mapped = mapped.drop_duplicates(subset=["peak"], keep="first")
        peak_to_idx = dict(zip(mapped["peak"], mapped["bin_index"].astype(int)))

    indices = np.full(len(centers), -1, dtype=np.int32)
    for i, peak in enumerate(centers["peak"]):
        if peak in peak_to_idx:
            indices[i] = peak_to_idx[peak]
    return indices


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="input .h5ad (cells x peaks)")
    p.add_argument("--bins", required=True, help="ChromBERT bin BED: chrom start end bin_index")
    p.add_argument("--output", required=True, help="output .h5ad with varm['indices']")
    p.add_argument("--drop-unmapped", action="store_true",
                   help="drop peaks with no bin overlap")
    args = p.parse_args()

    adata = sc.read_h5ad(args.input)
    centers = peaks_to_centers(adata.var_names)
    indices = map_with_bedtools(centers, Path(args.bins))
    n_miss = int((indices < 0).sum())
    print("mapped {}/{}; unmapped={}".format(len(indices) - n_miss, len(indices), n_miss))

    if args.drop_unmapped:
        keep = indices >= 0
        adata = adata[:, keep].copy()
        indices = indices[keep]
    elif n_miss:
        raise SystemExit(
            "Unmapped peaks present. Re-run with --drop-unmapped or fix bin BED."
        )

    adata.varm["indices"] = indices.astype(np.int32)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    print("wrote {}".format(args.output))


if __name__ == "__main__":
    main()
