#!/usr/bin/env python3
"""
Reproduce  figures/raw_forecast_OLS_vs_null_allmedia.png  from the raw-forecast CSVs,
using only the exported per-media CSV files (no cluster / torch / VAE / .npy needed).

Usage:
    python plot_raw_forecast_allmedia.py [CSV_DIR] [OUT_PNG]
Defaults:
    CSV_DIR = ./target_forecast_OLS_regression_raw
    OUT_PNG = ./figures/raw_forecast_OLS_vs_null_allmedia.png

Requires (local): pandas, numpy, matplotlib, seaborn, scipy.
"""
import sys
import os
import glob
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_rel

CSV_DIR = sys.argv[1] if len(sys.argv) > 1 else "./target_forecast_OLS_regression_raw"
OUT_PNG = sys.argv[2] if len(sys.argv) > 2 else "./figures/raw_forecast_OLS_vs_null_allmedia.png"

DISPLAY = {"OLS": "OLS (mu+last)", "OLS_raw": "OLS (raw window)",
           "endpoint": "endpoint", "persistence": "persistence"}
order = ["OLS (mu+last)", "OLS (raw window)", "endpoint", "persistence"]
box_palette   = {"OLS (mu+last)": "#00A400", "OLS (raw window)": "#1F77B4",
                 "endpoint": "#F08C00", "persistence": "#B0B0B0"}
strip_palette = {"OLS (mu+last)": "#006400", "OLS (raw window)": "#144870",
                 "endpoint": "#B35A00", "persistence": "#4d4d4d"}
SAMPLE_ORDER = ["Water-A", "Water-B", "Soil-A", "Soil-B", "Soil-C"]
YLIM_LO = -1.05
STEP = 0.085
box_off = lambda h: 0.2 * (h - 1.5)      # x-offset of each dodged box (4 hues, width=0.8)

# ---- load pooled-per-ASV R^2 (included ASVs only) ----
files = sorted(glob.glob(os.path.join(CSV_DIR, "*_foldR2.csv")))
if not files:
    raise SystemExit(f"No *_foldR2.csv files found in {CSV_DIR!r}")
df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
df = df[df["included"]].drop_duplicates(["sample", "delta_T", "asv", "model"]).copy()
df["R2"] = df["R2_pooled"]
df["model"] = df["model"].map(lambda m: DISPLAY.get(m, m))
samples = [s for s in SAMPLE_ORDER if s in set(df["sample"])] or sorted(df["sample"].unique())
DT = sorted(df["delta_T"].unique())


def sig_pairs(sample, dT):
    """Paired t-test across included ASVs (raw p); returns (hueA, hueB, p) with p < 0.05."""
    piv = df[(df["sample"] == sample) & (df["delta_T"] == dT)].pivot_table(
        index="asv", columns="model", values="R2")
    res = []
    for A, B in combinations(order, 2):
        if A in piv.columns and B in piv.columns:
            d = piv[[A, B]].dropna()
            if len(d) >= 2:
                _, p = ttest_rel(d[A].values, d[B].values)
                if p < 0.05:
                    res.append((order.index(A), order.index(B), float(p)))
    return sorted(res, key=lambda z: abs(z[0] - z[1]))


tops = []
for s in samples:
    for dT in DT:
        pr = sig_pairs(s, dT)
        if pr:
            gmax = df[(df["sample"] == s) & (df["delta_T"] == dT)]["R2"].max()
            tops.append(max(min(gmax, 1.0) + 0.05, 0.5) + (len(pr) - 1) * STEP + 0.06)
YLIM_HI = max([1.0] + tops)

fig, axes = plt.subplots(1, len(samples), figsize=(4.4 * len(samples), 4.8), sharey=True, squeeze=False)
axes = axes[0]
for ax, sample in zip(axes, samples):
    sub = df[df["sample"] == sample]
    sns.boxplot(data=sub, x="delta_T", y="R2", hue="model", hue_order=order, palette=box_palette,
                showfliers=False, width=0.8, ax=ax)
    sns.stripplot(data=sub, x="delta_T", y="R2", hue="model", hue_order=order, dodge=True,
                  palette=strip_palette, size=3.5, alpha=0.9, linewidth=0.3, edgecolor="white", ax=ax)
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_ylim(YLIM_LO, YLIM_HI)
    ax.set_title(sample); ax.set_xlabel(r"$\Delta T$ (days)")
    ax.set_ylabel(r"$R^2$" if sample == samples[0] else "")
    for xc, dT in enumerate(DT):
        pairs = sig_pairs(sample, dT)
        if not pairs:
            continue
        gmax = sub[sub["delta_T"] == dT]["R2"].max()
        base = max(min(gmax, 1.0) + 0.05, 0.5)
        for k, (hA, hB, p) in enumerate(pairs):
            y = base + k * STEP
            xL, xR = xc + box_off(min(hA, hB)), xc + box_off(max(hA, hB))
            ax.plot([xL, xL, xR, xR], [y - 0.02, y, y, y - 0.02], lw=0.8, color="k", clip_on=False)
            star = "***" if p < 0.001 else "**" if p < 0.01 else "*"
            ax.text((xL + xR) / 2, y + 0.005, star, ha="center", va="bottom", fontsize=8)
    n_off = int((sub[sub.model.isin(order[:2])]["R2"] < YLIM_LO).sum())
    if n_off:
        ax.text(0.5, 0.02, f"{n_off} OLS fold(s) off-scale", transform=ax.transAxes,
                ha="center", va="bottom", fontsize=7, color="#8B0000")
    h, l = ax.get_legend_handles_labels()
    (ax.legend(h[:4], l[:4], loc="lower left", fontsize=7, frameon=False)
     if sample == samples[0] else ax.legend_.remove())

fig.suptitle('Per-ASV OLS(mu+last) vs OLS(raw window) vs endpoint vs persistence — brackets: paired t-test across ASVs '
             '(* p<0.05, ** p<0.01, *** p<0.001, raw/uncorrected)')
fig.tight_layout()
os.makedirs(os.path.dirname(OUT_PNG) or ".", exist_ok=True)
# fig.savefig(OUT_PNG, dpi=300)  # figure-saving disabled for release
print(f"wrote {OUT_PNG}")
plt.show()
