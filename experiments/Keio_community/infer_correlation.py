import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import pearsonr
import math

# data shape: (n_sample, n_taxa), here n_taxa = 62
# data = your absolute abundance array

def fmt_p(p):
    # one-digit mantissa, no leading zero in exponent, e.g. 5e-3
    mant, exp = f"{p:.0e}".split("e")
    return f"{mant}e{int(exp)}"

def pairwise_pearson_filter_zeros(data, taxon_names=None, min_samples=5, threshold = 0):
    data = np.asarray(data)
    n_sample, n_taxa = data.shape

    if taxon_names is None:
        taxon_names = np.arange(n_taxa)+1

    r_mat = np.full((n_taxa, n_taxa), np.nan)
    p_mat = np.full((n_taxa, n_taxa), np.nan)
    n_mat = np.zeros((n_taxa, n_taxa), dtype=int)

    for i in range(n_taxa):
        for j in range(i, n_taxa):
            x = data[:, i]
            y = data[:, j]

            # Keep samples where both taxa are present
            mask = (x > threshold) & (y > threshold)

            x_filt = x[mask]
            y_filt = y[mask]
            n_valid = len(x_filt)

            n_mat[i, j] = n_valid
            n_mat[j, i] = n_valid

            # Need enough samples and non-constant values
            if (
                n_valid >= min_samples
                and np.std(x_filt) > 0
                and np.std(y_filt) > 0
            ):
                r, p = pearsonr(x_filt, y_filt)

                r_mat[i, j] = r
                r_mat[j, i] = r

                p_mat[i, j] = p
                p_mat[j, i] = p

    r_df = pd.DataFrame(r_mat, index=taxon_names, columns=taxon_names)
    p_df = pd.DataFrame(p_mat, index=taxon_names, columns=taxon_names)
    n_df = pd.DataFrame(n_mat, index=taxon_names, columns=taxon_names)

    return r_df, p_df, n_df


composition = np.loadtxt('./sequenced_data/sequence_composition.txt')

final_OD_normalized = np.load('./sequenced_data/sequenced_target_normal.npy')[:,4,-1]
norm_factor = np.loadtxt("./sequenced_data/max_normalization.txt")[-1]
OD_final = final_OD_normalized * norm_factor

absolute_abundance = composition * OD_final[:, np.newaxis]

fig, ax = plt.subplots(1,1,figsize=(9, 4))
indices = np.arange(absolute_abundance.shape[0])
bottom = np.zeros(absolute_abundance.shape[0])
for i in range(absolute_abundance.shape[1]):
    vals = absolute_abundance[:, i]
    ax.bar(indices, vals, bottom=bottom, width=1, edgecolor='none', label=f"k{i+1}" if i<4 else None)
    bottom += vals

absolute_abundance[absolute_abundance < 1e-3] = 1e-3
absolute_abundance = np.log10(absolute_abundance)
threshold = -3

r_df, p_df, n_df = pairwise_pearson_filter_zeros(absolute_abundance, min_samples=5, threshold=threshold)

# ---- Frequency of statistically significant pairs (p < 0.05) ----
# "enough samples" => a p-value was computed (finite); "significant" => p < 0.05
p_vals = p_df.values
n_taxa = p_vals.shape[0]

# Overall: all unique pairs (upper triangle, excluding the diagonal)
iu, ju = np.triu_indices(n_taxa, k=1)
p_upper = p_vals[iu, ju]
valid_all = np.isfinite(p_upper)
sig_all = valid_all & (p_upper < 0.05)
n_valid_all = int(valid_all.sum())
n_sig_all = int(sig_all.sum())
freq_all = n_sig_all / n_valid_all if n_valid_all else np.nan

print("=" * 56)
print("Frequency of statistically significant pairs (p < 0.05)")
print("among all pairs with enough samples")
print("-" * 56)
print(f"  Overall : {n_sig_all}/{n_valid_all} = {freq_all:.3f}")
print("-" * 56)

# Per target (taxa 0..3): partners are all other taxa j != i
per_target_freqs = []
per_target_counts = []
for i in range(4):
    row_p = p_vals[i].copy()
    row_p[i] = np.nan  # exclude self-pair
    valid_i = np.isfinite(row_p)
    sig_i = valid_i & (row_p < 0.05)
    n_valid_i = int(valid_i.sum())
    n_sig_i = int(sig_i.sum())
    freq_i = n_sig_i / n_valid_i if n_valid_i else np.nan
    per_target_freqs.append(freq_i)
    per_target_counts.append(n_sig_i)
    print(f"  K{i + 1:<2d}    : {n_sig_i}/{n_valid_i} = {freq_i:.3f}")

print("-" * 56)
print(f"  Mean freq of the 4 targets  : {np.nanmean(per_target_freqs):.3f}")
print(f"  Sig pairs per target (mean+std) : "
      f"{np.mean(per_target_counts):.1f} +/- {np.std(per_target_counts):.1f}")
print("=" * 56)

fig2, ax2 = plt.subplots(figsize=(5,5))

im = ax2.imshow(r_df.values, vmin=-1, vmax=1, cmap="coolwarm")

# Add p-values
for i in range(p_df.shape[0]):
    for j in range(p_df.shape[1]):
        p = p_df.iloc[i, j]

        if np.isfinite(p) and i != j and p<0.05:
            ax2.text(
                j, i,
                "*",
                ha="center",
                va="center",
                fontsize=6,
                color="black"
            )

n_strain = r_df.shape[0]
ax2.set_xticks([0, n_strain - 1])
ax2.set_yticks([0, n_strain - 1])
ax2.set_xticklabels([1, n_strain],fontsize=16)
ax2.set_yticklabels([1, n_strain],fontsize=16)
ax2.set_xlabel("strain id",fontsize=16)
ax2.set_ylabel("strain id",fontsize=16)
ax2.set_box_aspect(1)
cbar = fig2.colorbar(im, ax=ax2, shrink=0.6, ticks=[-1, -0.5, 0, 0.5, 1])
cbar.set_label("Pearson r",fontsize=16)

fig2.tight_layout()
# fig2.savefig("./figures/correlation_heatmap.png", dpi=300)  # figure-saving disabled for release

target_taxa = np.arange(4)   # taxa 0,1,2,3 internally

sig_pairs = {}

for i in target_taxa:
    for j in range(absolute_abundance.shape[1]):
        if i == j:
            continue

        p = p_df.iloc[i, j]
        r = r_df.iloc[i, j]
        n = n_df.iloc[i, j]

        if np.isfinite(p) and p < 0.05:
            if i not in sig_pairs:
                sig_pairs[i] = []
            sig_pairs[i].append((i, j, r, p, n))

n_cols = max([len(sig_pairs[i]) for i in target_taxa])
print(f"max number of sig pairs among the 4 targets: {n_cols}")

n_rows = 4

fig3, axes3 = plt.subplots(
    n_rows,
    n_cols,
    figsize=(1.1 * n_cols, 1.1 * n_rows),
)

for key, all_pairs in sig_pairs.items():
    for k, pair in enumerate(all_pairs):
        i, j, r, p, n = pair
        ax3 = axes3[key,k]

        x = absolute_abundance[:, j] # set the target abundances as the y axis
        y = absolute_abundance[:, i]

        mask = (x > threshold) & (y > threshold)

        x_plot = x[mask]
        y_plot = y[mask]

        x_plot = x_plot
        y_plot = y_plot

        ax3.scatter(x_plot, y_plot, s=20, alpha=0.8)

        slope, intercept = np.polyfit(x_plot, y_plot, 1)

        x_line = np.linspace(x_plot.min(), x_plot.max(), 100)
        y_line = slope * x_line + intercept

        ax3.plot(x_line, y_line, linewidth=1.5, c= 'k', zorder=10)

        ax3.set_xlim([-3,0])
        ax3.set_ylim([-3,0])
        
        ax3.set_xticks([-3, -2, -1, 0])
        ax3.set_yticks([-3, -2, -1, 0])

        # only the bottom-left subplot keeps tick labels
        if key == n_rows - 1 and k == 0:
            ax3.set_xticklabels([-3, -2, -1, 0])
            ax3.set_yticklabels([-3, -2, -1, 0])
        else:
            ax3.set_xticklabels([])
            ax3.set_yticklabels([])

        ax3.set_xlabel(rf"log$_{{10}}$K{j + 1}", fontsize=11,labelpad=0)
        ax3.set_ylabel(rf"log$_{{10}}$K{i + 1}", fontsize=11, labelpad=0)

        txt = f"r={r:.2f}\np={fmt_p(p)}\nn={n}"
        ax3.text(x=0.9, y=0.9, s=txt, fontsize=10, ha='right', va='top',transform=ax3.transAxes,
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.4, edgecolor='k'))
        
        ax3.set_aspect("equal", adjustable="box")

    # Remove empty axes
    for k in range(len(all_pairs), n_cols):
        axes3[key,k].set_visible(False)

fig3.subplots_adjust(wspace=0.18,hspace=0.45,left=0.04,right=0.99,
                     top=0.98,bottom=0.12)
# fig3.savefig("./figures/target_correlations.png", dpi=300)  # figure-saving disabled for release
plt.show()