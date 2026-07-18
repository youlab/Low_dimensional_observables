import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, leaves_list

data = np.loadtxt('./sequenced_data/sequence_composition.txt')
target_keio = {0: "EGFP", 1: "mTagBFP2", 2: "LSSmOrange", 3: "mCherry"}# target keio strains that carries different fp plasmids
tg_colors = ["#23A249","#74B3EB","#FF9000","#B53030"]

# NOTE -- column/channel permutation mismatch.
# The first four columns of sequence_composition.txt are NOT stored in the same
# order as the fluorescence channels of sequenced_target_normal.npy (channels
# 0-3 = EGFP, mTagBFP2, LSSmOrange, mCherry; see `idx_FP_map` in MLPVAE.py).
# Correlating the reconstructed absolute abundances (composition x final OD, as
# in infer_correlation.py) against the final-timepoint fluorescence identifies
# the raw column order as
#     col 0 -> mCherry, col 1 -> EGFP, col 2 -> mTagBFP2, col 3 -> LSSmOrange
# For visuliazation purpose, reorder the four target columns here
# so the plot follows the fluorescence-channel order given by `target_keio`.
# Richness, inverse Simpson and Bray-Curtis are permutation invariant, so only
# the stacked-bar colors and legend change.
# Caution: infer_correlation.py indexes the RAW columns, so its K1-K4 labels
# still refer to the unpermuted order.
fp_col_order = [1, 2, 3, 0]
data = data[:, fp_col_order + list(range(4, data.shape[1]))]

# normalize rows to relative abundances if not already
row_sums = data.sum(axis=1, keepdims=True)
nonzero = row_sums.squeeze() > 0
normed = np.zeros_like(data)
normed[nonzero] = data[nonzero] / row_sums[nonzero]

richness = (normed > 0).sum(axis=1)
order = np.argsort(richness)
normed = normed[order]

# Bray-Curtis pairwise distances between samples/rows
dist = pdist(normed, metric="braycurtis")

# Hierarchical clustering
Z = linkage(dist, method="average")  # average/UPGMA is common for Bray-Curtis

# Order of rows after clustering
order = leaves_list(Z)

data_ordered = normed[order]

n_samples, n_strains = data_ordered.shape

# colors for strains
cmap = plt.get_cmap("Set2")
colors = tg_colors+[cmap(i% cmap.N) for i in range(4,n_strains)]

fig, ax = plt.subplots(1,1,figsize=(10, 2.5))
indices = np.arange(n_samples)
# Richness / diversity axis
richness = np.sum(data_ordered > 0, axis=1)
inverse_simpson = 1 / np.sum(data_ordered**2, axis=1)
r_med  = np.median(richness)
r_q1, r_q3 = np.percentile(richness, [25, 75])
is_med = np.median(inverse_simpson)
is_q1, is_q3 = np.percentile(inverse_simpson, [25, 75])

print("=" * 45)
print(f"  Samples : {n_samples}")
print("-" * 45)
print(f"  Richness")
print(f"    Median         : {r_med:.1f}")
print(f"    IQR            : [{r_q1:.1f}, {r_q3:.1f}]")
print(f"    Min – Max      : [{int(richness.min())}, {int(richness.max())}]")
print("-" * 45)
print(f"  Inverse Simpson Diversity")
print(f"    Median         : {is_med:.2f}")
print(f"    IQR            : [{is_q1:.2f}, {is_q3:.2f}]")
print(f"    Min – Max      : [{inverse_simpson.min():.2f}, {inverse_simpson.max():.2f}]")
print("=" * 45)

ax.plot(indices, inverse_simpson, marker='o', 
         linewidth=1.5, markersize=3, color='k',zorder=20)
ax.set_xlim(-0.5, n_samples - 0.5)
ax.set_yticks([0,10,20])
ax.set_yticklabels([0,10,20],fontsize=20)
ax2 = ax.twinx()

bottom = np.zeros(n_samples)
for i in range(n_strains):
    vals = data_ordered[:, i]
    ax2.bar(indices, vals, bottom=bottom, width=1, color=colors[i], 
            edgecolor='none', label=f"k{i+1}" if i<4 else None)
    bottom += vals


ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_ylim(0, 1)

# put diversity curve above bars
ax.set_zorder(ax2.get_zorder() + 1)
ax.patch.set_visible(False)

fig.legend(loc='center right', bbox_to_anchor=(1, 0.5), fontsize=18,
           title_fontsize=18, ncol=1, frameon=True, title="strain")


#ax.set_ylabel('diversity',fontsize=16)
fig.subplots_adjust(left=0.1,right=0.86,bottom=0.12,top=0.95)  # make room for the legend
fig.savefig("./figures/community_composition.png",dpi=300)
fig.savefig("./figures/community_composition.svg",dpi=300)
plt.show()