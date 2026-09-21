# Levina-Bickel intrinsic-dimension estimates for the bounded gLV simulations (random background).
#
# Data: all 10,000 simulations of each dataset, rebuilt by concatenating the train (8,000) and
# test (2,000) splits in
#   simulations/VAE_gLV_simulation/saved_sims/bgLV/I*/bgLV_B*_T*_random_{train,test}.npy
# (restore with `bash zenodo/zenodo_download.sh 10.5281/zenodo.21368291 glv_saved_sims`).
# NOTE: this is the same sample size (10,000) as the original analysis, which read the raw MATLAB
# output. Unlike the raw output, the saved splits are normalized: each dataset is divided by a
# single global maximum (shared by train and test) and the sample order is shuffled. Neither
# changes the Levina-Bickel estimate, which depends only on ratios of neighbor distances.
#
# Ec estimates are cached in Ec_estimates/ and reused if present; delete a file to recompute it.
# The Ec_estimates/*.txt shipped with the repo were computed on the raw (unnormalized) 10,000 samples.
import os
from intrinsic_dimension import *
from scipy.optimize import curve_fit
from scipy import stats

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, "../../simulations/VAE_gLV_simulation/saved_sims/bgLV")
ec_dir = os.path.join(script_dir, "Ec_estimates")

colors = ['#03468F', '#FF9000','#00A400','#FF9797','#808080',]

glv_type="b"
sim_type="random"
n_ori=100
NT = np.arange(1,11,1)

k_list1 = np.arange(3,21,1)
k_list2 = np.arange(10,21,1)
B        = 100          # number of subsamples (replicate fits)
n_boot   = 2000          # rows per subsample, drawn without replacement

fig,axes=plt.subplots(len(NT),3,figsize=(6,1.2*len(NT)))
fig2,axes2=plt.subplots(1,3,figsize=(6,2.4))
for i,model_index in enumerate(np.arange(1,4,1)):
    ax=axes[0,i]
    ax.set_title("community %i"%(model_index))
    ax2=axes2[i]
    ax2.set_title("community %i"%(model_index))
    ec_file = os.path.join(ec_dir, f"{glv_type}gLV_I{model_index}_{sim_type}.txt")
    use_cache = os.path.exists(ec_file)
    if use_cache:
        Ec_LB = np.loadtxt(ec_file)
        assert np.array_equal(Ec_LB[0,:], NT), f"{ec_file} was computed for a different NT"
        print(f"model {model_index}: using cached {ec_file}")
    else:
        Ec_LB = np.zeros((3,len(NT)))
        Ec_LB[0,:] = NT
    for j,n_target in enumerate(NT):
        n_background=n_ori-n_target
        ax=axes[j,i]

        data_file = os.path.join(data_dir, "I%i/%sgLV_B%i_T%i_%s_%%s.npy"
                                 % (model_index, glv_type, n_background, n_target, sim_type, ))
        X=np.concatenate([np.load(data_file % "train"), np.load(data_file % "test")])  # 8,000 + 2,000 = 10,000
        X=X.reshape(X.shape[0], -1)  # (10000, n_target, 50) -> (10000, n_target*50), one row per sample

        dists, _, dimensions = fit(x=X, k_list=k_list1, n_jobs=8,)
        ax.scatter(k_list1, dimensions)#edgecolor=colors[j],

        if i==0:
            if n_target<8:
                ax.text(x=0.5, y=0.7, s="n = %i" % n_target, fontsize=12.5, transform=ax.transAxes, c="k")#c=colors[j],
            else:
                ax.text(x=0.5, y=0.1, s="n = %i" % n_target, fontsize=12.5, transform=ax.transAxes, c="k")#c=colors[j],

        ax.set_xlim([2,21])
        ax.set_xticks([3,10,20])
        ax.set_ylim([0,40])
        ax.set_yticks([0,20,40])

        if j!=len(NT)-1 or i!=0:
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        if not use_cache:
            # --- estimate the intrinsic dimension and the storage --------------------------------------------------------------
            boot_means = np.empty(B)

            # --- subsampling loop -----------------------------------------------------
            for b in range(B):
                # 1) draw n_boot distinct rows from X (*without* replacement, i.e. subsampling, not a bootstrap)
                idx = np.random.choice(np.arange(X.shape[0]), size=n_boot, replace=False)  # indices
                X_boot = X[idx]

                # 2) run the intrinsic-dimension estimator on this subsample
                _, _, dims = fit(x=X_boot, k_list=k_list2, n_jobs=8)

                # 3) store the single summary number you care about
                boot_means[b] = dims.mean()

            # --- mean and SD across subsamples ----------------------------------------
            mean_hat = boot_means.mean()
            sd_hat = boot_means.std(ddof=1)  # SD across subsamples (reported as the error bar)

            Ec_LB[1,j] = mean_hat
            Ec_LB[2,j] = sd_hat

        ax2.errorbar(NT[j], Ec_LB[1,j], yerr=Ec_LB[2,j],
                     color=colors[i], linewidth=0,
                     markersize=0, capsize=3, elinewidth=1.2, capthick=1.2, zorder=20)
        ax2.scatter(NT[j], Ec_LB[1,j], marker='o', facecolor="w", edgecolor = colors[i], s = 100, zorder = -10, linewidth=1.2)
    if not use_cache:
        os.makedirs(ec_dir, exist_ok=True)
        np.savetxt(ec_file,Ec_LB)
    # ax2.set_xlim([0,11])
    # ax2.set_ylim([0,11])
    # ax2.set_xticks([0, 5, 10])
    # ax2.set_yticks([0, 5, 10])
    if i != 0:
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
    # sigma-weighted least squares (sigma = SD across subsamples, treated as absolute), as in SM 1.1.5
    (slope, intercept), cov = curve_fit(lambda n, a, b: a*n + b, NT, Ec_LB[1,:], sigma=Ec_LB[2,:], absolute_sigma=True)
    slope_se = np.sqrt(cov[0,0])
    p_value = 2*stats.t.sf(abs(slope/slope_se), df=len(NT)-2)  # two-sided t-test, slope != 0
    w = 1/Ec_LB[2,:]**2
    resid = Ec_LB[1,:] - (slope*NT + intercept)
    r2_w = 1 - np.sum(w*resid**2)/np.sum(w*(Ec_LB[1,:] - np.average(Ec_LB[1,:], weights=w))**2)
    print(f"community {model_index}: slope {slope:.3f} ± {slope_se:.3f} SE, p = {p_value:.0e}; weighted R2 = {r2_w:.3f}")
    p2 = np.poly1d([slope, intercept])
    x=np.linspace(0,12,20)
    y=p2(x)
    ax2.plot(x,y,linewidth=1,c="k",zorder=-20)

fig.supxlabel('Number of Neighbors (k)')
fig.supylabel('Estimated Intrinsic Dimension')
fig.subplots_adjust(left=0.1,bottom=0.05,right=0.95,top=0.95,hspace=0.1,wspace=0.1)
axes2[0].set_xlabel("observables n",fontsize=13)
axes2[0].set_ylabel("intrinsic dimension",fontsize=13)
fig2.tight_layout()
plt.show()