# VAE compression for Fujita et al. Microbiome time series
from VAE_lab_microbiome_full import *
import sys
import time as timer
import pickle

sample = str(sys.argv[1])

data = np.loadtxt(f"./lab_microbiome_embedding_FUV_CV/{sample}_Ec.txt")
N_TG = data[0,:]

TRAIN_LOSS={}

for i,n_target in enumerate(N_TG):
    t0 = timer.perf_counter()
    n_target = int(n_target)
    Ec = int(np.ceil(data[1,i]))
    X_train = get_data(sample, n_target)
    model, train_loss = run(sample,n_target,Ec,X_train)
    TRAIN_LOSS[n_target]=train_loss
    train_time = timer.perf_counter() - t0
    print(f"{sample} {n_target}-target population Ec = {Ec}, train time: {train_time:.1f}")
filename = ("./vae_models_CV/%s_full_traindoc.pkl"
            % (sample))
with open(filename, 'wb') as f:
    pickle.dump(TRAIN_LOSS, f, protocol=pickle.HIGHEST_PROTOCOL)