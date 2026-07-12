# VAE compression for Ravel et al. Microbiome time series
from VAE_all_24subjects import *
import sys
import time as timer
import pickle

FUV = int(sys.argv[1])
n_target = 1
data = np.loadtxt(f"./vaginal_embedding_FUV_CV/Ec_data_FUV={FUV}.txt")
N_TG = data[0,:]
t0 = timer.perf_counter()
Ec = int(np.ceil(data[1,n_target-1]))
X_train = get_data(n_target)
model, train_loss = run(n_target,Ec,X_train)

train_time=timer.perf_counter()-t0
print(f"training finished, {n_target}-target; time used: {int(train_time)}s")