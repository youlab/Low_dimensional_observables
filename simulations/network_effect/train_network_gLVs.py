# VAE compression for neutral model
from VAE_model import *
import sys
import time as timer
import pickle

n_target = int(sys.argv[1])

model_index = str(sys.argv[2])

trial = int(sys.argv[3])

N_EBD = [1, 2, 3, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50]

X_train, X_test = get_data(n_target, model_index)
TRAIN_LOSS={}
TEST_LOSS={}
t0 = timer.perf_counter()
for n_embedding in N_EBD:
    train_loss, test_loss = run(n_target,n_embedding,X_train,X_test,model_index,trial)
    TRAIN_LOSS[n_embedding]=train_loss
    TEST_LOSS[n_embedding]=test_loss

loss_data = [TRAIN_LOSS,TEST_LOSS]
filename = (f"./vae_models/I{model_index}/network_T{n_target}_trial_{trial}_traindoc.pkl")
with open(filename, 'wb') as f:
    pickle.dump(loss_data, f, protocol=pickle.HIGHEST_PROTOCOL)
train_time=timer.perf_counter()-t0
print(f"neutral model, index: {model_index}, target population size: {n_target}")
print(f"training finished, time used: {train_time} s")