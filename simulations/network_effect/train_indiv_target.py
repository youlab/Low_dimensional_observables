# VAE compression for neutral model
from VAE_model import *
import sys
import time as timer
import pickle

target_id = int(sys.argv[1])

model_index = str(sys.argv[2])

trial = int(sys.argv[3])

N_EBD = [1, 2, 3, 4, 5, 6, 8, 10, 20]

X_train, X_test = get_indiv_target(target_id, model_index)
TRAIN_LOSS={}
TEST_LOSS={}
t0 = timer.perf_counter()
for n_embedding in N_EBD:
    train_loss, test_loss = run_indiv(target_id, n_embedding, X_train, X_test, model_index, trial)
    TRAIN_LOSS[n_embedding]=train_loss
    TEST_LOSS[n_embedding]=test_loss

loss_data = [TRAIN_LOSS,TEST_LOSS]
filename = (f"./vae_models/I{model_index}/indiv_s{target_id}_E{n_embedding}_trial_{trial}_traindoc.pkl")
with open(filename, 'wb') as f:
    pickle.dump(loss_data, f, protocol=pickle.HIGHEST_PROTOCOL)
train_time=timer.perf_counter()-t0
print(f"model index: {model_index}, target id: {target_id}, trial: {trial}")
print(f"training finished, time used: {train_time} s")