# VAE compression for mcrm
from VAE_model import *
import sys
import time as timer
import pickle

var_type = str(sys.argv[1])

n_target = int(sys.argv[2])

model_index = str(sys.argv[3])

trial = int(sys.argv[4])

if var_type == "species":
    N_EBD = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30]
elif var_type == "resource":
    N_EBD = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20]

X_train, X_test = get_data(n_target, var_type, model_index)
TRAIN_LOSS={}
TEST_LOSS={}
t0 = timer.perf_counter()
for n_embedding in N_EBD:
    train_loss, test_loss = run(var_type,n_target,n_embedding,X_train,X_test,model_index,trial)
    TRAIN_LOSS[n_embedding]=train_loss
    TEST_LOSS[n_embedding]=test_loss

loss_data = [TRAIN_LOSS,TEST_LOSS]
filename = (f"./vae_models/I{model_index}/{var_type}_T{n_target}_trial_{trial}_traindoc.pkl")
with open(filename, 'wb') as f:
    pickle.dump(loss_data, f, protocol=pickle.HIGHEST_PROTOCOL)
train_time=timer.perf_counter()-t0
print(f"mcrm model, target variable: {var_type}, index: {model_index}, target population size: {n_target}")
print(f"training finished, time used: {train_time} s")