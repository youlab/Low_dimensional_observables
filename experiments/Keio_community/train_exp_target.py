# VAE compression for experimental random community assembly
from VAE_model import *
import sys
import time as timer
import pickle

n_target = int(sys.argv[1])

trial = int(sys.argv[2])

N_EBD = [1,2,3,4,5,6,8,10,12,15,20]

X_train, X_test = get_data(n_target)
TRAIN_LOSS={}
TEST_LOSS={}
t0 = timer.perf_counter()
for n_embedding in N_EBD:
    train_loss, test_loss = run(n_target,n_embedding,X_train,X_test,trial)
    TRAIN_LOSS[n_embedding]=train_loss
    TEST_LOSS[n_embedding]=test_loss

loss_data = [TRAIN_LOSS,TEST_LOSS]
filename = ("./vae_models/T%i_trial_%i_traindoc.pkl"
            % (n_target, trial))
with open(filename, 'wb') as f:
    pickle.dump(loss_data, f, protocol=pickle.HIGHEST_PROTOCOL)
train_time=timer.perf_counter()-t0
print("training finished, time used: %i s"%train_time)