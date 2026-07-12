# VAE compression for Fujita et al. Microbiome time series
from VAE_lab_microbiome_CV import *
import sys
import time as timer
import pickle

sample = str(sys.argv[1])

n_target=int(sys.argv[2])

test_rep=int(sys.argv[3])

N_EBD=[1,2,3,4,5,6,8,10,12,15,20,30]
X_train, X_test = get_data(sample,n_target,test_rep)
TRAIN_LOSS={}
TEST_LOSS={}
t0 = timer.perf_counter()
for n_embedding in N_EBD:
    model, train_loss, test_loss = run(sample,n_target,n_embedding,X_train,X_test,test_rep)
    TRAIN_LOSS[n_embedding]=train_loss
    TEST_LOSS[n_embedding]=test_loss

loss_data = [TRAIN_LOSS,TEST_LOSS]
filename = ("./vae_models/%s_T%i_rep%i_traindoc.pkl"
            % (sample, n_target, test_rep))
with open(filename, 'wb') as f:
    pickle.dump(loss_data, f, protocol=pickle.HIGHEST_PROTOCOL)
train_time=timer.perf_counter()-t0
print("Sample: %s, target population size: %i, trial %i"%(sample,n_target,test_rep))
print("training finished, time used: %i s"%train_time)