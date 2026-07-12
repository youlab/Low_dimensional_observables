# VAE compression for Ravel et al. Microbiome time series
from VAE_for_vaginal_microbiome_CV import *
import sys
import time as timer
import pickle

n_target=int(sys.argv[1])

sub_id=int(sys.argv[2])

trial=int(sys.argv[3])

N_EBD=[1,2,3,4,5,6,8,10,12,15,20,30,40]
X_train, X_test = get_data(n_target,sub_id)
TRAIN_LOSS={}
TEST_LOSS={}
t0 = timer.perf_counter()

for n_embedding in N_EBD:
    model, train_loss, test_loss = run(n_target,n_embedding,X_train,X_test,sub_id,trial)
    TRAIN_LOSS[n_embedding]=train_loss
    TEST_LOSS[n_embedding]=test_loss

loss_data = [TRAIN_LOSS,TEST_LOSS]
filename = ("./vae_models_CV/T%i_sub%i_trial%i_traindoc.pkl"
            % (n_target, sub_id, trial))
with open(filename, 'wb') as f:
    pickle.dump(loss_data, f, protocol=pickle.HIGHEST_PROTOCOL)
train_time=timer.perf_counter()-t0
print("target population size: %i, test subject: %i, trial %i"%(n_target,sub_id,trial))
print("training finished, time used: %i s"%train_time)