# VAE compression for plasmid dynamics simulations
from VAE_model import *
import sys
import time as timer
import pickle

def get_data_plasmid(n_target, model_index):
    X_train = np.load("./saved_sims/pgLV/I%i/pgLV_B100_T%i_fixed_train.npy"
                               % (model_index, n_target ))
    X_test = np.load("./saved_sims/pgLV/I%i/pgLV_B100_T%i_fixed_test.npy"
                               % (model_index, n_target ))
    return X_train, X_test

n_tot=5
glv_type="p"
sim_type="fixed"
n_background=100

model_index = int(sys.argv[1])

n_target = int(sys.argv[2])

trial = int(sys.argv[3])

N_EBD=[1,2,3,5,6,8,10,12,15,20]

X_train,X_test = get_data_plasmid(n_target, model_index)

TRAIN_LOSS={}
TEST_LOSS={}
t0 = timer.perf_counter()
for n_embedding in N_EBD:
    train_loss, test_loss = run(n_background,n_target,n_embedding,X_train,X_test,glv_type,sim_type,model_index,trial)
    TRAIN_LOSS[n_embedding]=train_loss
    TEST_LOSS[n_embedding]=test_loss

loss_data = [TRAIN_LOSS,TEST_LOSS]
filename = ("./vae_models/%sgLV/I%i/%sgLV_B%i_T%i_%s_trial_%i_traindoc.pkl"
            % (glv_type, model_index, glv_type, n_background, n_target, sim_type, trial))
with open(filename, 'wb') as f:
    pickle.dump(loss_data, f, protocol=pickle.HIGHEST_PROTOCOL)
train_time=timer.perf_counter()-t0
print("model: %s gLV, simulation type: %s, index: %i, target population size: %i, community size: %i"%(glv_type,sim_type,model_index,n_target,n_tot))
print("training finished, time used: %i s"%train_time)