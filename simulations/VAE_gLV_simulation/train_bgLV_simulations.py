# VAE compression for bounded gLV simulations
from VAE_model import *
import sys
import time as timer
import pickle

model_index = int(sys.argv[1])

n_ori = int(sys.argv[2])

n_target = int(sys.argv[3])

sim_type = str(sys.argv[4])

trial = int(sys.argv[5])

n_background = n_ori-n_target

glv_type = "b"

if sim_type == "random":
    N_EBD = [1, 2, 3, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50]
if sim_type=="fixed":
    N_EBD = [1,2,3,5,6,8,10,12,15,20]
X_train, X_test = get_data(n_background, n_target, glv_type = glv_type, sim_type = sim_type, model_index = model_index)
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
print("model: %s gLV, simulation type: %s, index: %i, target population size: %i, community size: %i"%(glv_type,sim_type,model_index,n_target,n_ori))
print("training finished, time used: %i s"%train_time)