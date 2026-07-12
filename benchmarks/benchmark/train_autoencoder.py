# AutoEncoder compression for bounded gLV simulations
import sys
import time as timer
import pickle
from pathlib import Path
from AutoEncoder_model import *

model_index = int(sys.argv[1])

n_ori = int(sys.argv[2])

n_target = int(sys.argv[3])

sim_type = str(sys.argv[4])

trial = int(sys.argv[5])

n_background = n_ori-n_target

glv_type = "d"

N_EBD = [1,2,3,5,6,8,10,12,15,20,30]

X_train, X_test = get_data(n_background, n_target, glv_type = glv_type, sim_type = sim_type, model_index = model_index)
TRAIN_LOSS={}
TEST_LOSS={}
t0 = timer.perf_counter()
for n_embedding in N_EBD:
    train_loss, test_loss = run(n_background,n_target,n_embedding,X_train,X_test,glv_type,sim_type,model_index,trial)
    TRAIN_LOSS[n_embedding]=train_loss
    TEST_LOSS[n_embedding]=test_loss

loss_data = [TRAIN_LOSS,TEST_LOSS]

filename = Path(
    f"./AE_models/{glv_type}gLV/I{model_index}/{glv_type}gLV_B{n_background}_T{n_target}_{sim_type}_trial_{trial}_traindoc.pkl"
)
filename.parent.mkdir(parents=True, exist_ok=True)
with open(filename, 'wb') as f:
    pickle.dump(loss_data, f, protocol=pickle.HIGHEST_PROTOCOL)
train_time=timer.perf_counter()-t0
print("model: %s gLV, simulation type: %s, index: %i, target population size: %i, community size: %i"%(glv_type,sim_type,model_index,n_target,n_ori))
print("training finished, time used: %i s"%train_time)