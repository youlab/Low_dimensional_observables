# VAE compression for bounded gLV simulations
from VAE_model_retrain import *
import sys
import time as timer

model_index = int(sys.argv[1])

n_ori = int(sys.argv[2])

n_target = int(sys.argv[3])

sim_type = str(sys.argv[4])

n_background = n_ori-n_target

glv_type = "b"

X_train, X_test = get_data(n_background, n_target, glv_type = glv_type, sim_type = sim_type, model_index = model_index)

Ec_data = np.loadtxt(f"./saved_data/bgLV/Ec_bgLV_model{model_index}_{sim_type}.txt")
N_TG = Ec_data[0,:]
Ec = Ec_data[1,:]
k = np.where(N_TG==n_target)[0][0]
n_embedding = int(np.ceil(Ec[k]))

t0 = timer.perf_counter()
train_loss, test_loss = retrain(n_background,n_target,n_embedding,X_train,X_test,glv_type,sim_type,model_index)

train_time=timer.perf_counter()-t0
print("model: %s gLV, simulation type: %s, index: %i, target population size: %i, community size: %i"%(glv_type,sim_type,model_index,n_target,n_ori))
print("training finished, time used: %i s"%train_time)