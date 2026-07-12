from MLPVAE import *
import sys
import time as timer
import pickle

target = str(sys.argv[1])

trial = int(sys.argv[2])

t0 = timer.perf_counter()
train_loss, test_loss = run_cross_validation_growth(target,trial)

loss_data = [train_loss,test_loss]
filename = (f"./mlp_models/{target}/{target}_trial{trial}_traindoc.pkl")
with open(filename, 'wb') as f:
    pickle.dump(loss_data, f, protocol=pickle.HIGHEST_PROTOCOL)

train_time=timer.perf_counter()-t0
print("training finished, time used: %i s"%train_time)