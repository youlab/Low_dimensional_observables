# VAE compression for sparse (subsampled + interpolated) gLV simulations
# Method: Randomly subsample n_obs timepoints, linearly interpolate back to T=50
from VAE_model import *
import sys
import time as timer
import pickle
import numpy as np
import torch
from torch.utils.data import DataLoader
from pathlib import Path

model_index = int(sys.argv[1])
n_ori       = int(sys.argv[2])
n_target    = int(sys.argv[3])
sim_type    = str(sys.argv[4])
trial       = int(sys.argv[5])
n_obs       = int(sys.argv[6])

n_background = n_ori - n_target
glv_type = "b"

if sim_type == "random":
    N_EBD = [1, 2, 3, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50]
if sim_type == "fixed":
    N_EBD = [1, 2, 3, 5, 6, 8, 10, 12, 15, 20]

X_train, X_test = get_data(n_background, n_target, glv_type, sim_type, model_index)
T = X_train.shape[2]

# CHANGED: subsample + interpolate instead of truncate + pad
def subsample_and_interpolate(X, n_obs, seed=42):
    rng = np.random.default_rng(seed)
    t_all = np.arange(X.shape[2])
    t_obs = np.sort(rng.choice(X.shape[2], size=n_obs, replace=False))
    X_out = np.zeros_like(X)
    for b in range(X.shape[0]):
        for s in range(X.shape[1]):
            X_out[b, s, :] = np.interp(t_all, t_obs, X[b, s, t_obs])
    return X_out

X_train_sp = subsample_and_interpolate(X_train, n_obs)  # CHANGED
X_test_sp  = subsample_and_interpolate(X_test,  n_obs)  # CHANGED

# CHANGED: save_dir points to sparse/
save_dir = Path(f"./vae_models/{glv_type}gLV/I{model_index}/sparse")
save_dir.mkdir(parents=True, exist_ok=True)

TRAIN_LOSS, TEST_LOSS = {}, {}
t0 = timer.perf_counter()

for n_embedding in N_EBD:
    seed = 1000 * n_target + 10 * n_embedding + trial
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # ADDED
    np.random.seed(seed)              # ADDED

    lr = 1e-3
    batch_size = 64
    EPOCHS = 100

    X_train_t = torch.Tensor(X_train_sp).float()  # CHANGED
    X_test_t  = torch.Tensor(X_test_sp).float()   # CHANGED
    train_loader = DataLoader(X_train_t, batch_size=batch_size, shuffle=True)
    test_loader  = DataLoader(X_test_t,  batch_size=batch_size, shuffle=False)

    model = VAE(n_target, n_embedding)
    model.to("cuda:0")
    model.train()
    train_losses, test_losses = [], []
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    for ix_epoch in range(EPOCHS):
        train_err = train_model(model, train_loader, optimizer)
        test_err  = test_model(model, test_loader)
        train_losses.append(train_err)
        test_losses.append(test_err)

    # CHANGED: sparse/ subdirectory and sparse_ prefix in filename
    pth_path = save_dir / f"sparse_B{n_background}_T{n_target}_nobs{n_obs}_E{n_embedding}_trial{trial}.pth"
    torch.save(model.state_dict(), pth_path)

    TRAIN_LOSS[n_embedding] = train_losses
    TEST_LOSS[n_embedding]  = test_losses
    print(f"  K={n_embedding}: test={test_losses[-1]:.4e}")

# Save pkl — same format as before, just in sparse/
filename = save_dir / f"{glv_type}gLV_B{n_background}_T{n_target}_{sim_type}_nobs{n_obs}_trial{trial}_traindoc.pkl"
with open(filename, 'wb') as f:
    pickle.dump([TRAIN_LOSS, TEST_LOSS], f, protocol=pickle.HIGHEST_PROTOCOL)

train_time = timer.perf_counter() - t0
print(f"Done: model={model_index}, NT={n_target}, n_obs={n_obs}, trial={trial}, time={train_time:.0f}s")