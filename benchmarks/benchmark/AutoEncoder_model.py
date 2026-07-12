import numpy as np
import torch
from torch.utils.data import DataLoader
from torch.nn import functional as F
from torch import nn
from pathlib import Path

import pickle

device = "cuda:0"
class AutoEncoder(nn.Module):

    def __init__(self, n_target, latent_dim):
        # Deterministic autoencoder with MLP encoder/decoder blocks.
        super().__init__()
        self.latent_dim = latent_dim
        self.n_target = n_target
        self.T = 50
        self.hidden = 32

        input_dim = self.n_target * self.T

        # 3-layer encoder: each layer outputs 32 with LeakyReLU activations.
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, self.hidden),
            nn.LeakyReLU(),
            nn.Linear(self.hidden, self.hidden),
            nn.LeakyReLU(),
            nn.Linear(self.hidden, self.hidden),
            nn.LeakyReLU(),
            nn.Linear(self.hidden, self.latent_dim)
        )

        # 3-layer decoder: first two layers output 32; last maps back to the signal size.
        self.decoder = nn.Sequential(
            nn.Linear(self.latent_dim, self.hidden),
            nn.LeakyReLU(),
            nn.Linear(self.hidden, self.hidden),
            nn.LeakyReLU(),
            nn.Linear(self.hidden, self.hidden),
            nn.LeakyReLU(),
            nn.Linear(self.hidden, input_dim),
        )

    def forward(self, X):
        """Forward propagate through a deterministic autoencoder."""
        B = X.shape[0]
        x_flat = X.view(B, -1)

        code = self.encoder(x_flat)
        x_hat_flat = self.decoder(code)
        X_hat = x_hat_flat.view(B, self.n_target, self.T)

        return X_hat, code

def get_data(n_background, n_target, glv_type, sim_type, model_index):
    X_train = np.load("../../simulations/VAE_gLV_simulation/saved_sims/%sgLV/I%i/%sgLV_B%i_T%i_%s_train.npy"
                               % (glv_type, model_index, glv_type, n_background, n_target, sim_type, ))
    X_test = np.load("../../simulations/VAE_gLV_simulation/saved_sims/%sgLV/I%i/%sgLV_B%i_T%i_%s_test.npy"
                               % (glv_type, model_index, glv_type, n_background, n_target, sim_type, ))

    return X_train, X_test

def train_model(model, data_loader, optimizer):
    num_batches = len(data_loader)
    total_loss = 0
    model.train()

    for batch in data_loader:
        batch = batch.to("cuda:0")
        optimizer.zero_grad()
        pred, code = model(batch)
        loss = F.mse_loss(pred, batch)
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    avg_loss = total_loss / num_batches
    return avg_loss

def test_model(model, data_loader):

    num_batches = len(data_loader)
    total_loss = 0

    model.eval()
    with torch.no_grad():
        for batch in data_loader:
            batch = batch.to("cuda:0")
            pred, code = model(batch)
            loss = F.mse_loss(pred, batch)
            total_loss += loss.item()

    avg_loss = total_loss / num_batches
    return avg_loss

def run(n_background,n_target,n_embedding,X_train,X_test,glv_type,sim_type,model_index,trial=1):
    # set the random seed so the three trials have different initial model weights
    seed = 1000 * n_target + 10*n_embedding + trial
    torch.manual_seed(seed)
    
    lr = 1e-3
    batch_size = 64
    EPOCHS = 100
    
    X_train = torch.Tensor(X_train).float()
    X_test = torch.Tensor(X_test).float()
    train_loader = DataLoader(X_train, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(X_test, batch_size=batch_size, shuffle=False)

    model = AutoEncoder(n_target, n_embedding)
    model.to("cuda:0");
    model.train()
    train_losses = []
    test_losses = []
    
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    for ix_epoch in range(EPOCHS):
        train_err = train_model(model, train_loader, optimizer)
        test_err = test_model(model, test_loader)
        train_losses.append(train_err)
        test_losses.append(test_err)

    save_path = Path(
        f"./AE_models/{glv_type}gLV/I{model_index}/{sim_type}_B{n_background}_T{n_target}_E{n_embedding}_trial{trial}.pth"
    )
    save_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), save_path)

    print("training finished, with starting MSE %1.1e, and ending error %1.1e" % (test_losses[0], test_losses[-1]))
    return train_losses, test_losses