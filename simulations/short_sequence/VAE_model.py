import numpy as np
import torch
from torch.utils.data import DataLoader
from torch.nn import functional as F
from torch import nn
from pathlib import Path

class VAE(nn.Module):

    def __init__(self, n_target, latent_dim):
        # Call parent model constructor and store hidden state variables.
        super().__init__()
        self.latent_dim = latent_dim
        self.n_target = n_target
        self.T = 50
        self.channels = 32

        self.encoder = nn.Sequential(
            nn.Conv1d(in_channels=self.n_target, out_channels=self.channels, kernel_size=3, stride=1, padding=1),
            nn.LeakyReLU(),
            nn.Conv1d(in_channels=self.channels, out_channels=self.channels, kernel_size=3, stride=1, padding=1),
            nn.LeakyReLU(),
            nn.Conv1d(in_channels=self.channels, out_channels=self.channels, kernel_size=3, stride=1, padding=1),
            nn.LeakyReLU(),
        )

        self.mean_map = nn.Linear(self.channels * self.T, self.latent_dim)

        self.std_map = nn.Linear(self.channels * self.T, self.latent_dim)

        self.fc_decoder = nn.Sequential(
            nn.Linear(self.latent_dim, self.channels * self.T),
            nn.LeakyReLU(),
        )

        self.decoder = nn.Sequential(
            nn.ConvTranspose1d(in_channels=self.channels, out_channels=self.channels, kernel_size=3, stride=1,
                               padding=1),
            nn.LeakyReLU(),
            nn.ConvTranspose1d(in_channels=self.channels, out_channels=self.channels, kernel_size=3, stride=1,
                               padding=1),
            nn.LeakyReLU(),
            nn.ConvTranspose1d(in_channels=self.channels, out_channels=self.n_target, kernel_size=3, stride=1,
                               padding=1),
            nn.LeakyReLU()
        )

    def sample(self, mean, log_var):
        """Sample a given N(0,1) normal distribution given a mean and log of variance."""

        # First compute the variance from the log variance.
        var = torch.exp(0.5 * log_var)

        # Compute a scaled distribution
        eps = torch.randn_like(var)

        # Add the vectors
        z = mean + var * eps

        return z

    def forward(self, X):
        """Forward propogate through the model, return both the reconstruction and sampled mean and standard deviation
        for the system.
        """
        pre_code = self.encoder(X)
        B, C, L = pre_code.shape
        flattened = pre_code.view(B, C * L)

        mu = self.mean_map(flattened)
        log_var = self.std_map(flattened)

        code = self.sample(mu, log_var)

        # Pass through FC layers before decoding
        post_code = self.fc_decoder(code)

        X_hat = self.decoder(post_code.view(B, C, L))

        return X_hat, code, mu, log_var


def vae_loss(x_hat, x, mu, log_var):
    "Compute the sum of BCE and KL loss for the distribution."

    # weight for the KL divergence
    alpha = 1e-4

    # Compute the reconstruction loss
    BCE = F.mse_loss(x_hat, x)

    # Compute the KL divergence of the distribution.
    KLD = -0.5 * torch.mean(1 + log_var - mu.pow(2) - log_var.exp())

    return BCE + alpha * KLD

def get_data(n_background, n_target, glv_type, sim_type, model_index):
    import os
    # For NT=7-10, use Zhengqing's directory which has correct files
    ZQ_DIR = "/hpc/group/youlab/zz294/VAE_demo/saved_sims"
    OWN_DIR = "./saved_sims_bgLV"

    train_zq = "%s/I%i/%sgLV_B%i_T%i_%s_train.npy" % (ZQ_DIR, model_index, glv_type, n_background, n_target, sim_type)
    test_zq  = "%s/I%i/%sgLV_B%i_T%i_%s_test.npy"  % (ZQ_DIR, model_index, glv_type, n_background, n_target, sim_type)
    train_own= "%s/I%i/%sgLV_B%i_T%i_%s_train.npy" % (OWN_DIR, model_index, glv_type, n_background, n_target, sim_type)
    test_own = "%s/I%i/%sgLV_B%i_T%i_%s_test.npy"  % (OWN_DIR, model_index, glv_type, n_background, n_target, sim_type)

    # Use Zhengqing's dir if our own train file is missing or corrupted
    try:
        X_train = np.load(train_own)
        if X_train.shape != (8000, n_target, 50):
            raise ValueError("Wrong shape")
    except Exception:
        X_train = np.load(train_zq)

    try:
        X_test = np.load(test_own)
        if X_test.shape != (2000, n_target, 50):
            raise ValueError("Wrong shape")
    except Exception:
        X_test = np.load(test_zq)

    return X_train, X_test

def train_model(model, data_loader, optimizer):
    num_batches = len(data_loader)
    total_loss = 0
    model.train()

    for batch in data_loader:
        batch = batch.to("cuda:0")
        optimizer.zero_grad()
        pred, code, mu, log_var = model(batch)
        loss = vae_loss(pred, batch, mu, log_var)
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
            pred, code, mu, log_var = model(batch)
            loss = vae_loss(pred, batch, mu, log_var)
            total_loss += loss.item()

    avg_loss = total_loss / num_batches
    return avg_loss

def run(n_background,n_target,n_embedding,X_train,X_test,glv_type,sim_type,model_index,trial):
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

    model = VAE(n_target, n_embedding)
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
        for param_group in optimizer.param_groups:
            param_group['lr'] = lr


    save_dir = Path(f"./vae_models/{glv_type}gLV/I{model_index}")
    save_dir.mkdir(parents=True, exist_ok=True)

    save_path = save_dir / f"{sim_type}_B{n_background}_T{n_target}_E{n_embedding}_trial{trial}.pth"

    torch.save(model.state_dict(), save_path)

    print("training finished, with starting MSE %1.1e, and ending error %1.1e" % (test_losses[0], test_losses[-1]))
    return train_losses, test_losses