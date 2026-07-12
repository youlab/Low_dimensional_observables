import numpy as np
import torch
from torch.utils.data import DataLoader
from torch import nn

import time as timer
import pickle

import sys

# Configure GPU if available
if torch.cuda.is_available():

    device = "cuda:0"
else:
    device = "cpu"

class VAE(nn.Module):

    def __init__(self, n_target, latent_dim):
        # Call parent model constructor and store hidden state variables.
        super().__init__()
        self.latent_dim = latent_dim
        self.n_target = n_target
        self.channels = 16

        self.encoder = nn.Sequential(
            nn.Conv1d(in_channels=self.n_target, out_channels=self.channels, kernel_size=3, stride=1, padding=1),
            nn.LeakyReLU(),
            nn.Conv1d(in_channels=self.channels, out_channels=self.channels, kernel_size=3, stride=1, padding=1),
            nn.LeakyReLU(),
            nn.Conv1d(in_channels=self.channels, out_channels=self.channels, kernel_size=3, stride=1, padding=1),
            nn.LeakyReLU(),
        )

        self.mean_map = nn.Linear(self.channels * 30, self.latent_dim)

        self.std_map = nn.Linear(self.channels * 30, self.latent_dim)

        self.fc_decoder = nn.Sequential(
            nn.Linear(self.latent_dim, self.channels * 30),
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


class VAE_MLP(torch.nn.Module):

    def __init__(self, VAE, n_embedding, n_asv):
        super().__init__()

        self.channels = 32

        self.encoder = VAE.encoder.eval()
        self.mean_map = VAE.mean_map.eval()
        self.std_map = VAE.std_map.eval()
        self.mlp = nn.Sequential(
            nn.Linear(n_asv * n_embedding, self.channels),
            nn.LeakyReLU(),
            nn.Linear(self.channels, self.channels),
            nn.LeakyReLU(),
            nn.Linear(self.channels, 1),
        )

        # Freeze the parameters of the VAE
        for param in self.encoder.parameters():
            param.requires_grad = False
        for param in self.mean_map.parameters():
            param.requires_grad = False
        for param in self.std_map.parameters():
            param.requires_grad = False

    def sample(self, mean, log_var):
        """Sample a given N(0,1) normal distribution given a mean and log of variance."""

        # First compute the variance from the log variance.
        var = torch.exp(0.5 * log_var)

        # Compute a scaled distribution
        eps = torch.randn_like(var)

        # Add the vectors
        z = mean + var * eps

        return z

    def forward(self, Xs, n_asv):
        """Forward propogate through the model, return both the reconstruction and sampled mean and standard deviation
        for the system.
        """

        embedding = []
        for i in range(n_asv):
            X = Xs[:, i:i + 1, :]

            pre_code = self.encoder(X)
            B, C, L = pre_code.shape
            flattened = pre_code.view(B, C * L)

            mu = self.mean_map(flattened)
            log_var = self.std_map(flattened)

            code = self.sample(mu, log_var)
            embedding.append(code)

        embedding = torch.cat(embedding, dim=1)

        b_div = self.mlp(embedding).squeeze()

        return b_div


class CommShiftDataset(torch.utils.data.Dataset):

    def __init__(self, segment, b_div):
        self.X = segment
        self.Y = b_div

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return self.X[idx], self.Y[idx]


def get_collapse_data(sample, n_asv, threshold=0.3):
    ts_data = np.load(f"./lab_microbiome_datasets_CV/ranked/{sample}_ranked_relative.npy")
    abruptness = np.loadtxt(
        f"./community_shift/{sample}_abruptness.txt")
    Segments = []
    BDIV = []  # beta-diversity
    T = 30
    for start in range(ts_data.shape[1] - T - 5):
        seg = ts_data[:, start:start + T, 0:n_asv]
        bdiv = abruptness[:, start + T - 1]
        Segments.append(seg)
        BDIV.append(bdiv)
    Segments = np.array(Segments)  # shape (x,8,30,n_asv)
    Segments = np.transpose(Segments, axes=(1, 0, 3, 2))  # shape (8,x,n_asv,30)

    BDIV = np.array(BDIV)  # shape (x,8)
    BDIV = np.transpose(BDIV, axes=(1, 0))

    collapse = (BDIV > threshold).astype(np.int64)

    return Segments, collapse


def train_model(model, data_loader, optimizer, criterion, n_asv):
    model.train()
    total_loss = 0.0

    for x, y in data_loader:
        x, y = x.to(device), y.to(device).float()  # BCE expects float targets

        optimizer.zero_grad()
        logits = model(x, n_asv).squeeze(dim=-1)  # (batch,)  or (batch,1)
        loss = criterion(logits, y)  # BCEWithLogitsLoss
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    return total_loss / len(data_loader)


def test_model(model, data_loader, criterion, n_asv):
    model.eval()
    total_loss = 0.0

    with torch.no_grad():
        for x, y in data_loader:
            x, y = x.to(device), y.to(device).float()
            logits = model(x, n_asv).squeeze(dim=-1)
            loss = criterion(logits, y)
            total_loss += loss.item()

    return total_loss / len(data_loader)


def run_collapse_cv(sample, n_asv, threshold, trial):
    torch.manual_seed(10*n_asv+trial)
    lr0 = 1e-3
    EPOCHS = 40
    batch_sz = 16
    lr_decay = 0.98

    # ── Load pretrained VAE ──────────────────────────────────────────────
    Ec = int(np.ceil(np.loadtxt(f"./lab_microbiome_embedding_FUV_CV/{sample}_Ec.txt")[1, 0]))
    vae = VAE(1, Ec)
    vae.load_state_dict(torch.load(f"./vae_models_CV/{sample}_T1_E{Ec}_full.pth",
                                   weights_only=True))
    vae.eval().to(device)

    # ── Load full dataset ────────────────────────────────────────────────
    Segments, P = get_collapse_data(sample, n_asv, threshold)

    all_train_losses, all_test_losses = [], []

    for i in range(8):

        # ----- split -----
        X_test = torch.tensor(Segments[i]).float()
        y_test = torch.tensor(P[i]).float()
        X_train = torch.tensor(np.vstack([Segments[j] for j in range(8) if j != i])).float()
        y_train = torch.tensor(np.hstack([P[j] for j in range(8) if j != i])).float()

        train_ds = CommShiftDataset(X_train, y_train)
        test_ds = CommShiftDataset(X_test, y_test)
        train_ld = DataLoader(train_ds, batch_size=batch_sz, shuffle=True)
        test_ld = DataLoader(test_ds, batch_size=batch_sz, shuffle=False)

        # ----- pos_weight -----
        n_pos = (y_train == 1).sum()  # assumes labels are 0/1
        n_neg = y_train.numel() - n_pos
        pos_w = (n_neg / n_pos).clamp(max=30)  # optional clip for stability
        criterion = torch.nn.BCEWithLogitsLoss(
            pos_weight=pos_w.to(device)  # scalar tensor
        )

        # ----- model & optimiser -----
        model = VAE_MLP(vae, Ec, n_asv).float().to(device)
        optimiser = torch.optim.Adam(model.mlp.parameters(), lr=lr0)

        # ----- training -----
        train_losses, test_losses = [], []
        for epoch in range(EPOCHS):
            tr_loss = train_model(model, train_ld, optimiser, criterion, n_asv)
            te_loss = test_model(model, test_ld, criterion, n_asv)
            train_losses.append(tr_loss);
            test_losses.append(te_loss)

            # optional LR decay
            for g in optimiser.param_groups:
                g['lr'] *= lr_decay

        torch.save(model.state_dict(),
                   f"./mlp_models/collapse_classification/{sample}_{n_asv}ASV_fold{i+1}_threshold{threshold}_trial{trial}.pth")

        all_train_losses.append(train_losses)
        all_test_losses.append(test_losses)
    all_train_losses=np.array(all_train_losses)
    all_test_losses=np.array(all_test_losses)

    return all_train_losses, all_test_losses

if __name__=="__main__":
    sample = str(sys.argv[1])
    n_asv = int(sys.argv[2])
    threshold = float(sys.argv[3])
    trial = int(sys.argv[4])
    t0 = timer.perf_counter()

    train_losses, test_losses = run_collapse_cv(sample,n_asv,threshold,trial)

    loss_data = [train_losses, test_losses]
    filename = f"./mlp_models/collapse_classification/{sample}_{n_asv}ASV_threshold{threshold}_trial{trial}_traindoc.pth"
    with open(filename, 'wb') as f:
        pickle.dump(loss_data, f, protocol=pickle.HIGHEST_PROTOCOL)
    train_time = timer.perf_counter() - t0
    print(f"{sample} {n_asv}-ASV collapse forecasting MLP training finished, time used: {int(train_time)}s")