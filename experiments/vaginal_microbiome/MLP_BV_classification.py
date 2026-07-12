import numpy as np
import torch
from torch.utils.data import DataLoader,TensorDataset
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

        self.mean_map = nn.Linear(self.channels * 14, self.latent_dim)

        self.std_map = nn.Linear(self.channels * 14, self.latent_dim)

        self.fc_decoder = nn.Sequential(
            nn.Linear(self.latent_dim, self.channels * 14),
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


class VAE_MLP_Classifier(nn.Module):
    def __init__(self, VAE, n_embedding, n_asv, hidden_dim=32, n_classes=3):
        super().__init__()

        self.encoder = VAE.encoder.eval()
        self.mean_map = VAE.mean_map.eval()
        self.std_map = VAE.std_map.eval()

        self.mlp = nn.Sequential(
            nn.Linear(n_asv * n_embedding, hidden_dim),
            nn.LeakyReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LeakyReLU(),
            nn.Linear(hidden_dim, n_classes)  # logits for 3 classes
        )

        # Freeze VAE parameters
        for param in self.encoder.parameters():
            param.requires_grad = False
        for param in self.mean_map.parameters():
            param.requires_grad = False
        for param in self.std_map.parameters():
            param.requires_grad = False

    def sample(self, mean, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return mean + std * eps

    def forward(self, Xs, n_asv):
        embeddings = []
        for i in range(n_asv):
            X = Xs[:, i:i + 1, :]  # shape: (B, 1, T)
            pre_code = self.encoder(X)
            B, C, L = pre_code.shape
            flat = pre_code.view(B, C * L)

            mu = self.mean_map(flat)
            log_var = self.std_map(flat)
            z = self.sample(mu, log_var)
            embeddings.append(z)

        full_embedding = torch.cat(embeddings, dim=1)  # shape (B, n_asv * n_embedding)
        logits = self.mlp(full_embedding)  # shape (B, 3)

        return logits

def train_model(model, data_loader, optimizer, criterion, n_asv):
    model.train()
    total_loss = 0.0

    for x, y in data_loader:
        x = x.to(device)
        y = y.to(device).long()  # class labels as integers 0/1/2

        optimizer.zero_grad()
        logits = model(x, n_asv)  # shape: (B, 3)
        loss = criterion(logits, y)  # CrossEntropyLoss
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    return total_loss / len(data_loader)

def test_model(model, data_loader, criterion, n_asv):
    model.eval()
    total_loss = 0.0

    with torch.no_grad():
        for x, y in data_loader:
            x = x.to(device)
            y = y.to(device).long()

            logits = model(x, n_asv)
            loss = criterion(logits, y)
            total_loss += loss.item()

    return total_loss / len(data_loader)

def nugent_to_class(score):
    if score <= 3:
        return 0  # Normal
    elif score <= 6:
        return 1  # Intermediate
    else:
        return 2  # BV


def get_BV_next(n_target):
    # Load meta data
    composition_compiled = np.load(
        "./meta_data/microbiome_timeseries_ranked.npz")
    ngscore_compiled = np.load(
        "./meta_data/ngscore_timeseries.npz")

    subjects = composition_compiled.files
    subjects.remove("subject_49")

    Seg_all = {}
    BV_all = {}
    T = 14

    for name in subjects:
        ngscore_t = ngscore_compiled[name]
        composition_data = composition_compiled[name]
        Segments = []
        BV = []
        for start in range(composition_data.shape[0] - T - 1):
            seg = composition_data[start:start + T, 0:n_target]
            ngscore = ngscore_t[start + T]
            if np.isnan(ngscore):
                continue
            Segments.append(seg.T)
            BV.append(nugent_to_class(ngscore))
        BV = torch.Tensor(BV)
        Segments = torch.Tensor(np.array(Segments)).float()  # (batch, n_target, T)

        Seg_all[name] = Segments
        BV_all[name] = BV

    return Seg_all, BV_all

def run_menses_classification_cv(n_target,trial):
    torch.manual_seed(trial)
    lr0      = 1e-3
    EPOCHS   = 40
    batch_sz = 16
    lr_decay = 0.98

    # ── Load data ───────────────────────────────────────────────
    Seg_all, BV_all = get_BV_next(n_target)
    subjects = sorted(Seg_all.keys())

    Ec = 8
    vae = VAE(1, Ec)
    vae.load_state_dict(torch.load(f"./vae_models_CV/T1_E{Ec}_full.pth", weights_only=True))
    vae.eval().to(device)

    all_train_losses, all_test_losses = [], []

    for i, subject in enumerate(subjects):
        # Split
        X_test  = Seg_all[subject]                         # shape: (n_segs, n_target, T)
        y_test  = BV_all[subject].long()                   # shape: (n_segs,)
        X_train = torch.cat([Seg_all[s] for j,s in enumerate(subjects) if j != i])
        y_train = torch.cat([BV_all[s] for j,s in enumerate(subjects) if j != i]).long()

        train_ds = TensorDataset(X_train, y_train)
        test_ds  = TensorDataset(X_test, y_test)
        train_ld = DataLoader(train_ds, batch_size=batch_sz, shuffle=True)
        test_ld  = DataLoader(test_ds, batch_size=batch_sz, shuffle=False)

        # Class balance weights
        n_classes = 3
        class_counts = torch.bincount(y_train, minlength=n_classes).float()
        class_weights = (1.0 / class_counts).clamp(max=10)
        class_weights = class_weights / class_weights.sum() * n_classes
        criterion = torch.nn.CrossEntropyLoss(weight=class_weights.to(device))

        # Model
        model = VAE_MLP_Classifier(vae, Ec, n_target, hidden_dim=32, n_classes=3).to(device)
        optimiser = torch.optim.Adam(model.mlp.parameters(), lr=lr0)

        train_losses, test_losses = [], []
        for epoch in range(EPOCHS):
            tr_loss = train_model(model, train_ld, optimiser, criterion, n_target)
            te_loss = test_model(model, test_ld, criterion, n_target)
            train_losses.append(tr_loss)
            test_losses.append(te_loss)

            for g in optimiser.param_groups:
                g['lr'] *= lr_decay

        torch.save(model.state_dict(),
                   f"./mlp_models_CV/BV_forecast/{n_target}ASV_fold{i+1}_trial{trial}.pth")

        all_train_losses.append(train_losses)
        all_test_losses.append(test_losses)
    all_train_losses = np.array(all_train_losses)
    all_test_losses = np.array(all_test_losses)
    print(f"{n_target} classification result trial {trial}: train loss = {np.mean(all_train_losses[:,-1]):.3f} ± {np.std(all_train_losses[:,-1]):.3f} ")
    print(f"test loss = {np.mean(all_test_losses[:,-1]):.3f} ± {np.std(all_test_losses[:,-1]):.3f} ")
    return all_train_losses, all_test_losses

if __name__=="__main__":
    n_target = int(sys.argv[1])
    trial = int(sys.argv[2])
    t0 = timer.perf_counter()
    train_losses, test_losses = run_menses_classification_cv(n_target, trial)
    loss_data = [train_losses,test_losses]
    filename = f"./mlp_models_CV/BV_forecast/{n_target}ASV_trial{trial}_traindoc.pkl"
    with open(filename, 'wb') as f:
        pickle.dump(loss_data, f, protocol=pickle.HIGHEST_PROTOCOL)
    train_time = timer.perf_counter() - t0
    print(f"{n_target} taxa BV-forecasting MLP training finished, time used: {int(train_time)}s")
