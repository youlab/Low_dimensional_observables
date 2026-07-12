import numpy as np
import torch
from torch.utils.data import DataLoader
from torch.nn import functional as F
from torch import nn
from sklearn.model_selection import StratifiedKFold

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
        self.T = 168
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


class MLP_VAE(torch.nn.Module):

    def __init__(self, VAE, n_embedding, N):
        super().__init__()

        self.channels = 256

        self.mlp = nn.Sequential(
            nn.Linear(N, self.channels),
            nn.LeakyReLU(),
            nn.Linear(self.channels, self.channels),
            nn.LeakyReLU(),
            nn.Linear(self.channels, self.channels),
            nn.LeakyReLU(),
            nn.Linear(self.channels, n_embedding),
            nn.LeakyReLU(),
        )

        self.fc_decoder = VAE.fc_decoder.eval()
        self.decoder = VAE.decoder.eval()

        # Freeze the parameters of the VAE
        for param in self.fc_decoder.parameters():
            param.requires_grad = False
        for param in self.decoder.parameters():
            param.requires_grad = False

    def forward(self, X):

        B = X.shape[0]
        C = 32
        L = 168

        code = self.mlp(X)

        # Pass through FC layers before decoding
        post_code = self.fc_decoder(code)

        X_hat = self.decoder(post_code.view(B, C, L))

        return X_hat


class ABDCurvesDataSet(torch.utils.data.Dataset):

    def __init__(self, abd, curves):
        self.X = abd
        self.Y = curves

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return self.X[idx], self.Y[idx]


def get_abd_ts_data(target):
    idx_FP_map = {"EGFP": 0, "mTagBFP2": 1, "LSSmOrange": 2, "mCherry": 3,
                  "OD": 4}  # correspondence between the FP and their position in the time series data

    id_FP = idx_FP_map[target]
    Y = np.load("./sequenced_data/sequenced_target_normal.npy")[:, id_FP:id_FP + 1, :]

    composition = np.loadtxt("./sequenced_data/sequence_composition.txt")
    donor_abd = composition[:, 0:4]

    donor_abd = donor_abd / np.max(donor_abd, axis=0)

    Y = torch.Tensor(Y).float()
    donor_abd = torch.Tensor(donor_abd[:, np.newaxis]).float()

    full_dataset = ABDCurvesDataSet(donor_abd, Y)
    return full_dataset


# train and test function for mlp
def train_model(model, data_loader, optimizer):
    num_batches = len(data_loader)
    total_loss = 0
    model.mlp.train()

    for abd, curves in data_loader:
        abd = abd.to(device)
        curves = curves.to(device)

        optimizer.zero_grad()
        pred = model(abd)
        loss = F.mse_loss(pred, curves)
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    avg_loss = total_loss / num_batches
    return avg_loss


def test_model(model, data_loader):
    num_batches = len(data_loader)
    total_loss = 0

    model.mlp.eval()
    with torch.no_grad():
        for abd, curves in data_loader:
            abd = abd.to(device)
            curves = curves.to(device)
            pred = model(abd)
            loss = F.mse_loss(pred, curves)
            total_loss += loss.item()

    avg_loss = total_loss / num_batches
    return avg_loss

def run_cross_validation_growth(target,trial):
    seed = trial*42

    lr_0 = 3e-4

    EPOCHS = 100
    lr_decay = 0.99
    batch_size = 16

    n_embedding = 1
    n_target = 1
    # Load pretrained VAE
    vae_model = VAE(n_target, n_embedding)
    vae_model.load_state_dict(
        torch.load(
            f"./vae_models/{target}_E1.pth",
            weights_only=True
        )
    )
    vae_model.to("cuda:0")
    vae_model.eval()

    full_dataset = get_abd_ts_data(target)

    all_train_losses = []
    all_test_losses = []

    # stratified 10-fold CV based on maximum FP reading
    k_folds = 10
    y_np = full_dataset.Y.numpy()
    y_cont = np.max(y_np, axis=(1, 2))  # 1-D numpy array or list
    n_bins = 5  # 5 quantile bins  →  roughly equal counts
    bins = np.quantile(y_cont, np.linspace(0, 1, n_bins + 1))
    y_bins = np.digitize(y_cont, bins[1:-1])

    skf = StratifiedKFold(n_splits=k_folds, shuffle=True, random_state=seed)

    for fold, (train_idx, test_idx) in enumerate(skf.split(np.zeros(len(y_bins)), y_bins)):
        torch.manual_seed(seed+fold)
        lr = lr_0 * 1

        train_subset = torch.utils.data.Subset(full_dataset, train_idx)
        test_subset = torch.utils.data.Subset(full_dataset, test_idx)

        train_loader = DataLoader(train_subset, batch_size=batch_size, shuffle=True)
        test_loader = DataLoader(test_subset, batch_size=batch_size, shuffle=False)

        mlpvae = MLP_VAE(vae_model, n_embedding, 4).float().to(device)
        # print("MLP-VAE model parameters = %1.1e" % count_parameters(mlpvae))
        mlpvae.mlp.train()

        optimizer = torch.optim.Adam(mlpvae.mlp.parameters(), lr=lr)
        train_losses = []
        test_losses = []

        test_err = test_model(mlpvae, test_loader)
        # print(f"Initial test loss: {test_err:.2e}")

        for ix_epoch in range(EPOCHS):
            train_err = train_model(mlpvae, train_loader, optimizer)
            test_err = test_model(mlpvae, test_loader)
            train_losses.append(train_err)
            test_losses.append(test_err)
            # Exponential decay for learning rate
            lr *= lr_decay
            for param_group in optimizer.param_groups:
                param_group['lr'] = lr

        all_train_losses.append(train_losses)
        all_test_losses.append(test_losses)

        torch.save(mlpvae.state_dict(), f"./mlp_models/{target}/{target}_growth_curve_inference_{fold + 1}-fold_trial{trial}.pth")

    return all_train_losses, all_test_losses