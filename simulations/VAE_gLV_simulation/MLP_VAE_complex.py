import numpy as np
import torch
from torch.utils.data import DataLoader
from torch.nn import functional as F
from torch import nn
import sys
import time as timer

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


class MLP_VAE(torch.nn.Module):

    def __init__(self, VAE, n_embedding, N):
        super().__init__()

        self.channels = 128

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
        L = 50

        code = self.mlp(X)

        # Pass through FC layers before decoding
        post_code = self.fc_decoder(code)

        X_hat = self.decoder(post_code.view(B, C, L))

        return X_hat


class IVCurveDataSet(torch.utils.data.Dataset):

    def __init__(self, ivs, curves):
        self.X = ivs
        self.Y = curves

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return self.X[idx], self.Y[idx]


def get_data_100K(model_index, n_ori, n_target):
    n_background = n_ori - n_target

    X_train = np.load("./saved_sims/bgLV/I%i/bgLV_random_init_train_100K.npy"
                      % (model_index))
    X_test = np.load("./saved_sims/bgLV/I%i/bgLV_random_init_test_100K.npy"
                     % (model_index))
    X_train = torch.Tensor(X_train).float() * 5
    X_test = torch.Tensor(X_test).float() * 5

    y_train = np.load("./saved_sims/bgLV/I%i/bgLV_B%i_T%i_random_train_100K.npy"
                      % (model_index, n_background, n_target))
    y_test = np.load("./saved_sims/bgLV/I%i/bgLV_B%i_T%i_random_test_100K.npy"
                     % (model_index, n_background, n_target))

    y_train = torch.Tensor(y_train).float()
    y_test = torch.Tensor(y_test).float()

    data_train = IVCurveDataSet(X_train, y_train)
    data_test = IVCurveDataSet(X_test, y_test)

    return data_train, data_test

#train and test function for mlp probe
def train_model(model, data_loader, optimizer):
    num_batches = len(data_loader)
    total_loss = 0
    model.mlp.train()

    for iv,curves in data_loader:
        iv = iv.to("cuda:0")
        curves = curves.to("cuda:0")
        optimizer.zero_grad()
        pred = model(iv)
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
        for iv,curves in data_loader:
            iv = iv.to("cuda:0")
            curves = curves.to("cuda:0")
            pred = model(iv)
            loss = F.mse_loss(pred, curves)
            total_loss += loss.item()

    avg_loss = total_loss / num_batches
    return avg_loss

def train_VAEMLP(model_index,n_ori,n_target,n_embedding,X_train,X_test,trial):
    t0 = timer.perf_counter()
    # set the random seed so the three trials have different initial model weights
    seed = 1000 * n_target + 10*n_embedding + trial
    torch.manual_seed(seed)

    EPOCHS = 100
    lr = 2e-3
    lr_decay = 0.99
    batch_size = 64

    n_background = n_ori-n_target

    vae_model = VAE(n_target, n_embedding)
    vae_model.to("cuda:0");
    vae_model.load_state_dict(
        torch.load(
            "./vae_models/bgLV/I%i/random_B%i_T%i_E%i_retrain.pth"
            % (model_index, n_background, n_target, n_embedding), weights_only=True)
    )
    vae_model.eval()

    mlpvae = MLP_VAE(vae_model, n_embedding, n_ori)
    mlpvae.to("cuda:0");
    mlpvae.mlp.train()

    train_loader = DataLoader(X_train, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(X_test, batch_size=batch_size, shuffle=False)
    train_losses = []
    test_losses = []


    optimizer = torch.optim.Adam(mlpvae.mlp.parameters(), lr=lr)

    for ix_epoch in range(EPOCHS):
        train_err = train_model(mlpvae, train_loader, optimizer)
        test_err = test_model(mlpvae, test_loader)
        train_losses.append(train_err)
        test_losses.append(test_err)
        # Exponential decay for learning rate
        lr *= lr_decay
        for param_group in optimizer.param_groups:
            param_group['lr'] = lr

    torch.save(mlpvae.state_dict(),
               "./mlp_models/bgLV/I%i/MLP_T%i_E%i_trial%i.pth"
               % (model_index, n_target, n_embedding, trial)
               )
    train_time = timer.perf_counter() - t0
    print("training finished, with starting MSE %1.1e, and ending error %1.1e, training time %1.1f s" % (test_losses[0], test_losses[-1], train_time))
    return train_losses, test_losses

model_index = int(sys.argv[1])

n_target = int(sys.argv[2])

trial = int(sys.argv[3])

n_ori = 100

X_train, X_test = get_data_100K(model_index, n_ori, n_target)

Ec_data = np.loadtxt(f"./saved_data/bgLV/Ec_bgLV_model{model_index}_random.txt")
N_TG = Ec_data[0,:]
Ec = Ec_data[1,:]
k = np.where(N_TG==n_target)[0][0]
n_embedding = int(np.ceil(Ec[k]))

train_loss, test_loss = train_VAEMLP(model_index,n_ori,n_target,n_embedding,X_train,X_test,trial)

