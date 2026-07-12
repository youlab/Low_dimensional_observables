import numpy as np
import torch
from torch.utils.data import DataLoader
from torch.nn import functional as F
from torch import nn
from sklearn.model_selection import train_test_split


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


def vae_loss(x_hat, x, mu, log_var):
    "Compute the sum of BCE and KL loss for the distribution."

    # weight for the KL divergence
    alpha = 1e-4

    # Compute the reconstruction loss
    BCE = F.mse_loss(x_hat, x)

    # Compute the KL divergence of the distribution.
    KLD = -0.5 * torch.mean(1 + log_var - mu.pow(2) - log_var.exp())

    return BCE + alpha * KLD

def train_model(model, data_loader, optimizer, device = "cuda:0"):
    num_batches = len(data_loader)
    total_loss = 0
    model.train()

    for batch in data_loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        pred, code, mu, log_var = model(batch)
        loss = vae_loss(pred, batch, mu, log_var)
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    avg_loss = total_loss / num_batches
    return avg_loss

def test_model(model, data_loader, device = "cuda:0"):

    num_batches = len(data_loader)
    total_loss = 0

    model.eval()
    with torch.no_grad():
        for batch in data_loader:
            batch = batch.to(device)
            pred, code, mu, log_var = model(batch)
            loss = vae_loss(pred, batch, mu, log_var)
            total_loss += loss.item()

    avg_loss = total_loss / num_batches
    return avg_loss

def get_data(n_target,test_idx):
    arr = np.arange(1,25,1)
    excluded = np.array([test_idx])
    train_arr = np.setdiff1d(arr, excluded)
    X_train = []
    for i in train_arr:
        sub = np.load("./vaginal_dataset_CV/S%i_T%i.npy" % (i, n_target))
        X_train.append(sub)
    X_train = np.vstack(X_train)
    X_test = np.load("./vaginal_dataset_CV/S%i_T%i.npy" % (test_idx, n_target))
    return X_train, X_test

def run(n_target,n_embedding,X_train,X_test,sub_id,trial):
    # set the random seed so the three trials have different initial model weights
    seed = 1000 * n_target + 100 * n_embedding + 10*sub_id + trial
    torch.manual_seed(seed)
    lr = 5e-4
    lr_decay = 0.95
    batch_size = 64
    EPOCHS = 30
    X_train = torch.Tensor(X_train).float()
    X_test = torch.Tensor(X_test).float()
    train_loader = DataLoader(X_train,batch_size=batch_size,shuffle=True)
    test_loader = DataLoader(X_test,batch_size=batch_size,shuffle=False)

    model = VAE(n_target,n_embedding)
    model.to("cuda:0");
    model.train()
    train_losses=[]
    test_losses=[]
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    for ix_epoch in range(EPOCHS):
        train_err=train_model(model, train_loader, optimizer)
        test_err=test_model(model, test_loader)
        train_losses.append(train_err)
        test_losses.append(test_err)
        # Exponential decay for learning rate
        lr *= lr_decay
        for param_group in optimizer.param_groups:
            param_group['lr'] = lr
    print("%i-dim, %i-target model training finished, with starting MSE %1.1e, and ending error %1.1e"%(n_embedding,n_target,test_losses[0],test_losses[-1]))
    torch.save(model.state_dict(),
           "./vae_models_CV/T%i_E%i_S%i_trial%i.pth"%(n_target,n_embedding,sub_id,trial))
    return model,train_losses, test_losses