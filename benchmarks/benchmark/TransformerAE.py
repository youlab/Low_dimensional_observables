import numpy as np
import torch
from torch.utils.data import DataLoader
from torch.nn import functional as F
from torch import nn
from pathlib import Path
import pickle

device = "cuda:0"

class PositionalEncoding(nn.Module):
    def __init__(self, d_model: int, max_len: int = 50):
        super().__init__()

        pe = torch.zeros(max_len, d_model)  # (T, d_model)
        position = torch.arange(0, max_len, dtype=torch.float32).unsqueeze(1)  # (T, 1)
        div_term = torch.exp(
            torch.arange(0, d_model, 2, dtype=torch.float32)
            * (-torch.log(torch.tensor(10000.0)) / d_model)
        )

        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term[: pe[:, 1::2].shape[1]])
        pe = pe.unsqueeze(0)  # (1, T, d_model)

        self.register_buffer("pe", pe)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        x: (B, T, d_model)
        """
        return x + self.pe[:, :x.size(1)]


class TransformerAutoEncoder(nn.Module):
    """
    Transformer autoencoder for input of shape (B, n_target, T),
    treating each timepoint vector (n_target,) as one token.
    """
    def __init__(
        self,
        n_target: int,
        latent_dim: int,
        T: int = 50,
        d_model: int = 32,
        nhead: int = 4,
        num_encoder_layers: int = 2,
        num_decoder_layers: int = 2,
        dim_feedforward: int = 32,
    ):
        super().__init__()

        self.n_target = n_target
        self.T = T
        self.d_model = d_model
        self.latent_dim = latent_dim

        # Token embedding
        self.input_proj = nn.Linear(n_target, d_model)
        self.pos_encoder = PositionalEncoding(d_model, max_len=T)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=0.0,
            activation=nn.LeakyReLU(),
            batch_first=True,
        )
        self.encoder = nn.TransformerEncoder(
            encoder_layer,
            num_layers=num_encoder_layers,
            enable_nested_tensor=False,
        )

        # Latent bottleneck
        self.to_latent = nn.Linear(T * d_model, latent_dim)

        self.from_latent = nn.Linear(latent_dim, T * d_model)

        self.pos_decoder = PositionalEncoding(d_model, max_len=T)

        decoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=0.0,
            activation=nn.LeakyReLU(),
            batch_first=True,
        )
        self.decoder = nn.TransformerEncoder(
            decoder_layer,
            num_layers=num_decoder_layers,
            enable_nested_tensor=False,
        )

        self.output_proj = nn.Linear(d_model, n_target)

    def encode_tokens(self, x: torch.Tensor) -> torch.Tensor:
        """
        x: (B, n_target, T)
        returns token features: (B, T, d_model)
        """
        x = x.transpose(1, 2)   # (B, T, n_target)
        x = self.input_proj(x)  # (B, T, d_model)
        x = self.pos_encoder(x)
        x = self.encoder(x)
        return x

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """
        x: (B, n_target, T)
        returns latent z: (B, latent_dim)
        """
        h = self.encode_tokens(x)         # (B, T, d_model)
        h_flat = h.reshape(h.size(0), -1) # (B, T*d_model)
        z = self.to_latent(h_flat)        # (B, latent_dim)
        return z

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """
        z: (B, latent_dim)
        returns reconstruction: (B, n_target, T)
        """
        h = self.from_latent(z)                             # (B, T*d_model)
        h = h.view(z.size(0), self.T, self.d_model)        # (B, T, d_model)
        h = self.pos_decoder(h)
        h = self.decoder(h)
        out = self.output_proj(h)                          # (B, T, n_target)
        out = out.transpose(1, 2)                          # (B, n_target, T)
        return out

    def forward(self, x: torch.Tensor):
        """
        x: (B, n_target, T)
        returns:
            x_recon: (B, n_target, T)
            z:       (B, latent_dim)
        """
        z = self.encode(x)
        x_recon = self.decode(z)
        return x_recon, z
    
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

    model = TransformerAutoEncoder(n_target, n_embedding)
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
        f"./transformer/{glv_type}gLV/I{model_index}/{sim_type}_B{n_background}_T{n_target}_E{n_embedding}_trial{trial}.pth"
    )
    save_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), save_path)

    print("training finished, with starting MSE %1.1e, and ending error %1.1e" % (test_losses[0], test_losses[-1]))
    return train_losses, test_losses