import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from scipy.stats import pearsonr
import shap
import matplotlib.pyplot as plt

# Load your consolidated CSV
df = pd.read_csv("combined_df.csv")

# Drop missing values and shuffle
df = df.dropna().sample(frac=1, random_state=42).reset_index(drop=True)

# Select omics data ranges without applying log transformation
omics_ranges = {
    "Metabolomics": (12, 268),
    "MGS": (2790, 3043),
    "KEGG_1": (282, 2575),
    "KEGG_2": (2575, 2740),
    "KEGG_3": (2740, 2790),
    "GMM": (268, 282),
}
for key, (start, end) in omics_ranges.items():
    df.iloc[:, start:end] = df.iloc[:, start:end].astype(np.float64)

# First split: 60% train, 40% temp
train_df, temp_df = train_test_split(df, test_size=0.4, random_state=42)

# Second split: 50% of temp = 20% of total for val and test
val_df, test_df = train_test_split(temp_df, test_size=0.5, random_state=42)
# Save validation set for consistent reuse
val_df.to_csv("val_df.csv", index=False)


# Define column groups
covariate_cols = df.columns[2:7].tolist()  
x1_cols = [df.columns[8]]
x2_cols = [df.columns[9]]         
y1_col = [df.columns[10]]          
y2_col = [df.columns[11]]         

# Extract tensors
def df_to_tensor(df, cols):
    return torch.tensor(df[cols].values, dtype=torch.float32)

# Training tensors
omics_train = {key: df_to_tensor(train_df, df.columns[start:end]) for key, (start, end) in omics_ranges.items()}
x1_train = df_to_tensor(train_df, x1_cols)
x2_train = df_to_tensor(train_df, x2_cols)
y1_train = df_to_tensor(train_df, y1_col)
y2_train = df_to_tensor(train_df, y2_col)
covariates_train = df_to_tensor(train_df, covariate_cols)

# Validation tensors
omics_val = {key: df_to_tensor(val_df, df.columns[start:end]) for key, (start, end) in omics_ranges.items()}
x1_val = df_to_tensor(val_df, x1_cols)
x2_val = df_to_tensor(val_df, x2_cols)
y1_val = df_to_tensor(val_df, y1_col)
y2_val = df_to_tensor(val_df, y2_col)
covariates_val = df_to_tensor(val_df, covariate_cols)

# Test tensors
omics_test = {key: df_to_tensor(test_df, df.columns[start:end]) for key, (start, end) in omics_ranges.items()}
x1_test = df_to_tensor(test_df, x1_cols)
x2_test = df_to_tensor(test_df, x2_cols)
y1_test = df_to_tensor(test_df, y1_col)
y2_test = df_to_tensor(test_df, y2_col)
covariates_test = df_to_tensor(test_df, covariate_cols)

# Custom Dataset
class MultimodalDataset(Dataset):
    def __init__(self, omics, x1, x2, y1, y2, covariates):
        self.omics = omics
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.covariates = covariates

    def __len__(self):
        return len(next(iter(self.omics.values())))

    def __getitem__(self, idx):
        omics_data = {key: value[idx] for key, value in self.omics.items()}
        return omics_data, self.x1[idx], self.x2[idx], self.y1[idx], self.y2[idx], self.covariates[idx]

# DataLoaders
train_dataset = MultimodalDataset(omics_train, x1_train, x2_train, y1_train, y2_train, covariates_train)
val_dataset = MultimodalDataset(omics_val, x1_val, x2_val, y1_val, y2_val, covariates_val)
test_dataset = MultimodalDataset(omics_test, x1_test, x2_test, y1_test, y2_test, covariates_test)

train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False)
test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

# Omics Encoder with VAE-style latent sampling
class OmicsEncoder(nn.Module):
    def __init__(self, input_dim, latent_dim):
        super(OmicsEncoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 128),
            nn.ReLU(),
            nn.Linear(128, latent_dim * 2)  # outputs mean and log-variance
        )

    def forward(self, x):
        encoded = self.encoder(x)
        mean, log_var = encoded[:, :latent_dim], encoded[:, latent_dim:]
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        z = mean + eps * std  # reparameterization trick
        return z, mean, log_var

# Attention mechanism over modalities
class Attention(nn.Module):
    def __init__(self, latent_dim, num_modalities):
        super(Attention, self).__init__()
        self.attention = nn.Sequential(
            nn.Linear(latent_dim * num_modalities, 128),
            nn.ReLU(),
            nn.Linear(128, num_modalities),
            nn.Softmax(dim=1)
        )

    def forward(self, z):
        weights = self.attention(z)
        weighted_z = torch.stack([
            weights[:, i].unsqueeze(1) * z[:, i * latent_dim:(i + 1) * latent_dim]
            for i in range(weights.size(1))
        ], dim=1)
        return weighted_z.sum(dim=1)

# Y1 Predictor
class Y1Predictor(nn.Module):
    def __init__(self, latent_dim, x1_dim, x2_dim, covariate_dim, output_dim):
        super(Y1Predictor, self).__init__()
        self.predictor = nn.Sequential(
            nn.Linear(latent_dim + x1_dim + x2_dim + covariate_dim, 128),
            nn.ReLU(),
            nn.Linear(128, output_dim)
        )

    def forward(self, m, x1, x2, covariates):
        combined = torch.cat((m, x1, x2, covariates), dim=1)
        return self.predictor(combined)

# Y2 Predictor
class Y2Predictor(nn.Module):
    def __init__(self, latent_dim, x1_dim, x2_dim, y1_dim, covariate_dim, output_dim):
        super(Y2Predictor, self).__init__()
        self.predictor = nn.Sequential(
            nn.Linear(latent_dim + x1_dim + x2_dim + y1_dim + covariate_dim, 128),
            nn.ReLU(),
            nn.Linear(128, output_dim)
        )

    def forward(self, m, x1, x2, y1, covariates):
        combined = torch.cat((m, x1, x2, y1, covariates), dim=1)
        return self.predictor(combined)

# Full Model with VAE + Attention
class MultimodalModel(nn.Module):
    def __init__(self, omics_dims, latent_dim, x1_dim, x2_dim, y1_dim, y2_dim, covariate_dim):
        super(MultimodalModel, self).__init__()
        self.encoders = nn.ModuleDict({
            key: OmicsEncoder(dim, latent_dim) for key, dim in omics_dims.items()
        })
        self.attention = Attention(latent_dim, len(omics_dims))
        self.y1_predictor = Y1Predictor(latent_dim, x1_dim, x2_dim, covariate_dim, y1_dim)
        self.y2_predictor = Y2Predictor(latent_dim, x1_dim, x2_dim, y1_dim, covariate_dim, y2_dim)

    def forward(self, omics, x1, x2, covariates):
        z, means, log_vars = [], [], []
        for key, encoder in self.encoders.items():
            latent, mean, log_var = encoder(omics[key])
            z.append(latent)
            means.append(mean)
            log_vars.append(log_var)
        z = torch.cat(z, dim=1)
        z = self.attention(z)
        y1 = self.y1_predictor(z, x1, x2, covariates)
        y2 = self.y2_predictor(z, x1, x2, y1, covariates)
        return y1, y2, z, means, log_vars

# KL divergence for VAE loss
def kl_divergence(mean, log_var):
    return -0.5 * torch.sum(1 + log_var - mean.pow(2) - log_var.exp())

# Initialize model
omics_dims = {key: end - start for key, (start, end) in omics_ranges.items()}
latent_dim = 64
x1_dim = len(x1_cols)
x2_dim = len(x2_cols)
y1_dim = len(y1_col)
y2_dim = len(y2_col)
covariate_dim = len(covariate_cols)

model = MultimodalModel(omics_dims, latent_dim, x1_dim, x2_dim, y1_dim, y2_dim, covariate_dim)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

# Loss and optimizer
criterion = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

# Evaluate model performance
def evaluate(model, dataloader, criterion):
    model.eval()
    total_loss = 0
    with torch.no_grad():
        for omics, x1, x2, y1_true, y2_true, covariates in dataloader:
            y1_pred, y2_pred, _, _, _ = model(omics, x1, x2, covariates)
            loss = criterion(y1_pred, y1_true) + criterion(y2_pred, y2_true)
            total_loss += loss.item()
    return total_loss / len(dataloader)

# Training loop with optional KL loss
def train(model, train_loader, val_loader, optimizer, criterion, epochs=10, beta=0.001):
    model.train()
    for epoch in range(epochs):
        train_loss = 0
        for omics, x1, x2, y1_true, y2_true, covariates in train_loader:
            optimizer.zero_grad()
            y1_pred, y2_pred, _, means, log_vars = model(omics, x1, x2, covariates)
            recon_loss = criterion(y1_pred, y1_true) + criterion(y2_pred, y2_true)
            kl_loss = sum(kl_divergence(m, lv) for m, lv in zip(means, log_vars)) / len(means)
            loss = recon_loss + beta * kl_loss
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
        avg_train_loss = train_loss / len(train_loader)
        val_loss = evaluate(model, val_loader, criterion)
        print(f"Epoch {epoch+1}, Training Loss: {avg_train_loss:.4f}, Validation Loss: {val_loss:.4f}")
# Train the model
train(model, train_loader, val_loader, optimizer, criterion, epochs=10)
model.eval()

# Export latent space
def export_latent_space(model, dataloader, filename="latent_space.csv"):
    model.eval()
    latent_space = []
    for omics, x1, x2, y1_true, y2_true, covariates in dataloader:
        omics = {k: v.to(device) for k, v in omics.items()}
        x1, x2, covariates = x1.to(device), x2.to(device), covariates.to(device)
        with torch.no_grad():
            _, _, m, _, _ = model(omics, x1, x2, covariates)
        latent_space.append(m.cpu().numpy())
    latent_space = np.concatenate(latent_space, axis=0)
    pd.DataFrame(latent_space).to_csv(filename, index=False)

export_latent_space(model, val_loader)

# Load latent space and Y1/Y2
latent_space = pd.read_csv("latent_space.csv")
val_df = pd.read_csv("val_df.csv")
y1_val = val_df.iloc[:, 10].values
y2_val = val_df.iloc[:, 11].values

def rank_latent_dimensions(latent_space, y1, y2):
    correlations_y1 = [pearsonr(latent_space[:, i], y1)[0] for i in range(latent_space.shape[1])]
    correlations_y2 = [pearsonr(latent_space[:, i], y2)[0] for i in range(latent_space.shape[1])]
    importance_df = pd.DataFrame({
        'Latent Dimension': range(latent_space.shape[1]),
        'Correlation with Y1': correlations_y1,
        'Correlation with Y2': correlations_y2,
        'Absolute Correlation with Y1': np.abs(correlations_y1),
        'Absolute Correlation with Y2': np.abs(correlations_y2)
    })
    return importance_df.sort_values(by=['Absolute Correlation with Y1', 'Absolute Correlation with Y2'], ascending=False)

importance_df = rank_latent_dimensions(latent_space.values, y1_val, y2_val)
importance_df.to_csv("latent_importance_ranking.csv", index=False)

# Prepare input matrix for SHAP
def prepare_input_matrix(dataloader):
    inputs = []
    for omics, x1, x2, y1_true, y2_true, covariates in dataloader:
        omics = {k: v.to(device) for k, v in omics.items()}
        x1, x2, covariates = x1.to(device), x2.to(device), covariates.to(device)
        with torch.no_grad():
            _, _, m, _, _ = model(omics, x1, x2, covariates)
        combined = torch.cat((m, x1, x2, covariates), dim=1)
        inputs.append(combined.cpu().numpy())
    return torch.tensor(np.concatenate(inputs, axis=0), dtype=torch.float32)

input_tensor = prepare_input_matrix(val_loader).to(device)

# SHAP Wrappers
class Y1Wrapper(nn.Module):
    def __init__(self, model):
        super().__init__()
        self.model = model

    def forward(self, x):
        m = x[:, :latent_dim]
        x1 = x[:, latent_dim:latent_dim + x1_dim]
        x2 = x[:, latent_dim + x1_dim:latent_dim + x1_dim + x2_dim]
        covariates = x[:, -covariate_dim:]
        return self.model.y1_predictor(m, x1, x2, covariates)

class Y2Wrapper(nn.Module):
    def __init__(self, model):
        super().__init__()
        self.model = model

    def forward(self, x):
        m = x[:, :latent_dim]
        x1 = x[:, latent_dim:latent_dim + x1_dim]
        x2 = x[:, latent_dim + x1_dim:latent_dim + x1_dim + x2_dim]
        covariates = x[:, -covariate_dim:]
        y1 = self.model.y1_predictor(m, x1, x2, covariates)
        return self.model.y2_predictor(m, x1, x2, y1, covariates)

# SHAP Analysis
wrapped_model_y1 = Y1Wrapper(model).to(device)
wrapped_model_y2 = Y2Wrapper(model).to(device)

explainer_y1 = shap.DeepExplainer(wrapped_model_y1, input_tensor)
shap_values_y1 = explainer_y1.shap_values(input_tensor)

explainer_y2 = shap.DeepExplainer(wrapped_model_y2, input_tensor)
shap_values_y2 = explainer_y2.shap_values(input_tensor)

shap.summary_plot(shap_values_y1, input_tensor.cpu().numpy(), show=False)
plt.savefig("shap_summary_plot_y1.png")
plt.show()

shap.summary_plot(shap_values_y2, input_tensor.cpu().numpy(), show=False)
plt.savefig("shap_summary_plot_y2.png")
plt.show()

np.save("shap_values_y1.npy", shap_values_y1)
np.save("shap_values_y2.npy", shap_values_y2)

# Final performance evaluation
train_loss = evaluate(model, train_loader, criterion)
val_loss = evaluate(model, val_loader, criterion)
test_loss = evaluate(model, test_loader, criterion)

performance_stats = {
    "Training Loss": train_loss,
    "Validation Loss": val_loss,
    "Test Loss": test_loss
}
pd.DataFrame([performance_stats]).to_csv("model_performance.csv", index=False)

print("Model training and evaluation complete. Latent space and performance statistics exported.")


