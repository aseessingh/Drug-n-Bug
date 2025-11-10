
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import numpy as np

# Load your consolidated CSV
df = pd.read_csv("combined_df.csv")

# Optional: check for missing values or shuffle
df = df.dropna().sample(frac=1, random_state=42).reset_index(drop=True)

# First split: 60% train, 40% temp
train_df, temp_df = train_test_split(df, test_size=0.4, random_state=42)

# Second split: 50% of temp = 20% of total for val and test
val_df, test_df = train_test_split(temp_df, test_size=0.5, random_state=42)

# Define column groups
omics_ranges = {
    "Metabolomics": (12, 268),
    "MGS": (2790,3043),
    "KEGG_1": (282,2575),
    "KEGG_2": (2575,2740),
    "KEGG_3": (2740,2790),
    "GMM": (268,282),
}
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

# Omics Encoder
class OmicsEncoder(nn.Module):
    def __init__(self, input_dim, latent_dim):
        super(OmicsEncoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 128),
            nn.ReLU(),
            nn.Linear(128, latent_dim)
        )

    def forward(self, x):
        return self.encoder(x)

# Y1 Predictor: uses M + X1 + X2 + Covariates
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

# Y2 Predictor: uses M + X1 + X2 + Y1 + Covariates
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

# Full Model
class MultimodalModel(nn.Module):
    def __init__(self, omics_dims, latent_dim, x1_dim, x2_dim, y1_dim, y2_dim, covariate_dim):
        super(MultimodalModel, self).__init__()
        self.encoders = nn.ModuleDict({key: OmicsEncoder(dim, latent_dim) for key, dim in omics_dims.items()})
        self.y1_predictor = Y1Predictor(latent_dim * len(omics_dims), x1_dim, x2_dim, covariate_dim, y1_dim)
        self.y2_predictor = Y2Predictor(latent_dim * len(omics_dims), x1_dim, x2_dim, y1_dim, covariate_dim, y2_dim)

    def forward(self, omics, x1, x2, covariates):
        m = torch.cat([encoder(omics[key]) for key, encoder in self.encoders.items()], dim=1)
        y1 = self.y1_predictor(m, x1, x2, covariates)
        y2 = self.y2_predictor(m, x1, x2, y1, covariates)
        return y1, y2, m

# Initialize model
omics_dims = {key: end - start for key, (start, end) in omics_ranges.items()}
latent_dim = 64
x1_dim = len(x1_cols)
x2_dim = len(x2_cols)
y1_dim = len(y1_col)
y2_dim = len(y2_col)
covariate_dim = len(covariate_cols)

model = MultimodalModel(omics_dims, latent_dim, x1_dim, x2_dim, y1_dim, y2_dim, covariate_dim)

# Loss and optimizer
criterion = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

# Training loop
def train(model, dataloader, optimizer, criterion, epochs=10):
    model.train()
    for epoch in range(epochs):
        for omics, x1, x2, y1_true, y2_true, covariates in dataloader:
            optimizer.zero_grad()
            y1_pred, y2_pred, _ = model(omics, x1, x2, covariates)
            loss = criterion(y1_pred, y1_true) + criterion(y2_pred, y2_true)
            loss.backward()
            optimizer.step()
        print(f"Epoch {epoch+1}, Loss: {loss.item():.4f}")

# Train the model
train(model, train_loader, optimizer, criterion, epochs=10)

# Export latent space
def export_latent_space(model, dataloader, filename="latent_space.csv"):
    model.eval()
    latent_space = []
    for omics, x1, x2, y1_true, y2_true, covariates in dataloader:
        _, _, m = model(omics, x1, x2, covariates)
        latent_space.append(m.detach().cpu().numpy())
    latent_space = np.concatenate(latent_space, axis=0)
    pd.DataFrame(latent_space).to_csv(filename, index=False)

# Export latent space for validation set
export_latent_space(model, val_loader)

# Visualize latent space with PCA and t-SNE
latent_space = pd.read_csv("latent_space.csv")

# PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(latent_space)

# t-SNE
tsne = TSNE(n_components=2, random_state=42)
tsne_result = tsne.fit_transform(latent_space)

# Plot
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.scatter(pca_result[:, 0], pca_result[:, 1], s=50)
plt.title("PCA of Latent Space")
plt.xlabel("PC1")
plt.ylabel("PC2")

plt.subplot(1, 2, 2)
plt.scatter(tsne_result[:, 0], tsne_result[:, 1], s=50)
plt.title("t-SNE of Latent Space")
plt.xlabel("t-SNE1")
plt.ylabel("t-SNE2")

plt.tight_layout()
plt.savefig("latent_space_visualization.png")
plt.show()

# Compute permutation feature importance
def compute_permutation_importance(model, dataloader, n_repeats=5):
    model.eval()
    baseline_latent_space = []
    for omics, x1, x2, y1_true, y2_true, covariates in dataloader:
        _, _, m = model(omics, x1, x2, covariates)
        baseline_latent_space.append(m.detach().cpu().numpy())
    baseline_latent_space = np.concatenate(baseline_latent_space, axis=0)

    feature_importances = {key: np.zeros(omics_val[key].shape[1]) for key in omics_val.keys()}

    for key in omics_val.keys():
        for i in range(omics_val[key].shape[1]):
            permuted_latent_space = []
            for omics_batch, x1, x2, y1_true, y2_true, covariates in dataloader:
                omics_copy = {k: v.clone() for k, v in omics_batch.items()}
                omics_copy[key][:, i] = omics_copy[key][torch.randperm(omics_copy[key].size(0)), i]
                _, _, m = model(omics_copy, x1, x2, covariates)
                permuted_latent_space.append(m.detach().cpu().numpy())
            permuted_latent_space = np.concatenate(permuted_latent_space, axis=0)
            feature_importances[key][i] = np.mean(np.abs(baseline_latent_space - permuted_latent_space))

    return feature_importances

# Visualize top contributing features
def plot_feature_importances(feature_importances, top_n=30):
    all_importances = []
    for key in feature_importances.keys():
        for i, importance in enumerate(feature_importances[key]):
            all_importances.append((f"{key}_{i}", importance))
    all_importances = sorted(all_importances, key=lambda x: x[1], reverse=True)[:top_n]

    feature_names, importances = zip(*all_importances)
    plt.figure(figsize=(12, 6))
    plt.barh(feature_names, importances)
    plt.xlabel("Permutation Feature Importance")
    plt.ylabel("Feature")
    plt.title("Top Contributing Features to Latent Space")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()

# Run importance analysis
feature_importances = compute_permutation_importance(model, val_loader)
plot_feature_importances(feature_importances)

# Evaluate model performance
def evaluate(model, dataloader, criterion):
    model.eval()
    total_loss = 0
    with torch.no_grad():
        for omics, x1, x2, y1_true, y2_true, covariates in dataloader:
            y1_pred, y2_pred, _ = model(omics, x1, x2, covariates)
            loss = criterion(y1_pred, y1_true) + criterion(y2_pred, y2_true)
            total_loss += loss.item()
    return total_loss / len(dataloader)

val_loss = evaluate(model, val_loader, criterion)
test_loss = evaluate(model, test_loader, criterion)

# Export model performance statistics
performance_stats = {
    "Validation Loss": val_loss,
    "Test Loss": test_loss
}
pd.DataFrame([performance_stats]).to_csv("model_performance.csv", index=False)

print("Model training and evaluation complete. Latent space and performance statistics exported.")
