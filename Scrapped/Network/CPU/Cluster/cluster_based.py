import pandas as pd
import numpy as np
import os
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer, KNNImputer
from tensorflow.keras.models import Model
from tensorflow.keras.layers import (
    Input, Dense, Concatenate, MultiHeadAttention, LayerNormalization,
    Add, Lambda, Flatten, Dropout, BatchNormalization
)
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.optimizers import Adam
import tensorflow.keras.backend as K
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
import seaborn as sns
import shap
import pingouin as pg
from scipy.stats import pearsonr
import tensorflow as tf

# Load data
cluster1 = pd.read_csv("cluster1_features.csv", index_col=0)
cluster2 = pd.read_csv("cluster2_features.csv", index_col=0)
meta = pd.read_csv("meta_data.csv", index_col=0)

# Extract covariates, x, y
covariates = meta.iloc[:, 0:6]
x = meta.iloc[:, 6:8]
y = meta.iloc[:, 8:10]

# Impute missing values in y using KNN
from sklearn.impute import KNNImputer
knn_imputer = KNNImputer(n_neighbors=5)
y.iloc[:, :] = knn_imputer.fit_transform(y)

train_idx, test_idx = train_test_split(meta.index, test_size=0.3, random_state=42)
val_idx, test_idx = train_test_split(test_idx, test_size=0.5, random_state=42)

# Normalize covariates, x, y using training set statistics
scalers = {}
scaled_data = {}
for name, data in zip(['covariates', 'x', 'y'], [covariates, x, y]):
    scaler = StandardScaler()
    scaler.fit(data.loc[train_idx])
    scaled_data[name] = scaler.transform(data)
    scalers[name] = scaler


omics_data = {}
omics_inputs = {}
encoded_outputs = []

for i, cluster_df in enumerate([cluster1, cluster2], start=1):
    cluster_name = f"cluster_{i}"
    valid_idx = cluster_df.index.intersection(train_idx)
    if len(valid_idx) == 0:
        print(f"⚠️ Skipping {cluster_name}: no training samples.")
        continue

    scaler = StandardScaler()
    cluster_scaled = scaler.fit_transform(cluster_df.loc[valid_idx])
    omics_data[cluster_name] = scaler.transform(cluster_df)

    input_layer = Input(shape=(cluster_scaled.shape[1],), name=f"{cluster_name}_input")
    dense = Dense(64, activation='relu', name=f"{cluster_name}_dense")(input_layer)
    omics_inputs[cluster_name] = input_layer
    encoded_outputs.append(dense)

# === Stack encoded outputs for cross-modal attention ===
stacked = Lambda(lambda x: K.stack(x, axis=1))(encoded_outputs)  # Shape: (batch, num_blocks, features)

# === Cross-modal attention block ===
attention = MultiHeadAttention(num_heads=4, key_dim=16)(stacked, stacked)
attention = Add()([stacked, attention])  # Residual connection
attention = LayerNormalization()(attention)

# === Flatten and project to latent space ===
flattened = Flatten()(attention)
latent = Dense(64, activation='relu', name="attention_latent_space")(flattened)

# === Mediation model inputs ===
x_input = Input(shape=(x.shape[1],), name="x_input")
y_input = Input(shape=(y.shape[1],), name="y_input")
cov_input = Input(shape=(covariates.shape[1],), name="covariates_input")

# === Mediation layers ===
mediator = Concatenate()([latent, x_input, y_input, cov_input])
mediator = Dense(64, activation='relu')(mediator)
mediator = Dense(32, activation='relu')(mediator)

# === Output layer ===
output = Dense(latent.shape[1], activation='linear', name="reconstructed_latent")(mediator)

# === Prepare input data ===
def get_inputs(indices):
    return (
        [omics_data[name][indices] for name in omics_data] +
        [x.loc[indices].values, y.loc[indices].values, covariates.loc[indices].values]
    )

train_inputs = get_inputs(train_idx)
val_inputs = get_inputs(val_idx)
test_inputs = get_inputs(test_idx)

# === Cache predictions before training ===
train_preds = model.predict(train_inputs, batch_size=64)
val_preds = model.predict(val_inputs, batch_size=64)

# === Train the model ===
model.fit(
    train_inputs, train_preds,
    validation_data=(val_inputs, val_preds),
    epochs=100,
    batch_size=64,
    callbacks=[EarlyStopping(patience=10, restore_best_weights=True)]
)

# === Extract latent space representations ===
latent_model = Model(inputs=all_omics_inputs, outputs=latent)

latent_full = latent_model.predict(
    [omics_data[name] for name in omics_data]
)

latent_train = latent_full[train_idx]
latent_test = latent_full[test_idx]


# === Use training data for SEM ===
x_array = x.loc[train_idx].values
y_array = y.loc[train_idx].values
results = []

for i in range(x_array.shape[1]):
    for j in range(y_array.shape[1]):
        x_var = x_array[:, i]
        y_var = y_array[:, j]
        with pm.Model() as sem_model:
            beta_latent = pm.Normal("beta_latent", mu=0, sigma=1, shape=latent_train.shape[1])
            beta_x = pm.Normal("beta_x", mu=0, sigma=1)
            intercept = pm.Normal("intercept", mu=0, sigma=1)
            sigma = pm.HalfNormal("sigma", sigma=1)
            mu = intercept + pm.math.dot(latent_train, beta_latent) + beta_x * x_var
            y_obs = pm.Normal("y_obs", mu=mu, sigma=sigma, observed=y_var)
            trace = pm.sample(1000, tune=500, target_accept=0.9, progressbar=False)
        
        summary = az.summary(trace, var_names=["beta_x", "beta_latent"], hdi_prob=0.95)
        direct = summary.loc["beta_x", "mean"]
        direct_hdi = summary.loc["beta_x", ["hdi_2.5%", "hdi_97.5%"]].values
        indirect = summary.filter(like="beta_latent", axis=0)["mean"].sum()
        indirect_hdi = summary.filter(like="beta_latent", axis=0)[["hdi_2.5%", "hdi_97.5%"]].sum(axis=0)
        total = direct + indirect
        total_hdi = direct_hdi + indirect_hdi

        results.append({
            "x_var": f"x{i+1}",
            "y_var": f"y{j+1}",
            "direct_effect": direct,
            "direct_hdi_2.5%": direct_hdi[0],
            "direct_hdi_97.5%": direct_hdi[1],
            "indirect_effect": indirect,
            "indirect_hdi_2.5%": indirect_hdi[0],
            "indirect_hdi_97.5%": indirect_hdi[1],
            "total_effect": total,
            "total_hdi_2.5%": total_hdi[0],
            "total_hdi_97.5%": total_hdi[1]
        })

# === Save Bayesian diagnostics ===
os.makedirs("outputs_cluster", exist_ok=True)
az_summary = az.summary(trace, var_names=["beta_x", "beta_latent"], hdi_prob=0.95)
az_summary.to_csv("outputs_cluster/bayesian_summary.csv")
idata = az.convert_to_inference_data(trace)
waic = az.waic(idata)
loo = az.loo(idata)
with open("outputs_cluster/model_diagnostics.txt", "w") as f:
    f.write(f"WAIC:\n{waic}\n\n")
    f.write(f"LOO:\n{loo}\n")

# === Correlation between latent space and y ===
latent_df = pd.DataFrame(latent_full, columns=[f"latent_{i}" for i in range(latent_full.shape[1])])
y_df = pd.DataFrame(y.values, columns=[f"y{i+1}" for i in range(y.shape[1])])
correlation_matrix = pd.DataFrame(index=latent_df.columns, columns=y_df.columns)
pval_matrix = pd.DataFrame(index=latent_df.columns, columns=y_df.columns)

for y_col in y_df.columns:
    for latent_col in latent_df.columns:
        r, p = pearsonr(latent_df[latent_col], y_df[y_col])
        correlation_matrix.loc[latent_col, y_col] = r
        pval_matrix.loc[latent_col, y_col] = p

# Annotate with asterisks for significance
annot_matrix = correlation_matrix.copy().astype(str)
for row in annot_matrix.index:
    for col in annot_matrix.columns:
        p = pval_matrix.loc[row, col]
        stars = '*' if p < 0.05 else ''
        annot_matrix.loc[row, col] = f"{float(correlation_matrix.loc[row, col]):.2f}{stars}"

correlation_matrix.to_csv("outputs_cluster/latent_y_correlation.csv")
pval_matrix.to_csv("outputs_cluster/latent_y_pvalues.csv")

# === Heatmap with asterisks ===
plt.figure(figsize=(10, 6))
sns.heatmap(correlation_matrix.astype(float), annot=annot_matrix, fmt='', cmap="coolwarm", center=0)
plt.title("Correlation between Latent Dimensions and y (with significance)")
plt.tight_layout()
plt.savefig("outputs_cluster/latent_y_correlation_heatmap.png")
plt.close()

# === Partial correlation adjusted for covariates ===
partial_corrs = []
for y_col in y_df.columns:
    for latent_col in latent_df.columns:
        df_partial = pd.concat([latent_df[latent_col], y_df[y_col], covariates.reset_index(drop=True)], axis=1)
        df_partial.columns = ['latent', 'y'] + [f'cov{i}' for i in range(covariates.shape[1])]
        pcorr = pg.partial_corr(data=df_partial, x='latent', y='y', covar=df_partial.columns[2:], method='pearson')
        partial_corrs.append({
            'latent': latent_col,
            'y': y_col,
            'r': pcorr['r'].values[0],
            'p-val': pcorr['p-val'].values[0]
        })

partial_corr_df = pd.DataFrame(partial_corrs)
partial_corr_df.to_csv("outputs_cluster/partial_correlations.csv", index=False)

# === Save latent space and SEM results ===
np.savetxt("outputs_cluster/latent_space.csv", latent_full, delimiter=",")
results_df = pd.DataFrame(results)
results_df.to_csv("outputs_cluster/bayesian_sem_effects.csv", index=False)

# === SHAP Analysis ===

# Combine all omics data into a single array for SHAP
omics_combined = np.concatenate(
    [omics_data[name] for name in omics_data],
    axis=1
)

# Combine feature names
feature_names = []
for name in omics_data:
    feature_names.extend([f"{name}_{i}" for i in range(omics_data[name].shape[1])])

# Run SHAP
explainer = shap.DeepExplainer(latent_model, omics_combined)
shap_values = explainer.shap_values(omics_combined)

# Save SHAP summary plots
os.makedirs("outputs_cluster/shap_plots", exist_ok=True)
for i in range(len(shap_values)):
    shap.summary_plot(shap_values[i], omics_combined, feature_names=feature_names, show=False)
    plt.title(f"Latent Dimension {i}")
    plt.savefig(f"outputs_cluster/shap_plots/shap_summary_latent_dim_{i}.png")
    plt.close()

print("✅ All SHAP summary plots saved in 'outputs_cluster/shap_plots'")



