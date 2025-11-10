import pandas as pd
import numpy as np
import os
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.impute import SimpleImputer
from sklearn.impute import KNNImputer
from tensorflow.keras.models import Model
from tensorflow.keras.layers import (
    Input, Dense, Concatenate, MultiHeadAttention, LayerNormalization,
    Add, Lambda, Flatten, Dropout, BatchNormalization
)
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.optimizers import Adam
import shap
import tensorflow.keras.backend as K
import pymc as pm
import arviz as az
import seaborn as sns
import matplotlib.pyplot as plt
import shap

# Load data
df = pd.read_csv("combined_df_new.csv")
imputer = KNNImputer(n_neighbors=5)
df.iloc[:, 10:11] = imputer.fit_transform(df.iloc[:, 10:11])

# Define column ranges
omics_ranges = {
    "Metabolomics": (12, 176),
    "KEGG_KO": (176, 250),
    "KEGG_Module": (250,254),
    "KEGG_pathway" : (254,255),
}

# Define covariates, x, y
covariates = df.iloc[:, 2:8]
x = df.iloc[:, 8:10]
y = df.iloc[:, 10:12]

# Split data
train_idx, test_idx = train_test_split(df.index, test_size=0.3, random_state=42)
val_idx, test_idx = train_test_split(test_idx, test_size=0.5, random_state=42)

# Normalize covariates, x, y using training set statistics
scalers = {}
scaled_data = {}

for name, data in zip(['covariates', 'x', 'y'], [covariates, x, y]):
    scaler = StandardScaler()
    scaler.fit(data.loc[train_idx])  # Fit only on training data
    scaled_data[name] = scaler.transform(data)  # Transform the full dataset
    scalers[name] = scaler

omics_data = {}
omics_inputs = {}
encoded_outputs = []
input_feature_names = []

# Process non-KEGG omics blocks
for name, (start, end) in omics_ranges.items():
    if name.startswith("KEGG"):
        continue
    block = df.iloc[:, start:end]
    feature_names = block.columns.tolist()
    input_feature_names.extend(feature_names)

    scaler = StandardScaler()
    block_scaled = scaler.fit_transform(block.loc[train_idx])
    omics_data[name] = StandardScaler().fit_transform(block)

    input_layer = Input(shape=(block_scaled.shape[1],), name=f"{name}_input")
    dense = Dense(64, activation='relu', name=f"{name}_dense")(input_layer)
    omics_inputs[name] = input_layer
    encoded_outputs.append(dense)
omics_data = {}
omics_inputs = {}
encoded_outputs = []
input_feature_names = []

# 🧬 Process non-KEGG omics blocks
for name, (start, end) in omics_ranges.items():
    if name.startswith("KEGG"):
        continue
    block = df.iloc[:, start:end]
    feature_names = block.columns.tolist()
    input_feature_names.extend(feature_names)

    scaler = StandardScaler()
    scaler.fit(block.loc[train_idx])  # Fit only on training data
    block_scaled = scaler.transform(block.loc[train_idx])  # For input shape
    omics_data[name] = scaler.transform(block)  # Transform full dataset

    input_layer = Input(shape=(block_scaled.shape[1],), name=f"{name}_input")
    dense = Dense(64, activation='relu', name=f"{name}_dense")(input_layer)
    omics_inputs[name] = input_layer
    encoded_outputs.append(dense)

# 🧬 Combined KEGG block
kegg_block = pd.concat([
    df.iloc[:, omics_ranges['KEGG_KO'][0]:omics_ranges['KEGG_KO'][1]],
    df.iloc[:, omics_ranges['KEGG_Module'][0]:omics_ranges['KEGG_Module'][1]],
    df.iloc[:, omics_ranges['KEGG_pathway'][0]:omics_ranges['KEGG_pathway'][1]]
], axis=1)

# Impute and scale
from sklearn.impute import SimpleImputer
kegg_imputer = SimpleImputer(strategy='mean')
kegg_block_imputed = kegg_imputer.fit_transform(kegg_block.loc[train_idx])
kegg_scaler = StandardScaler()
kegg_block_scaled = kegg_scaler.fit_transform(kegg_block_imputed)
omics_data['KEGG'] = kegg_scaler.transform(kegg_imputer.transform(kegg_block))

# Input and encoding
kegg_input = Input(shape=(kegg_block_scaled.shape[1],), name="KEGG_combined_input")
kegg_dense = Dense(64, activation='relu', name="KEGG_combined_dense")(kegg_input)
kegg_dense = BatchNormalization()(kegg_dense)
kegg_dense = Dropout(0.2)(kegg_dense)

omics_inputs['KEGG'] = kegg_input
encoded_outputs.append(kegg_dense)

# === Stack encoded outputs for cross-modal attention ===
stacked = Lambda(lambda x: tf.stack(x, axis=1))(encoded_outputs)


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

# === Build and compile the model ===
all_omics_inputs = []
for v in omics_inputs.values():
    all_omics_inputs.extend(v if isinstance(v, list) else [v])

model = Model(inputs=all_omics_inputs + [x_input, y_input, cov_input], outputs=output)
model.compile(optimizer=Adam(1e-3), loss='mse')

# === Prepare input data ===
def get_inputs(indices):
    return (
        [omics_data[name][indices] for name in omics_ranges if not name.startswith("KEGG")] +
        [omics_data['KEGG'][indices]] +
        [x.loc[indices], y.loc[indices], covariates.loc[indices]]
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
    epochs=100, batch_size=64,
    callbacks=[EarlyStopping(patience=10, restore_best_weights=True)]
)

# === Extract latent space representations ===
latent_model = Model(inputs=all_omics_inputs, outputs=latent)
latent_full = latent_model.predict(
    [omics_data[name] for name in omics_ranges if not name.startswith("KEGG")] +
    [omics_data['KEGG']]
)
latent_train = latent_full[train_idx]
latent_test = latent_full[test_idx]
model.save("outputs4/latent_model.keras")

# === Use training data ===
x_array = x.loc[train_idx].values
y_array = y.loc[train_idx].values

results = []

for i in range(x_array.shape[1]):
    for j in range(y_array.shape[1]):
        x_var = x_array[:, i]
        y_var = y_array[:, j]

        with pm.Model() as model:
            beta_latent = pm.Normal("beta_latent", mu=0, sigma=1, shape=latent_train.shape[1])
            beta_x = pm.Normal("beta_x", mu=0, sigma=1)
            intercept = pm.Normal("intercept", mu=0, sigma=1)
            sigma = pm.HalfNormal("sigma", sigma=1)

            mu = intercept + pm.math.dot(latent_train, beta_latent) + beta_x * x_var
            y_obs = pm.StudentT("y_obs", mu=mu, sigma=sigma, nu=4, observed=y_var)

            trace = pm.sample(
                1000,
                tune=500,
                target_accept=0.9,
                progressbar=False,
                return_inferencedata=True,
                idata_kwargs={"log_likelihood": True}  # ✅ Ensures WAIC/LOO compatibility
            )

        summary = az.summary(trace, var_names=["beta_x", "beta_latent"], hdi_prob=0.95)
        direct = summary.loc["beta_x", "mean"]
        direct_hdi = pd.Series(summary.loc["beta_x", ["hdi_2.5%", "hdi_97.5%"]])
        indirect = summary.filter(like="beta_latent", axis=0)["mean"].sum()
        indirect_hdi = pd.Series(summary.filter(like="beta_latent", axis=0)[["hdi_2.5%", "hdi_97.5%"]].sum(axis=0))
        total = direct + indirect
        total_hdi = direct_hdi + indirect_hdi

        results.append({
            "x_var": f"x{i+1}",
            "y_var": f"y{j+1}",
            "direct_effect": direct,
            "direct_hdi_2.5%": direct_hdi.iloc[0],
            "direct_hdi_97.5%": direct_hdi.iloc[1],
            "indirect_effect": indirect,
            "indirect_hdi_2.5%": indirect_hdi.iloc[0],
            "indirect_hdi_97.5%": indirect_hdi.iloc[1],
            "total_effect": total,
            "total_hdi_2.5%": total_hdi.iloc[0],
            "total_hdi_97.5%": total_hdi.iloc[1]
        })



# === 1. Bayesian diagnostics ===
os.makedirs("outputs4", exist_ok=True)
az_summary = az.summary(trace, var_names=["beta_x", "beta_latent"], hdi_prob=0.95)
az_summary.to_csv("outputs4/bayesian_summary.csv")

idata = trace
waic = az.waic(idata)
loo = az.loo(idata)
with open("outputs4/model_diagnostics.txt", "w") as f:
    f.write(f"WAIC:\n{waic}\n\n")
    f.write(f"LOO:\n{loo}\n")

# === 2. Correlation between latent space and y ===
latent_df = pd.DataFrame(latent_full, columns=[f"latent_{i}" for i in range(latent_full.shape[1])])
y_df = pd.DataFrame(y.values, columns=[f"y{i+1}" for i in range(y.shape[1])])

correlation_matrix = latent_df.corrwith(y_df.iloc[:, 0]).to_frame(name="y1")
if y.shape[1] > 1:
    correlation_matrix["y2"] = latent_df.corrwith(y_df.iloc[:, 1])

correlation_matrix.to_csv("outputs4/latent_y_correlation.csv")

# === 3. Heatmap ===
plt.figure(figsize=(10, 6))
sns.heatmap(correlation_matrix, annot=True, cmap="coolwarm", center=0)
plt.title("Correlation between Latent Dimensions and y")
plt.tight_layout()
plt.savefig("outputs4/latent_y_correlation_heatmap.png")
plt.close()

# === 4. Save other outputs ===
np.savetxt("outputs4/latent_space.csv", latent_full, delimiter=",")

# Combine all omics data into a single array for SHAP
omics_combined = np.concatenate(
    [omics_data[name] for name in omics_ranges if not name.startswith("KEGG")] +
    [omics_data["KEGG"]],
    axis=1
)

# Combine feature names
feature_names = []
for name, (start, end) in omics_ranges.items():
    if name.startswith("KEGG"):
        continue
    feature_names.extend(df.columns[start:end])
feature_names.extend(kegg_block.columns)

# Run SHAP
explainer = shap.DeepExplainer(latent_model, omics_combined)
shap_values = explainer.shap_values(omics_combined)

# Save SHAP summary plots
os.makedirs("outputs4", exist_ok=True)
for i in range(len(shap_values)):
    shap.summary_plot(shap_values[i], omics_combined, feature_names=feature_names, show=False)
    plt.title(f"Latent Dimension {i}")
    plt.savefig(f"outputs4/shap_summary_plot_latent_dim_{i}.png")
    plt.close()


results_df = pd.DataFrame(results)
results_df.to_csv("outputs4/bayesian_sem_effects.csv", index=False)
print("✅ All Bayesian SEM outputs saved in the 'outputs4' directory.")
print("✅ SHAP analysis completed and summary plots saved.")


