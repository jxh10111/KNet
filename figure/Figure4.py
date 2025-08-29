import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import normalize
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import hdbscan
from scipy.stats import rankdata


# Load data
df_chemspace_28M = pd.read_csv('ChEMBL34_KKBQ224_dc_keras_model_121924_328targets_cutoff7_KNet_Pred_ChemSpace28M_actives_only.csv')
df_FDA = pd.read_csv('ChEMBL34_KKBQ224_dc_keras_model_328targets_cutoff7_KNet_Pred_Enamine_FDA_approved_Drugs_1122cmpds.csv')
df_kinase = pd.read_csv('ChEMBL34_KKBQ224_dc_keras_model_328targets_cutoff7_KNet_Pred_Enamine_Kinase_Targeted_Library_64880_cpds.csv')
df_np = pd.read_csv('ChEMBL34_KKBQ224_dc_keras_model_328targets_cutoff7_KNet_Pred_Natural_Products_MEGx_NATx_41820_compounds.csv')

df_chemspace_28M.rename(columns={'smiles': 'Canonical_Smiles'}, inplace=True)

# Keep compounds with any prediction >= 0.5
df_FDA_data = df_FDA.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')
df_kinase_data = df_kinase.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')
df_np_data = df_np.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')
df_chemspace_28M_data = df_chemspace_28M.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')


# Using predicted probability as input for UMAP
# Keep compounds with at least 1 predicted kinase (p >= 0.5)
mask = (df_chemspace_28M_data >= 0.5).any(axis=1) # change input here
filtered_preds = df_chemspace_28M_data[mask] # change input here

binary_profiles = (filtered_preds >= 0.5).astype(int)
filtered_preds['MaxProb'] = filtered_preds.max(axis=1)
filtered_preds['HitCount'] = binary_profiles.sum(axis=1)


# Compute Gini
def gini(x):
    x = np.sort(np.array(x, dtype=float))
    if x.sum() == 0:
        return 0.0
    n = len(x)
    idx = np.arange(1, n + 1)
    return (2 * np.sum(idx * x) / (n * x.sum())) - (n + 1) / n

filtered_preds['Gini'] = filtered_preds.drop(columns=['HitCount', 'MaxProb']).apply(lambda r: gini(r.values), axis=1)

# UMAP embedding
umap_model = umap.UMAP(n_neighbors=15, min_dist=0.01, metric='euclidean', random_state=42)
coords = umap_model.fit_transform(filtered_preds.drop(columns=['HitCount', 'MaxProb', 'Gini']).values)
filtered_preds['UMAP1'], filtered_preds['UMAP2'] = coords[:,0], coords[:,1]

# Add Rank for Gini
filtered_preds['Rank'] = filtered_preds['Gini'] \
    .rank(method='dense', ascending=False) \
    .astype(int)


labels = {
    1: "Single-kinase hit",
    2: "Dual-kinase hits",
    3: "Triple-kinase hits",
    4: "Quadruple-kinase hits"
}
df_plot = filtered_preds[filtered_preds['HitCount'].isin(labels)].copy()
df_plot['Hit Count'] = df_plot['HitCount'].map(labels)

# Compute min and max
gmin, gmax = df_plot['Gini'].min(), df_plot['Gini'].max()

# Create a new normalized column
df_plot['Gini_norm'] = (df_plot['Gini'] - gmin) / (gmax - gmin)
print(df_plot['Gini_norm'].min(), df_plot['Gini_norm'].max())
# 0.0 1.0

# Plot UMAP colored by Gini
plt.figure(figsize=(8, 6))
sc = plt.scatter(
    df_plot['UMAP1'],
    df_plot['UMAP2'],
    c=df_plot['Gini_norm'],
    cmap='plasma', # viridis_r or plasma_r
    s=12,
    edgecolor='none'
)
plt.colorbar(sc, label='Gini Index (Normalized)')
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
# plt.title("FDA Approved (1K) - UMAP of Predicted Active Kinase Profiles")
# plt.title("Kinase Targeted (65K) - UMAP of Predicted Active Kinase Profiles")
# plt.title("Natural Products (42K) - UMAP of Predicted Active Kinase Profiles")
plt.title('ChemSpace (28M) - UMAP of Predicted Active Kinase Profiles')
plt.grid(True)
plt.tight_layout()
plt.show()

# plt.savefig("UMAP_FDA_PredictedProbability_Clutering_ColorByGINI.png", dpi=300, bbox_inches="tight")
# plt.savefig("UMAP_Enamine_Kinase_PredictedProbability_Clutering_ColorByGINI.png", dpi=300, bbox_inches="tight")
# plt.savefig("UMAP_NaturalProduct_PredictedProbability_Clutering_ColorByGINI.png", dpi=300, bbox_inches="tight")
# plt.savefig("UMAP_ChemSpace_28M_PredictedProbability_Clutering_ColorByGINI.png", dpi=300, bbox_inches="tight")


# Save dfs with GINI and Rank and HitCount information
# df_plot.to_csv("FDA_umap_GiniRanked_hitcount.csv", index=True)
# df_plot.to_csv("Enamine_Kinase_umap_GiniRanked_hitcount.csv", index=True)
# df_plot.to_csv("NaturalProduct_umap_GiniRanked_hitcount.csv", index=True)
# df_plot.to_csv("ChemSpace_28M_umap_GiniRanked_hitcount.csv", index=True)
