from scipy.stats import rankdata

# Using predicted probability as input for UMAP
# Keep compounds with at least 1 predicted kinase (p >= 0.5)
mask = (df_chemspace_28M_data >= 0.5).any(axis=1) # change input here
filtered_preds = df_chemspace_28M_data[mask] # change input here

binary_profiles = (filtered_preds >= 0.5).astype(int)
filtered_preds['MaxProb'] = filtered_preds.max(axis=1)
filtered_preds['HitCount'] = binary_profiles.sum(axis=1)


# Compute Gini on the rank‑normalized values
def gini(x):
    x = np.sort(np.array(x, dtype=float))
    if x.sum() == 0:
        return 0.0
    n = len(x)
    idx = np.arange(1, n + 1)
    return (2 * np.sum(idx * x) / (n * x.sum())) - (n + 1) / n

filtered_preds['Gini'] = filtered_preds.drop(columns=['HitCount', 'MaxProb']).apply(lambda r: gini(r.values), axis=1)

# UMAP embedding on raw profiles
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

order = [
    "Single-kinase hit",
    "Dual-kinase hits",
    "Triple-kinase hits",
    "Quadruple-kinase hits"
]

# Plot UMAP
plt.figure(figsize=(8, 6))
sns.scatterplot(data=df_plot, x="UMAP1", y="UMAP2", hue="Hit Count", hue_order=order, palette='tab10', s=10)
# plt.title("FDA Approved (1K) - UMAP of Predicted Active Kinase Profiles")
# plt.title("Kinase Targeted (65K) - UMAP of Predicted Active Kinase Profiles")
# plt.title("Natural Products (42K) - UMAP of Predicted Active Kinase Profiles")
# plt.title("ChemSpace (28M) - UMAP of Predicted Active Kinase Profiles")
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.tight_layout()
plt.grid(True)
# plt.show()

# plt.savefig("UMAP_FDA_PredictedProbability_Clutering_ColorByHitCount.png", dpi=300, bbox_inches="tight")
# plt.savefig("UMAP_Enamine_Kinase_PredictedProbability_Clutering_ColorByHitCount.png", dpi=300, bbox_inches="tight")
# plt.savefig("UMAP_NaturalProduct_PredictedProbability_Clutering_ColorByHitCount.png", dpi=300, bbox_inches="tight")
# plt.savefig("UMAP_ChemSpace_28M_PredictedProbability_Clutering_ColorByHitCount.png", dpi=300, bbox_inches="tight")
