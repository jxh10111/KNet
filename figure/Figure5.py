# Using cutoff of 0.5
threshold = 0.5
active_mask = (df_chemspace_28M_data >= threshold).any(axis=1) # change input here
filtered_preds = df_chemspace_28M_data[active_mask] # change input here
filtered_indices = np.where(active_mask)[0]    

# Gini coefficient
def gini(array): 
    if np.all(array == 0): return 0
    sorted_array = np.sort(array)
    n = len(array)
    cumvals = np.cumsum(sorted_array)
    gini_coeff = (n + 1 - 2 * np.sum(cumvals) / cumvals[-1]) / n
    return gini_coeff

# Compute selectivity profile metrics
ginis = np.apply_along_axis(gini, axis=1, arr=filtered_preds)
mean_pred = filtered_preds.mean(axis=1)
max_pred = filtered_preds.max(axis=1)
count_above_05 = (filtered_preds >= 0.5).sum(axis=1)

# Create feature matrix & normalize
features = np.stack([ginis, mean_pred, max_pred, count_above_05], axis=1)
features_norm = normalize(features)

# UMAP projection
embedding = umap.UMAP(n_neighbors=15, min_dist=0.1, metric='euclidean', random_state=42).fit_transform(features_norm)


df = pd.DataFrame({
    "CompoundIndex": filtered_indices,
    "Gini": ginis,
    "MeanPred": mean_pred,
    "MaxPred": max_pred,
    "CountAbove0.5": count_above_05,
    "UMAP 1": embedding[:, 0],
    "UMAP 2": embedding[:, 1]
})

hit_count_labels = {
    1: "single-kinase hit",
    2: "dual-kinase hits",
    3: "triple-kinase hits"
}

df = df[df["CountAbove0.5"].isin(hit_count_labels.keys())]
df["HitCountLabel"] = df["CountAbove0.5"].map(hit_count_labels)

custom_palette = {
    "single-kinase hit": "#1f77b4",   # blue
    "dual-kinase hits": "#2ca02c",   # green
    "triple-kinase hits": "#d62728"  # red
}

# Plot UMAP projection
plt.figure(figsize=(8, 6))
sns.scatterplot(data=df, x="UMAP 1", y="UMAP 2", hue="HitCountLabel", palette=custom_palette, s=10, linewidth=0)
# plt.title("FDA Approved (1K) - UMAP of GINI Index From Predicted Active Kinase Profiles")
# plt.title("Kinase Targeted (65K) - UMAP of GINI Index From Predicted Active Kinase Profiles")
# plt.title("Natural Products (42K) - UMAP of GINI Index From Predicted Active Kinase Profiles")
# plt.title("ChemSpace (28M) - UMAP of GINI Index From Predicted Active Kinase Profiles")
# plt.legend(title="Predicted Hits", loc='upper right')
# For chemspace and FDA-Approved
plt.legend(title="Predicted Hits", loc="lower right")
plt.grid(True)
plt.tight_layout()
# plt.show()

# plt.savefig("UMAP_FDA_GINI_PredictedProbability_Clutering_ColorByHitCount.png", dpi=300, bbox_inches="tight")
# plt.savefig("UMAP_Enamine_Kinase_GINI_PredictedProbability_Clutering_ColorByHitCount.png", dpi=300, bbox_inches="tight")
# plt.savefig("UMAP_NaturalProduct_GINI_PredictedProbability_Clutering_ColorByHitCount.png", dpi=300, bbox_inches="tight")
# plt.savefig("UMAP_ChemSpace_28M_GINI_PredictedProbability_Clutering_ColorByHitCount.png", dpi=300, bbox_inches="tight")
