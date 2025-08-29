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
