import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
df_chemspace_28M = pd.read_csv('ChemSpace_28M_umap_GiniRanked_hitcount.csv', index_col=0)
df_FDA = pd.read_csv('FDA_umap_GiniRanked_hitcount.csv', index_col=0)
df_kinase = pd.read_csv('Enamine_Kinase_umap_GiniRanked_hitcount.csv', index_col=0)
df_np = pd.read_csv('NaturalProduct_umap_GiniRanked_hitcount.csv', index_col=0)

print(df_chemspace_28M.columns)
# Index(['PLK1', 'ATM', 'GSK3B', 'GSK3A', 'SRC', 'PFKFB3', 'MAPK14', 'ABL1',
#        'CDK5', 'MAPKAPK2',
#        ...
#        'PGK1', 'HitCount', 'MaxProb', 'Gini', 'UMAP1', 'UMAP2', 'Rank',
#        'Hit Count', 'RankNorm_forward', 'RankNorm_inverse'],
#       dtype='object', length=337)

target_cols = [col for col in df_chemspace_28M.columns if col not in ['HitCount', 'MaxProb', 'Gini', 'UMAP1', 'UMAP2', 'Rank', 'Hit Count', 'RankNorm_forward', 'RankNorm_inverse']]
df_chemspace_28M['MaxProb'] = df_chemspace_28M[target_cols].max(axis=1)
df_FDA['MaxProb'] = df_FDA[target_cols].max(axis=1)
df_kinase['MaxProb'] = df_kinase[target_cols].max(axis=1)
df_np['MaxProb'] = df_np[target_cols].max(axis=1)

df_chemspace_28M_subset = df_chemspace_28M[['Gini', 'MaxProb']].copy()
df_FDA_subset = df_FDA[['Gini', 'MaxProb']].copy()
df_kinase_subset = df_kinase[['Gini', 'MaxProb']].copy()
df_np_subset = df_np[['Gini', 'MaxProb']].copy()


# Normalize Gini and MaxProb
for df in [df_chemspace_28M_subset, df_FDA_subset, df_kinase_subset, df_np_subset]:
    df['Gini_norm'] = (df['Gini'] - df['Gini'].min()) / (df['Gini'].max() - df['Gini'].min())
    df['MaxProb_norm'] = (df['MaxProb'] - df['MaxProb'].min()) / (df['MaxProb'].max() - df['MaxProb'].min())

# Composite activity selectivity score
for df in [df_chemspace_28M_subset, df_FDA_subset, df_kinase_subset, df_np_subset]:
    df['SelectivityScore'] = 0.5 * df['Gini_norm'] + 0.5 * df['MaxProb_norm']

# Plot
def plot_library(df, library_name, filename=None):
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(
        df['Gini_norm'],
        df['MaxProb_norm'],
        c=df['SelectivityScore'],
        cmap='plasma_r', # viridis, viridis_r
        vmin=0,  # ensure colorbar starts at 0
        vmax=1,  # ensure colorbar ends at 1
        s=15,
        edgecolor='none'
    )
    plt.colorbar(scatter, label='Composite Activity Selectivity Score')
    plt.xlabel('Gini Coefficient')
    # plt.ylabel('Max Predicted Probability')
    plt.ylabel('Max Predicted Activity Score')
    plt.title(f'{library_name}')
    plt.grid(True)
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')

    # plt.show()

plot_library(df_chemspace_28M_subset, "ChemSpace (28M)", "chemspace_28M_landscape_red_reverse.png")
plot_library(df_FDA_subset, "FDA Approved (1K)", "FDA_landscape_red_reverse.png")
plot_library(df_kinase_subset, "Kinase Targeted (65K)", "kinase_landscape_red_reverse.png")
plot_library(df_np_subset, "Natural Products (42K)", "natural_products_landscape_red_everse.png")


# # save dfs with GINI and MaxProb
# df_chemspace_28M_subset.to_csv("ChemSpace_28M_GINI_MaxProb.csv", index=True)
# df_FDA_subset.to_csv("FDA_GINI_MaxProb.csv", index=True)
# df_kinase_subset.to_csv("Enamine_Kinase_GINI_MaxProb.csv", index=True)
# df_np_subset.to_csv("NaturalProduct_GINI_MaxProb.csv", index=True)















