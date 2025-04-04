import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from rdkit.Chem import Draw
from sklearn.cluster import KMeans
from sklearn.mixture import BayesianGaussianMixture
from sklearn.preprocessing import StandardScaler
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

# -------------------------------
# Step 1: Load Three CSVs and Combine Them
# -------------------------------
files = {
    'Natural Products': 'data/ChEMBL31_KKBQ122_dc_keras_model_328targets_cutoff7_KNet_Pred_Natural_Products_MEGx_NATx_41820_compounds.csv',
    'FDA-approved Drugs': 'data/ChEMBL31_KKBQ122_dc_keras_model_328targets_cutoff7_KNet_Pred_Enamine_FDA_approved_Drugs_1122cmpds.csv',
    'Kinase-targeted Library': 'data/ChEMBL31_KKBQ122_dc_keras_model_328targets_cutoff7_KNet_Pred_Enamine_Kinase_Targeted_Library_64880_cpds.csv'
}

dfs = []
for lib, file in files.items():
    df = pd.read_csv(file)
    df['library'] = lib  # Add library label
    dfs.append(df)
# Combine data from all libraries
data = pd.concat(dfs, ignore_index=True)

# Identify gene (kinase) columns – exclude 'Canonical_Smiles' and 'library'
gene_cols = [col for col in data.columns if col.lower() not in ['canonical_smiles', 'library']]

# Convert gene columns to numeric (coerce errors to NaN) and drop rows with missing values
data[gene_cols] = data[gene_cols].apply(pd.to_numeric, errors='coerce')
data.dropna(subset=gene_cols, inplace=True)

# -------------------------------
# Step 2: Compute the GINI Index & Library-Specific Selectivity
# -------------------------------
def gini_coefficient(x):
    """Compute the GINI coefficient for a numpy array of non-negative values."""
    arr = np.array(x, dtype=float)
    if np.sum(arr) == 0:
        return 0
    sorted_arr = np.sort(arr)
    n = len(arr)
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * sorted_arr) - (n + 1) * np.sum(sorted_arr)) / (n * np.sum(sorted_arr))

# Compute the GINI index for each compound (row-wise across gene columns)
data['gini'] = data[gene_cols].apply(lambda row: gini_coefficient(row), axis=1)

# Compute the median GINI index for each library
library_threshold = data.groupby('library')['gini'].median().to_dict()

print("Library GINI thresholds:")
for lib, thresh in library_threshold.items():
    print(f"{lib}: {thresh}")

# Assign selectivity per compound using the library-specific threshold
data['selectivity'] = data.apply(
    lambda row: 'Selective' if row['gini'] >= library_threshold[row['library']] else 'Unselective',
    axis=1
)

# (Optional) Defragment the DataFrame to improve performance warnings
data = data.copy()

# -------------------------------
# Step 3: Generate Chemical Features and Perform Clustering
# -------------------------------
# Create a Morgan fingerprint generator using the new RDKit API
mg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

def smiles_to_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # Generate fingerprint using the new generator API
    fp = mg.GetFingerprint(mol)
    # Use the fingerprint's GetNumBits() to obtain its size
    arr = np.zeros((fp.GetNumBits(),), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

# Use "Canonical_Smiles" column for SMILES strings
data['features'] = data['Canonical_Smiles'].apply(smiles_to_fp)
data = data[data['features'].notnull()].copy()

# Combine features from all libraries into one matrix for clustering
features_matrix = np.vstack(data['features'].values)

# Run k-means clustering on the combined data (using n_init='auto' for compatibility with latest scikit-learn)
n_clusters = 5  # Adjust as needed
kmeans = KMeans(n_clusters=n_clusters, n_init='auto', random_state=42)
data['cluster'] = kmeans.fit_predict(features_matrix)

# Check that the total number of rows equals the sum of the cluster counts
total_rows = len(data)
total_cluster_count = data['cluster'].value_counts().sum()
print("Total rows in data:", total_rows)
print("Sum of cluster counts:", total_cluster_count)
assert total_rows == total_cluster_count, "Mismatch between total rows and cluster counts!"

print("Overall cluster counts:")
print(data['cluster'].value_counts())

print("\nCluster counts by library:")
print(data.groupby(['cluster', 'library']).size())

print("\nCluster counts by library and selectivity:")
print(data.groupby(['cluster', 'library', 'selectivity']).size())

for clust in [2, 3]:
    for lib in data['library'].unique():
         subset = data[(data['cluster'] == clust) & (data['library'] == lib)]
         if len(subset) > 0:
             print(f"Cluster {clust}, Library {lib}: count = {len(subset)}")
             print(f"  Min GINI: {subset['gini'].min():.3f}, Max GINI: {subset['gini'].max():.3f}, Median GINI: {subset['gini'].median():.3f}")
         else:
             print(f"Cluster {clust}, Library {lib}: no compounds")

# -------------------------------
# Step 4: Create the Barplot (with Representative Chemotype Structures)
# -------------------------------
# Aggregate data to compute the mean GINI index for each combination of cluster, library, and selectivity
agg = data.groupby(['cluster', 'library', 'selectivity'])['gini'].mean().reset_index()
agg['lib_sel'] = agg['library'] + '_' + agg['selectivity']

plt.figure(figsize=(12, 6))
ax = sns.barplot(x='cluster', y='gini', hue='lib_sel', data=agg, palette='viridis')
plt.xlabel('Chemotype Cluster')
plt.ylabel('Average Selectivity (GINI) Index')
plt.title('Selective vs. Unselective Compounds by Chemotype and Library')
plt.legend(title='Library_Selectivity', bbox_to_anchor=(1.05, 1), loc='upper left')

# Compute cluster centers (average fingerprint for each cluster)
cluster_centers = data.groupby('cluster')['features'].apply(lambda x: np.mean(np.vstack(x), axis=0)).reset_index()

# For each cluster, find the compound whose fingerprint is closest to the cluster center.
rep_mols = {}
for cl in sorted(data['cluster'].unique()):
    center = cluster_centers.loc[cluster_centers['cluster'] == cl, 'features'].iloc[0]
    cluster_data = data[data['cluster'] == cl]
    distances = cluster_data['features'].apply(lambda x: np.linalg.norm(x - center))
    rep_index = distances.idxmin()  # index of the compound closest to the center
    rep_mol = Chem.MolFromSmiles(data.loc[rep_index, 'Canonical_Smiles'])
    rep_mols[cl] = rep_mol

# Annotate the barplot with representative molecule images for each cluster
for cl in sorted(rep_mols.keys()):
    mol = rep_mols[cl]
    img = Draw.MolToImage(mol, size=(150, 150))
    imagebox = OffsetImage(img, zoom=0.5)
    ab = AnnotationBbox(imagebox, (cl, ax.get_ylim()[0]), frameon=False, box_alignment=(0.5, -0.1))
    ax.add_artist(ab)

plt.tight_layout()
plt.show()

# -------------------------------
# Step 5: Create the Scatterplot (Cluster Similarity vs. Selectivity using BGMM)
# -------------------------------
# Compute cluster centers (average fingerprint for each cluster) again
cluster_centers = data.groupby('cluster')['features'].apply(lambda x: np.mean(np.vstack(x), axis=0)).reset_index()
centers_matrix = np.vstack(cluster_centers['features'])

# Fit a Bayesian Gaussian Mixture Model on the cluster centers.
# Here we use 2 components as an example.
bgmm = BayesianGaussianMixture(n_components=2, random_state=42)
bgmm.fit(centers_matrix)

# Compute the posterior probability for component 0 as the similarity measure.
similarity = bgmm.predict_proba(centers_matrix)[:, 0]
cluster_centers['similarity'] = similarity

# Merge the mean GINI index for each cluster (aggregated across libraries)
mean_gini = data.groupby('cluster')['gini'].mean().reset_index()
cluster_centers = cluster_centers.merge(mean_gini, on='cluster')

plt.figure(figsize=(8, 6))
sns.scatterplot(x='similarity', y='gini', hue='cluster', data=cluster_centers, palette='viridis', s=100)
plt.xlabel('Similarity (BGMM probability for component 0)')
plt.ylabel('Mean Selectivity (GINI) Index')
plt.title('Cluster Similarity vs. Selectivity (via BGMM)')
plt.tight_layout()
plt.show()
