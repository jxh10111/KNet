import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem import rdMolDescriptors
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
from sklearn.mixture import BayesianGaussianMixture


FDA_approved_lib = pd.read_csv('ChEMBL34_KKBQ224_dc_keras_model_328targets_cutoff7_KNet_Pred_Enamine_FDA_approved_Drugs_1122cmpds.csv')
Kinase_targeted_lib = pd.read_csv('ChEMBL34_KKBQ224_dc_keras_model_328targets_cutoff7_KNet_Pred_Enamine_Kinase_Targeted_Library_64880_cpds.csv')
Natural_products_lib = pd.read_csv('ChEMBL34_KKBQ224_dc_keras_model_328targets_cutoff7_KNet_Pred_Natural_Products_MEGx_NATx_41820_compounds.csv')
ChemSpace_86K = pd.read_csv('ChEMBL34_KKBQ224_dc_keras_model_121924_328targets_cutoff7_KNet_Pred_ChemSpace0003_86K.csv', sep='\t')
ChemSpace_86K.rename(columns={'smiles': 'Canonical_Smiles'}, inplace=True)

FDA_approved_lib['Library'] = 'FDA Approved (1K)' 
Kinase_targeted_lib['Library'] = 'Kinase Targeted (65K)'
Natural_products_lib['Library'] = 'Natural Products (42K)'
ChemSpace_86K['Library'] = 'ChemSpace (86K)'

# Concatenate the dataframes into one
combined_df = pd.concat([
    FDA_approved_lib[['Canonical_Smiles', 'Library']],
    Kinase_targeted_lib[['Canonical_Smiles', 'Library']],
    Natural_products_lib[['Canonical_Smiles', 'Library']],
    ChemSpace_86K[['Canonical_Smiles', 'Library']]
], ignore_index=True)

# Save combined file
combined_df.to_csv('ChEMBL34_KKBQ224_dc_keras_model_4Libs_combined_df.csv', index=False)

# Group by the 'Canonical_Smiles'
overlap_df = combined_df.groupby('Canonical_Smiles')['Library'].unique().reset_index()

# Filter to only keep compounds that appear in more than one library
overlap_df = overlap_df[overlap_df['Library'].apply(lambda libs: len(libs) > 1)]
print(overlap_df)



# BayesianGaussianMixture - probabilistic clustering approach
def calc_fp(smiles, fp_size=1024, radius=2):
    mol = Chem.MolFromSmiles(smiles)
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=fp_size)
    a = np.zeros((0,), dtype=np.float32)
    Chem.DataStructs.ConvertToNumpyArray(fp, a)
    return a

combined_df['fingerprint'] = combined_df['Canonical_Smiles'].apply(calc_fp)
fingerprint_array = np.array(list(combined_df['fingerprint']))



# Train Separate Bayesian Models per Library
libraries = combined_df['Library'].unique()
models = {}

for lib in libraries:
    indices = combined_df['Library'] == lib
    X_lib = fingerprint_array[indices].astype(float)
    
    # Train a BayesianGaussianMixture model on the fingerprints for this library
    bgmm = BayesianGaussianMixture(n_components=10, random_state=42)
    bgmm.fit(X_lib)
    models[lib] = bgmm


# Function to compute log-likelihood in batches
def compute_log_probs_in_batches(model, X, batch_size=10000, dtype=np.float32):
    n_samples = X.shape[0]
    log_probs_list = []
    for start in range(0, n_samples, batch_size):
        end = min(start + batch_size, n_samples)
        batch = X[start:end].astype(dtype)
        log_probs_list.append(model.score_samples(batch))
    return np.concatenate(log_probs_list)


score_matrix = pd.DataFrame(index=combined_df.index, columns=libraries)
    
for lib, model in models.items():
    # Compute log-likelihood scores in batches to avoid memory issues.
    log_probs = compute_log_probs_in_batches(model, fingerprint_array, batch_size=10000, dtype=np.float32)
    score_matrix[lib] = log_probs

# Normalize Scores so that scores are comparable across models.
scaler = MinMaxScaler()
score_matrix_scaled = pd.DataFrame(scaler.fit_transform(score_matrix), 
                                   index=score_matrix.index, 
                                   columns=score_matrix.columns)
score_matrix_scaled.to_csv('score_matrix_scaled_FourdLibs.csv', index=False)



# Apply t-SNE to the normalized score matrix
tsne = TSNE(n_components=2, random_state=42)
embedding_tsne = tsne.fit_transform(score_matrix_scaled)

# Add the t-SNE components back to the DataFrame
combined_df['tsne_score1'] = embedding_tsne[:, 0]
combined_df['tsne_score2'] = embedding_tsne[:, 1]

combined_df['Library'] = combined_df['Library'].replace({
    'ChemSpace_86K': 'ChemSpace (86K)',
    'FDA_approved': 'FDA Approved (1K)',
    'Kinase_targeted': 'Kinase Targeted (65K)',
    'Natural_products': 'Natural Products (42K)'
})


Paired = sns.color_palette("Paired").as_hex()
custom_palette = {
    "FDA Approved (1K)":      Paired[9],  # purple
    "Kinase Targeted (65K)":         Paired[7],  # orange
    "Natural Products (42K)":  Paired[5],  # red
    "ChemSpace (86K)":   Paired[1],  # blue
}

plt.figure(figsize=(8, 6))
sns.scatterplot(data=combined_df, 
                x='tsne_score1', 
                y='tsne_score2', 
                hue='Library', 
                palette=custom_palette, 
                s=10, 
                alpha=0.7)
# plt.title("t-SNE Visualization on Normalized Score Matrix")
plt.xlabel("t-SNE Component 1")
plt.ylabel("t-SNE Component 2")
plt.legend(loc='best')
# plt.legend(loc='lower right')
# plt.show()

# Save the plot
plt.savefig("bayesian_ecfp4_tsne_visualization_FourLibs.png", dpi=300, bbox_inches='tight')


