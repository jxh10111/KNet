import csv
import threading
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, Draw
from sklearn.cluster import KMeans
from sklearn.mixture import BayesianGaussianMixture
from k_means_constrained import KMeansConstrained
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

# -------------------------------
# Config
# -------------------------------
files = {
    'Library A': 'data/ChEMBL31_KKBQ122_dc_keras_model_328targets_cutoff7_KNet_Pred_Enamine_FDA_approved_Drugs_1122cmpds.csv',
    'Library B': 'data/ChEMBL31_KKBQ122_dc_keras_model_328targets_cutoff7_KNet_Pred_Enamine_Kinase_Targeted_Library_64880_cpds.csv',
    'Library C': 'data/ChEMBL31_KKBQ122_dc_keras_model_328targets_cutoff7_KNet_Pred_Natural_Products_MEGx_NATx_41820_compounds.csv',
    'Library D': 'data/ChEMBL34_KKBQ224_dc_keras_model_121924_328targets_cutoff7_KNet_Pred_ChemSpace0001_28K.csv'
}
activity_cutoff = 5
micro_size = 10
n_macro = 10

# -------------------------------
# Utility Functions
# -------------------------------
def gini_coefficient(arr):
    x = np.sort(np.array(arr, dtype=float))
    if x.sum() == 0:
        return 0.0
    n = x.size
    idx = np.arange(1, n + 1)
    return (2 * (idx * x).sum() - (n + 1) * x.sum()) / (n * x.sum())

mg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
def smiles_to_fp(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = mg.GetFingerprint(mol)
    arr = np.zeros((fp.GetNumBits(),), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

# -------------------------------
# Spinner for long steps
# -------------------------------
class Spinner:
    def __init__(self, interval=10):
        self.interval = interval
        self.running = False

    def start(self):
        self.running = True
        threading.Thread(target=self._spin, daemon=True).start()

    def _spin(self):
        while self.running:
            print(".", end="", flush=True)
            time.sleep(self.interval)

    def stop(self):
        self.running = False

# -------------------------------
# Enable tqdm for pandas
# -------------------------------
tqdm.pandas(desc="SMILES→Fingerprint")

# -------------------------------
# Load & Preprocess
# -------------------------------
records = []
for lib_name, path in files.items():
    # detect delimiter
    with open(path, newline='') as f:
        sample = f.read(4096)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=[',','\t',';'])
            delim = dialect.delimiter
        except csv.Error:
            delim = ','
    df = pd.read_csv(path, sep=delim)
    df['library'] = lib_name

    # rename first column to Canonical_Smiles
    first_col = df.columns[0]
    df = df.rename(columns={first_col: 'Canonical_Smiles'})

    # identify kinase columns
    kin_cols = [c for c in df.columns if c not in ('Canonical_Smiles','library')]

    # coerce to numeric & drop invalid
    kin_num = df[kin_cols].apply(pd.to_numeric, errors='coerce')
    valid = kin_num.notnull().all(axis=1)
    df = df.loc[valid].reset_index(drop=True)
    kin_num = kin_num.loc[valid].reset_index(drop=True)

    # binarize & compute GINI
    kin_bin = (kin_num >= activity_cutoff).astype(int)
    gini_s = kin_bin.apply(gini_coefficient, axis=1)

    # assemble cleaned record
    rec = pd.concat([
        df[['Canonical_Smiles','library']].reset_index(drop=True),
        kin_bin.reset_index(drop=True),
        pd.Series(gini_s, name='gini')
    ], axis=1)
    records.append(rec)

data = pd.concat(records, ignore_index=True).copy()

# -------------------------------
# Selectivity
# -------------------------------
medians = data.groupby('library')['gini'].median().to_dict()
data['selectivity'] = data.apply(
    lambda r: 'Selective' if r['gini'] >= medians[r['library']] else 'Unselective',
    axis=1
)

# -------------------------------
# Fingerprints with progress bar
# -------------------------------
data['fp'] = data['Canonical_Smiles'].progress_map(smiles_to_fp)
data = data[data['fp'].notnull()].reset_index(drop=True)
fp_matrix = np.vstack(data['fp'].values)
N = fp_matrix.shape[0]

# -------------------------------
# Stage 1: Micro‑clustering (exactly 10 compounds each)
# -------------------------------
print("Starting micro‑clustering (size=10); this may take a while", end="", flush=True)
spinner = Spinner(interval=10)  
spinner.start()

k_micro = KMeansConstrained(
    n_clusters=N // micro_size,
    size_min=micro_size,
    size_max=micro_size,
    random_state=42,
    verbose=1
)
micro_labels = k_micro.fit_predict(fp_matrix)

spinner.stop()
print("\nMicro‑clustering complete.\n")
data['micro_cluster'] = micro_labels

# compute micro‑centroids & GINI
centroids = np.vstack([
    fp_matrix[micro_labels == i].mean(axis=0)
    for i in range(N // micro_size)
])
micro_gini = np.array([
    data.loc[data['micro_cluster'] == i, 'gini'].mean()
    for i in range(N // micro_size)
])

# -------------------------------
# Stage 2: Macro‑clustering (exactly 10 clusters)
# -------------------------------
print("Starting macro‑clustering (10 clusters)…")
k_macro = KMeans(
    n_clusters=n_macro,
    random_state=42,
    n_init='auto',
    verbose=1
)
macro_labels = k_macro.fit_predict(centroids)
data['cluster'] = macro_labels[micro_labels]
print("Macro‑clustering complete.\n")

# -------------------------------
# Plots 1–4: Grouped Barplots per Library
# -------------------------------
for lib in files:
    sub = data[data['library'] == lib]
    agg = sub.groupby(['cluster','selectivity'])['gini'] \
             .mean().unstack(fill_value=0).sort_index()

    fig, ax = plt.subplots(figsize=(8,4))
    agg.plot(kind='bar', ax=ax)
    ax.set_title(f'{lib}: Mean GINI by Macro‑cluster & Selectivity')
    ax.set_xlabel('Macro‑cluster')
    ax.set_ylabel('Mean GINI')
    ax.legend(title='Selectivity', loc='upper right')

    for c in agg.index:
        fps = np.vstack(sub[sub['cluster'] == c]['fp'].tolist())
        center_fp = fps.mean(axis=0)
        dists = sub[sub['cluster'] == c]['fp'].apply(lambda x: np.linalg.norm(x - center_fp))
        idx = dists.idxmin()
        mol = Chem.MolFromSmiles(sub.loc[idx, 'Canonical_Smiles'])
        img = Draw.MolToImage(mol, size=(100,100))
        ab = AnnotationBbox(
            OffsetImage(img, zoom=1.0),
            (c, 0),
            frameon=False,
            box_alignment=(0.5, -0.2)
        )
        ax.add_artist(ab)

    plt.tight_layout()
    plt.show()

# -------------------------------
# Plot 5: Similarity vs. GINI Scatterplot
# -------------------------------
bgmm = BayesianGaussianMixture(n_components=2, random_state=42)
proba = bgmm.fit_predict_proba(centroids)
sim = proba[:, 0]

cc = pd.DataFrame({
    'micro': np.arange(N // micro_size),
    'macro': macro_labels,
    'gini': micro_gini,
    'sim': sim
})
plot_df = cc.groupby('macro').agg({'gini':'mean','sim':'mean'}).reset_index()

fig, ax = plt.subplots(figsize=(6,6))
ax.scatter(
    plot_df['sim'], plot_df['gini'],
    c=plot_df['macro'], cmap='tab10', s=100
)
for _, row in plot_df.iterrows():
    ax.text(
        row['sim'], row['gini'],
        f"C{int(row['macro'])}",
        ha='center', va='center',
        fontsize=8, color='white'
    )

ax.set_title('Macro‑cluster Similarity vs. Mean GINI')
ax.set_xlabel('Similarity (BGMM P component 0)')
ax.set_ylabel('Mean GINI')
plt.tight_layout()
plt.show()
