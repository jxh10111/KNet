import csv
import threading
import time
import pickle
import multiprocessing

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
# Utility Functions
# -------------------------------
def gini_coefficient(arr):
    """Compute the GINI coefficient of a 1D array."""
    x = np.sort(np.array(arr, dtype=float))
    if x.sum() == 0:
        return 0.0
    n = x.size
    idx = np.arange(1, n + 1)
    return (2 * (idx * x).sum() - (n + 1) * x.sum()) / (n * x.sum())


mg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
def smiles_to_fp(smi):
    """Convert a SMILES string to a 1024‑bit Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = mg.GetFingerprint(mol)
    arr = np.zeros((fp.GetNumBits(),), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


class Spinner:
    """Simple console spinner for long‑running steps."""
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


def load_and_preprocess_single(lib_name, path, activity_cutoff):
    """
    Load one CSV, binarize actives/inactives, compute per‑compound GINI and selectivity,
    and compute fingerprints. Returns a DataFrame `data` for just that library,
    plus the fingerprint matrix.
    """
    # 1) detect delimiter
    with open(path, newline='') as f:
        sample = f.read(4096)
        try:
            delim = csv.Sniffer().sniff(sample, delimiters=[',','\t',';']).delimiter
        except csv.Error:
            delim = ','
    df = pd.read_csv(path, sep=delim)
    df['library'] = lib_name

    # 2) rename first column to Canonical_Smiles
    first_col = df.columns[0]
    df = df.rename(columns={first_col: 'Canonical_Smiles'})

    # 3) identify kinase columns & coerce to numeric
    kin_cols = [c for c in df.columns if c not in ('Canonical_Smiles','library')]
    kin_num = df[kin_cols].apply(pd.to_numeric, errors='coerce')

    # ──────────────────────────────────────────────────────────────────────────────
    # DEBUG: print out any kinase values ≥ activity_cutoff
    mask = kin_num >= activity_cutoff
    if mask.values.any():
        print(f"\nFound {mask.values.sum()} kin_num entries ≥ {activity_cutoff}:")
        # stack() will give you (row_index, column) → value
        high_vals = kin_num[mask].stack()
        print(high_vals.head(20))   # show up to 20 hits
        if len(high_vals) > 20:
            print("…")
    else:
        print(f"\nNo kin_num entries ≥ {activity_cutoff} were found.")
    # ──────────────────────────────────────────────────────────────────────────────

    valid = kin_num.notnull().all(axis=1)
    df = df.loc[valid].reset_index(drop=True)
    kin_num = kin_num.loc[valid].reset_index(drop=True)

    # 4) binarize & compute GINI
    kin_bin = (kin_num >= activity_cutoff).astype(int)
    gini_s = kin_bin.apply(gini_coefficient, axis=1)

    # 5) assemble into one DataFrame
    data = pd.concat([
        df[['Canonical_Smiles','library']].reset_index(drop=True),
        kin_bin.reset_index(drop=True),
        pd.Series(gini_s, name='gini')
    ], axis=1)

    # 6) label selectivity using the _median_ GINI for this library
    threshold = 0.5 # data['gini'].median()
    data['selectivity'] = data['gini'].apply(lambda g: 'Selective' if g >= threshold else 'Unselective')

    # 7) compute fingerprints
    tqdm.pandas(desc=f"{lib_name}: SMILES→FP")
    data['fp'] = data['Canonical_Smiles'].progress_map(smiles_to_fp)
    data = data[data['fp'].notnull()].reset_index(drop=True)

    # 8) build fingerprint matrix
    fp_matrix = np.vstack(data['fp'].values)

    return data, fp_matrix


def do_micro_clustering(fp_matrix, micro_size):
    """Perform size‑constrained micro‑clustering and return labels & centroids."""
    N = fp_matrix.shape[0]
    n_micro = max(1, N // micro_size)
    size_max = int(np.ceil(N / n_micro))

    print(f" → Micro‑clustering {N} compounds into {n_micro} clusters", end="", flush=True)
    spinner = Spinner(interval=10)
    spinner.start()

    k_micro = KMeansConstrained(
        n_clusters=n_micro,
        size_min=micro_size,
        size_max=size_max,
        random_state=42,
        verbose=1,
        n_jobs=1
    )
    labels = k_micro.fit_predict(fp_matrix)

    spinner.stop()
    print("  done.")

    centroids = np.vstack([
        fp_matrix[labels == i].mean(axis=0)
        for i in range(n_micro)
    ])
    return labels, centroids


def do_macro_clustering(centroids, n_macro):
    """Cluster the micro‑centroids into exactly `n_macro` macro‑clusters."""
    print(f" → Macro‑clustering into {n_macro} clusters", end="", flush=True)
    k_macro = KMeans(n_clusters=n_macro, random_state=42, n_init='auto', verbose=1)
    labels = k_macro.fit_predict(centroids)
    print("  done.\n")
    return labels


def prepare_plot1to4_data(data):
    """Aggregate mean GINI and pick representative SMILES per macro‑cluster & selectivity."""
    agg = (
        data.groupby(['cluster','selectivity'])['gini']
            .mean()
            .unstack(fill_value=0)
            .sort_index()
    )

    rep_smiles = {}
    for c in agg.index:
        fps = np.vstack(data.loc[data['cluster']==c, 'fp'].tolist())
        center = fps.mean(axis=0)
        dists = data['fp'].apply(lambda x: np.linalg.norm(x - center))
        idx = dists.idxmin()
        rep_smiles[c] = data.loc[idx, 'Canonical_Smiles']

    return agg, rep_smiles


def plot_grouped_barplots(lib_name, data, agg, rep_smiles):
    """Grouped barplot with auto‐sized insets and outside legend."""
    fig, ax = plt.subplots(figsize=(8,4))

    # leave extra room on the right for the legend
    fig.subplots_adjust(right=0.8)

    # draw the bars
    agg.plot(kind='bar', ax=ax)
    ax.set_title(f'{lib_name}: Mean GINI by Macro‑cluster & Selectivity')
    ax.set_xlabel('Macro‑cluster')
    ax.set_ylabel('Mean GINI')

    # put legend outside at upper left of the padded area
    ax.legend(
        title='Selectivity',
        loc='upper left',
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0
    )

    # compute zoom so that for N bars, zoom ~ 2/N but never below 0.3
    N = len(agg.index)
    zoom = max(0.3, min(0.6, 2.0 / N))

    for c in agg.index:
        smiles = rep_smiles[c]
        mol    = Chem.MolFromSmiles(smiles)
        img    = Draw.MolToImage(mol, size=(80,80))  # slightly smaller base

        ab = AnnotationBbox(
            OffsetImage(img, zoom=zoom),
            (c, 0),
            frameon=False,
            # lift them further above the axis so they clear the bars
            box_alignment=(0.5, -0.2)
        )
        ax.add_artist(ab)

    plt.tight_layout()
    plt.show()




def prepare_plot5_data(centroids, micro_gini, macro_labels):
    """Compute similarity vs. mean GINI for a single library."""
    bgmm = BayesianGaussianMixture(n_components=2, random_state=42)
    bgmm.fit(centroids)
    sim = bgmm.predict_proba(centroids)[:, 0]

    cc = pd.DataFrame({
        'macro': macro_labels,
        'gini' : micro_gini,
        'sim'  : sim
    })
    return cc.groupby('macro').agg({'gini':'mean','sim':'mean'}).reset_index()


def plot_scatter(lib_name, plot_df):
    """Draw the similarity vs. GINI scatter for a single library."""
    fig, ax = plt.subplots(figsize=(6,6))
    ax.scatter(
        plot_df['sim'], plot_df['gini'], c=plot_df['macro'], cmap='tab10', s=100
    )
    for _, row in plot_df.iterrows():
        ax.text(
            row['sim'], row['gini'], f"C{int(row['macro'])}",
            ha='center', va='center', fontsize=8, color='white'
        )

    ax.set_title(f'{lib_name}: Similarity vs. Mean GINI')
    ax.set_xlabel('Similarity (BGMM P comp 0)')
    ax.set_ylabel('Mean GINI')
    plt.tight_layout()
    plt.show()


def main():
    multiprocessing.set_start_method('spawn', force=True)

    # -------------------------------
    # Config
    # -------------------------------
    files = {
        'Library A': 'data/ChEMBL31_KKBQ122_dc_keras_model_328targets_cutoff7_KNet_Pred_Enamine_FDA_approved_Drugs_1122cmpds.csv',
        'Library B': 'data/ChEMBL31_KKBQ122_dc_keras_model_328targets_cutoff7_KNet_Pred_Enamine_Kinase_Targeted_Library_64880_cpds.csv',
        'Library C': 'data/ChEMBL31_KKBQ122_dc_keras_model_328targets_cutoff7_KNet_Pred_Natural_Products_MEGx_NATx_41820_compounds.csv',
        'Library D': 'data/ChEMBL34_KKBQ224_dc_keras_model_121924_328targets_cutoff7_KNet_Pred_ChemSpace0001_28K.csv'
    }
    activity_cutoff = 0.5
    micro_size      = 10
    n_macro         = 10

    for lib_name, path in files.items():
        print(f"\n=== Processing {lib_name} ===")

        # 1) Load & preprocess
        data, fp_matrix = load_and_preprocess_single(lib_name, path, activity_cutoff)

        # print the library‐level median GINI to the console:
        median_gini = data['gini'].median()
        print(f"{lib_name} median GINI = {median_gini:.3f}")

        # 2) Micro‑clustering
        micro_labels, centroids = do_micro_clustering(fp_matrix, micro_size)
        data['micro_cluster'] = micro_labels

        micro_gini = data.groupby('micro_cluster')['gini'].mean().values

        # 3) Macro‑clustering
        macro_labels = do_macro_clustering(centroids, n_macro)
        data['cluster'] = macro_labels[micro_labels]

        # 4) Prepare & pickle Plot1-4 data
        agg, rep_smiles = prepare_plot1to4_data(data)
        with open(f'data/{lib_name}_plot1to4.pkl', 'wb') as f:
            pickle.dump({'agg': agg, 'smiles': rep_smiles}, f)

        # 5) Draw Plot1-4
        plot_grouped_barplots(lib_name, data, agg, rep_smiles)

        # 6) Prepare & pickle Plot5 data
        plot5_df = prepare_plot5_data(centroids, micro_gini, macro_labels)
        with open(f'data/{lib_name}_plot5.pkl', 'wb') as f:
            pickle.dump(plot5_df, f)

        # 7) Draw Plot5
        plot_scatter(lib_name, plot5_df)


if __name__ == "__main__":
    main()
