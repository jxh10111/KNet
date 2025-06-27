import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity


#### Data processing ####

# Load data
df = pd.read_csv('ChEMBL34_KKBQ224_dc_keras_model_121924_328targets_cutoff7_KNet_Pred_ChemSpace0003_86K.csv', index_col=0)
df = pd.read_csv('ChEMBL34_KKBQ224_dc_keras_model_328targets_cutoff7_KNet_Pred_Enamine_FDA_approved_Drugs_1122cmpds.csv', index_col=0)
df = pd.read_csv('ChEMBL34_KKBQ224_dc_keras_model_328targets_cutoff7_KNet_Pred_Enamine_Kinase_Targeted_Library_64880_cpds.csv', index_col=0)
df = pd.read_csv('ChEMBL34_KKBQ224_dc_keras_model_328targets_cutoff7_KNet_Pred_Natural_Products_MEGx_NATx_41820_compounds.csv', index_col=0)

smiles = df.index.tolist()
binary = (df >= 0.5)
labels = np.where(binary.any(axis=1), 'selective', 'unselective')

# Generate 1024-bit ECFP4 fingerprints
fps = []
valid = []
for i, smi in enumerate(smiles):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        arr = np.array(list(map(int, fp.ToBitString())), dtype=int)
        fps.append(arr)
        valid.append(i)
fps = np.vstack(fps)
labels = labels[valid]
valid_smiles = [smiles[i] for i in valid]

# Global KMeans → chemotype clusters
n_clusters = max(2, len(fps) // 10)
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
cluster_ids = kmeans.fit_predict(fps)

# Compute per-cluster selectivity fraction
cluster_props = {
    cid: (labels[cluster_ids==cid] == 'selective').mean()
    for cid in range(n_clusters)
}



#### Panel 1 - Barplot ####

binary_valid = binary.loc[valid_smiles]

records = []
for cid in np.unique(cluster_ids):
    mask = (cluster_ids == cid)
    smiles_c = np.array(valid_smiles)[mask]
    n_actives = (labels[mask] == 'selective').sum()
    if n_actives == 0:
        continue
    
    sub = binary_valid.loc[smiles_c]
    n_distinct_kinases = sub.any(axis=0).sum()
    n_total = mask.sum()
    n_inactives = n_total - n_actives
    sub = binary_valid.loc[smiles_c]
    kinase_hits_bool = sub.any(axis=0)
    n_distinct       = int(kinase_hits_bool.sum())
    kinase_list = list(kinase_hits_bool.index[kinase_hits_bool])
    
    records.append({
        'cluster':      cid,
        'N_actives':    n_actives,
        'N_inactives':   n_inactives,
        'N_total':      n_total,
        'N_kinases':    n_distinct_kinases,
        'Ratio':        n_distinct_kinases / n_actives,
        'Kinase_List': kinase_list
    })

df_stats = pd.DataFrame(records).set_index('cluster')

# Sort by number of kinases predicted by cluster (high→low)
df_stats = df_stats.sort_values('N_kinases', ascending=False)

# Plot bars + line on twin y‐axes ---
fig, ax1 = plt.subplots(figsize=(14,8), dpi=300)

ax1.bar(
    df_stats.index.astype(str),
    df_stats['N_actives'],
    color='tomato',
    label='Predicted active'
)
ax1.bar(
    df_stats.index.astype(str),
    df_stats['N_inactives'],
    bottom=df_stats['N_actives'],
    color='skyblue',
    label='Predicted inactive'
)
ax1.set_xlabel('Cluster ID')
ax1.set_ylabel('Number of compounds')
ax1.tick_params(axis='x', rotation=45, labelsize=4)

ax1.set_yticks(np.arange(0, df_stats[['N_actives', 'N_inactives']].sum(axis=1).max() + 5, 5))
ax1.yaxis.grid(True, linestyle='--', linewidth=0.5)

from matplotlib.ticker import MaxNLocator
ax2 = ax1.twinx()
ln = ax2.plot(
    df_stats.index.astype(str),
    df_stats['N_kinases'],
    color='navy',
    marker='o',
    markersize=2, 
    linewidth=1,
    label='Kinases predicted by cluster',
    zorder=1  # Lower zorder
)
ax2.set_ylabel('Total number of kinases predicted by cluster')
ax2.yaxis.set_major_locator(MaxNLocator(integer=True))

# Combine legends
bars, bar_labels = ax1.get_legend_handles_labels()
lines, line_labels = ax2.get_legend_handles_labels()
legend = ax2.legend(
    bars + lines,
    bar_labels + line_labels,
    loc='upper right',
    fontsize=10
    )

# plt.title('Cluster promiscuity: # actives vs. kinase/active ratio', fontsize=12)
# plt.title('ChemSpace (86K)') 
# plt.title('FDA Approved (1K)')
# plt.title('Kinase Targeted (65K)')
# plt.title('Natural Products (42K)')
plt.tight_layout(pad=2)
# plt.show()

# plt.savefig('ChemSpace_86K_selectivity_N_kinases_barplot.png', bbox_inches='tight', pad_inches=0)
# plt.savefig('FDA_Approved_selectivity_N_kinases_barplot.png', bbox_inches='tight', pad_inches=0)
# plt.savefig('Kinase_Targeted_selectivity_N_kinases_barplot.png', bbox_inches='tight', pad_inches=0)
# plt.savefig('Natural_Product_selectivity_N_kinases_barplot.png', bbox_inches='tight', pad_inches=0)

# df_stats.to_csv('ChemSpace_86K_cluster_stats.csv')
# df_stats.to_csv('FDA_Approved_cluster_stats.csv')
# df_stats.to_csv('Kinase_Targeted_cluster_stats.csv')
# df_stats.to_csv('Natural_Product_cluster_stats.csv')




#### Panel 2 - TSNE plot ####

n_bits = fps.shape[1]
cluster_bitfps = []
for cid in range(n_clusters):
    mask = (cluster_ids == cid)
    if not mask.any():
        cluster_bitfps.append(np.zeros(n_bits, dtype=int))
    else:
        or_fp = (fps[mask].sum(axis=0) > 0).astype(int)
        cluster_bitfps.append(or_fp)
cluster_bitfps = np.vstack(cluster_bitfps)  # shape = (n_clusters, 1024)

# Run t-SNE on bit‐vectors with Jaccard (Tanimoto) distance
tsne = TSNE(
    n_components=2,
    metric='jaccard',
    init='pca',
    learning_rate='auto',
    random_state=42
)
cent_emb = tsne.fit_transform(cluster_bitfps)


# Randomily annotate 6 red circles with diversity in the TSNE plot
from matplotlib import cm, colors as mcolors
from adjustText import adjust_text  # pip install adjustText

fracs = np.array([cluster_props[cid] for cid in range(n_clusters)])
sizes = np.array([(cluster_ids==cid).sum() for cid in range(n_clusters)])*5

# Set up colormap and normalization
cmap = cm.get_cmap('Reds')
norm = mcolors.Normalize(vmin=0, vmax=1)
edge_colors = cmap(norm(fracs))  # fracs = array of selectivity fractions

# Plot hollow circles
fig, ax = plt.subplots(figsize=(6,6), dpi=300)
ax.scatter( 
    cent_emb[:,0], cent_emb[:,1],
    s=sizes,
    facecolors='none',       # no fill
    edgecolors=edge_colors,  # colored outlines
    linewidths=0.4,
    alpha=0.9,
    zorder=1
)

# Add a colorbar
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # dummy mappable
cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Selectivity Fraction')

# Arrow‐point and label each all‐selective cluster or clusters with at least one active cmpd
eligible_cids = [cid for cid, f in cluster_props.items() if f>=0.5] 

np.random.seed(42) 
selected_cids = np.random.choice(eligible_cids, size=6, replace=False)

# compute a reference split on the x-axis
xs = cent_emb[:,0]
median_x = np.median(xs)
dx = 5   # horizontal label offset
dy = -6  # vertical label offset

# dx=5, dy=-6 (kinase-targetd and natural products and chemspace 86k)
# dx=0.8, dy=-0.6 (FDA-approved)

for cid in selected_cids: # all_sel_cids
    x, y = cent_emb[cid]
    if x > median_x:
        text_x = x + dx
        ha = 'left'
    else:
        text_x = x - dx
        ha = 'right'
        
    ax.annotate(
        f"Cluster {cid}",
        xy=(x, y),
        xytext=(text_x, y + dy),
        ha=ha,
        va='center',
        arrowprops=dict(
            arrowstyle="->",
            color="black",
            linewidth=1.2
        ),
        fontsize=5,
        backgroundcolor="none",
        zorder=2
    )

ax.set_xlabel('t-SNE 1')
ax.set_ylabel('t-SNE 2')
plt.tight_layout(pad=2)
# plt.show()

# fig.savefig('ChemSpace_86K_centroid_tsne_hollow_tanimoto_similarity.png', bbox_inches='tight', pad_inches=0)
# fig.savefig('FDA_Approved_centroid_tsne_hollow_tanimoto_similarity.png', bbox_inches='tight', pad_inches=0)
# fig.savefig('Kinase_Targeted_centroid_tsne_hollow_tanimoto_similarity.png', bbox_inches='tight', pad_inches=0)
# fig.savefig('Natural_Product_centroid_tsne_hollow_tanimoto_similarity.png', bbox_inches='tight', pad_inches=0)
# plt.close(fig)




#### Panel 3 - chemotype structures ####

# Draw chemotype for all_selective clusters
from matplotlib.patches import Rectangle
from pubchempy import get_compounds

# For the randomly selected 6 clusters
for cid in selected_cids:
    members = np.where(cluster_ids == cid)[0]
    centroid = fps[members].mean(axis=0)
    dists = np.linalg.norm(fps[members] - centroid, axis=1)
    medoid_idx = members[np.argmin(dists)]
    
    medoid_smi = smiles[medoid_idx]
    
    mol = Chem.MolFromSmiles(smiles[medoid_idx])
    img = Draw.MolToImage(mol, size=(800, 800)) 

    fig, ax = plt.subplots(figsize=(5, 5.5), dpi=600)
    ax.imshow(img)
    ax.axis('off')
    
    ax.set_xlim(0, img.size[0])
    ax.set_ylim(img.size[1], 0)
    
    # Add black frame around image
    rect = Rectangle((0, 0), 1, 1, linewidth=2, edgecolor='black', facecolor='none', transform=ax.transAxes, zorder=10)
    ax.add_patch(rect)
    
    # title + subtitle
    ax.set_title(f'Cluster {cid}', fontsize=14)
    # Get stats from df_stats
    n_act = df_stats.loc[cid, 'N_actives']
    n_inact = df_stats.loc[cid, 'N_inactives']
    n_kin = df_stats.loc[cid, 'N_kinases']
    
    plt.figtext(
        0.5, 0.04,
        f'# of Actives: {n_act}   |   # of Inactives: {n_inact}   |   # of Targeted Kinases: {n_kin}',
        horizontalalignment='center',
        va='bottom',
        transform=ax.transAxes,
        fontsize=10
    )
    
    plt.tight_layout()
    filename = f'cluster_{cid:02d}_medoid_v2.png'
    plt.savefig(filename, dpi=600, bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)

    print(f"Saved {filename}")
    
    

