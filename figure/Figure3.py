from matplotlib.ticker import FuncFormatter

def hit_distribution_from_long_with_totals(df_long: pd.DataFrame, cutoff: float, total_compounds: int) -> pd.DataFrame:
    df_active = df_long[df_long['pred_value'] >= cutoff]
    if df_active.empty:
        return pd.DataFrame({
            "Hits": hits_order,
            "Count": [0]*len(hits_order),
            "PercentAmongActives": [0.0]*len(hits_order),
            "PercentOfAll": [0.0]*len(hits_order),
            "Cutoff": [cutoff]*len(hits_order),
            "ActiveCount": [0]*len(hits_order),
            "TotalCompounds": [total_compounds]*len(hits_order),
            "ActiveRate": [0.0]*len(hits_order),
        })

    # per compound, how many unique kinases pass cutoff
    hit_counts = (
        df_active.groupby('Canonical_Smiles')['kinase']
                 .nunique()
                 .rename('hit_count')
    )

    # keep compounds with >=1 hit (actives)
    active_counts = hit_counts[hit_counts > 0]

    # bin into 1,2,3,4,5,6+ hit count
    bins = active_counts.apply(lambda x: str(x) if x <= 5 else "6+")
    vc = bins.value_counts().reindex(hits_order_raw, fill_value=0)

    total_active = int(vc.sum())
    perc_among_actives = (vc / total_active * 100.0) if total_active > 0 else (vc * 0.0)
    perc_of_all = (vc / total_compounds * 100.0) if total_compounds > 0 else (vc * 0.0)
    active_rate = (total_active / total_compounds * 100.0) if total_compounds > 0 else 0.0

    out = pd.DataFrame({
        "Hits": [hit_label_map[h] for h in vc.index],
        "Count": vc.values.astype(int),
        "PercentAmongActives": perc_among_actives.values,
        "PercentOfAll": perc_of_all.values,
        "Cutoff": cutoff,
        "ActiveCount": total_active,
        "TotalCompounds": total_compounds,
        "ActiveRate": active_rate,
    })
    return out


all_dist = []
for lib_name, df_long in libraries_long.items():
    df_long = df_long.copy()
    df_long['pred_value'] = pd.to_numeric(df_long['pred_value'], errors='coerce')
    df_long = df_long.dropna(subset=['pred_value', 'Canonical_Smiles', 'kinase'])
    total_compounds = df_long['Canonical_Smiles'].nunique()

    for c in cutoffs:
        dist = hit_distribution_from_long_with_totals(df_long, c, total_compounds)
        dist['Library'] = lib_name
        all_dist.append(dist)

dist_df = pd.concat(all_dist, ignore_index=True)
dist_df.to_csv("Hit_Distribution_with_ActualRates.csv", index=False)


# Choose y-axis mode: 'percent' for actual hit rate by bin (% of ALL compounds), or 'count' for # compounds
Y_MODE = 'count'   # 'percent' or 'count'
y_col  = 'PercentOfAll' if Y_MODE == 'percent' else 'Count'
y_label = "Percent of all compounds (hit rate)" if Y_MODE == 'percent' else "Number of Predicted Active Compounds"

libs_in_order = list(libraries_long.keys())
n_libs = len(libs_in_order)

fig, axes = plt.subplots(1, n_libs, figsize=(4.6*n_libs, 5.8), sharey=(Y_MODE=='percent'))
if n_libs == 1:
    axes = [axes]

for ax, lib in zip(axes, libs_in_order):
    sub = dist_df[dist_df['Library'] == lib].copy()

    # pivot to have rows=cutoff, cols=Hits, values=y_col (PercentOfAll or Count)
    piv = sub.pivot_table(index='Cutoff', columns='Hits', values=y_col, aggfunc='sum').reindex(index=cutoffs)

    # ensure all hit bins exist as columns
    for h in hits_order:
        if h not in piv.columns:
            piv[h] = 0.0
    piv = piv[hits_order]  # column order

    # data for annotation: overall active count & rate per cutoff
    agg = (sub.groupby('Cutoff')[['ActiveCount','TotalCompounds','ActiveRate']]
              .first()
              .reindex(index=cutoffs))

    # stacked bars
    bottoms = np.zeros(len(piv))
    x = np.arange(len(piv.index))
    bar_containers = []
    for h in hits_order:
        bars = ax.bar(
            x, piv[h].values, bottom=bottoms,
            label=h, color=hits_colors[h], edgecolor='white', linewidth=0.4
        )
        bar_containers.append(bars)
        bottoms += piv[h].values

    # annotate each bar with overall hit rate and (actives/total)
    for i, cutoff_val in enumerate(piv.index):
        ar = agg.loc[cutoff_val, 'ActiveRate']
        ac = int(agg.loc[cutoff_val, 'ActiveCount'])
        tc = int(agg.loc[cutoff_val, 'TotalCompounds'])
        label = f"{ar:.1f}% ({ac}/{tc})" if Y_MODE=='percent' else f"{ar:.1f}%"
        ax.text(i, bottoms[i] + (0.8 if Y_MODE=='percent' else max(1, bottoms.max()*0.01)),
                label, ha='center', va='bottom', fontsize=8)

    ax.set_title(lib)
    ax.set_xticks(x)
    ax.set_xticklabels([str(c) for c in cutoffs])
    ax.set_xlabel("Cutoff")
    ax.set_ylabel(y_label)
    if lib.startswith("ChemSpace"):
        ax.set_ylabel("Number of Predicted Active Compounds (in Millions)")
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f"{y/1_000_000:.1f}"))
        ax.yaxis.offsetText.set_visible(False)

    ax.grid(axis='y', linestyle='--', alpha=0.35)
    if Y_MODE == 'percent':
        ax.set_ylim(0, max(105, bottoms.max()*1.15))

handles, labels = axes[-1].get_legend_handles_labels()
fig.legend(handles, labels, title="Hit count", loc='center left', bbox_to_anchor=(1.0, 0.5))
fig.tight_layout()
plt.savefig(f"Figure3_ActualHitRate_{'Percent' if Y_MODE=='percent' else 'Counts'}.png", dpi=300, bbox_inches="tight")
plt.show()
