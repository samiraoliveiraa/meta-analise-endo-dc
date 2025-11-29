import pandas as pd
import numpy as np
import shap
from scipy.spatial.distance import pdist, squareform
from skbio.stats.ordination import pcoa
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import warnings

from sklearn.linear_model import LogisticRegression
from interpret.glassbox import ExplainableBoostingClassifier
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, average_precision_score, roc_auc_score, precision_recall_curve, roc_curve, auc, classification_report, confusion_matrix
from sklearn.utils.class_weight import compute_sample_weight
from skbio.stats.composition import clr
from sklearn.ensemble import RandomForestClassifier

def pcoa_comparativa(estudo):

    tax_original = pd.read_csv(f"taxonomia_{estudo}.csv")
    tax_corrigida = pd.read_csv(f"taxonomia_{estudo}_ajustada.csv")
    metadata = pd.read_csv(f"{estudo}_metadata.tsv", sep='\t')
    
    metadata = metadata.rename(columns={"sample-id": "sample_id"})
    if "Unnamed: 0" in tax_corrigida.columns:
        tax_corrigida = tax_corrigida.rename(columns={"Unnamed: 0": "index"})
    
    tax_original_numeric = tax_original.iloc[:, 1:-3].apply(pd.to_numeric, errors='coerce').fillna(0)
    tax_corrigida_numeric = tax_corrigida.iloc[:, 1:].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    tax_original_numeric.index = tax_original["index"]
    tax_corrigida_numeric.index = tax_corrigida["index"]
    
    samples_common = tax_original_numeric.index.intersection(
    tax_corrigida_numeric.index).intersection(metadata['sample_id'])
    tax_original_numeric = tax_original_numeric.loc[samples_common]
    tax_corrigida_numeric = tax_corrigida_numeric.loc[samples_common]
    metadata_sub = metadata[metadata['sample_id'].isin(samples_common)]
    
    def bray_curtis_matrix(df):
        return pd.DataFrame(
            squareform(pdist(df.values, metric='braycurtis')),
            index=df.index,
            columns=df.index
        )
    
    bray_original = bray_curtis_matrix(tax_original_numeric)
    bray_corrigida = bray_curtis_matrix(tax_corrigida_numeric)
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pcoa_original_res = pcoa(bray_original)
        pcoa_corrigida_res = pcoa(bray_corrigida)
    
    df_original_plot = pd.DataFrame({
        "PCoA1": pcoa_original_res.samples['PC1'],
        "PCoA2": pcoa_original_res.samples['PC2'],
        "sample_id": tax_original_numeric.index
    }).merge(metadata_sub, on="sample_id")
    
    df_corrigida_plot = pd.DataFrame({
        "PCoA1": pcoa_corrigida_res.samples['PC1'],
        "PCoA2": pcoa_corrigida_res.samples['PC2'],
        "sample_id": tax_corrigida_numeric.index
    }).merge(metadata_sub, on="sample_id")
    
    var_exp_orig = pcoa_original_res.proportion_explained * 100
    var_exp_corr = pcoa_corrigida_res.proportion_explained * 100
    
    unique_groups = df_original_plot['group'].nunique()
    palette = sns.color_palette("tab10", n_colors=unique_groups)
    
    sns.set(style="white", font_scale=1.2)
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    
    df_original_plot = df_original_plot.rename(columns={'group': 'Grupo', 'study_id': 'Estudo'})
    df_corrigida_plot = df_corrigida_plot.rename(columns={'group': 'Grupo', 'study_id': 'Estudo'})
    
    sns.scatterplot(
        data=df_original_plot,
        x="PCoA1", y="PCoA2",
        hue="Grupo", style="Estudo",
        s=120, alpha=0.9,
        palette=palette,
        ax=axes[0],
        legend=False
    )
    axes[0].set_title("PCoA - Original", fontsize=14)
    axes[0].set_xlabel(f"PCoA1 [{var_exp_orig.iloc[0]:.1f}%]")
    axes[0].set_ylabel(f"PCoA2 [{var_exp_orig.iloc[1]:.1f}%]")

    sns.scatterplot(
        data=df_corrigida_plot,
        x="PCoA1", y="PCoA2",
        hue="Grupo", style="Estudo",
        s=120, alpha=0.9,
        palette=palette,
        ax=axes[1],
        legend=False
    )
    axes[1].set_title("PCoA - Com correção de batch effects", fontsize=14)
    axes[1].set_xlabel(f"PCoA1 [{var_exp_corr.iloc[0]:.1f}%]")
    axes[1].set_ylabel(f"PCoA2 [{var_exp_corr.iloc[1]:.1f}%]")

    handles_group = [
        mlines.Line2D([], [], color=palette[i], marker='o', linestyle='None', markersize=10,
                      label=g) for i, g in enumerate(df_original_plot['Grupo'].unique())
    ]
    markers = ['o', 's', '^', 'D', 'P', 'X', '*']
    handles_study = [
        mlines.Line2D([], [], color='black', marker=markers[i % len(markers)], linestyle='None', markersize=10,
                      label=s) for i, s in enumerate(df_original_plot['Estudo'].unique())
    ]

    legend_elements = []
    legend_elements.append(mlines.Line2D([], [], color='none', label="Grupo"))
    legend_elements.extend(handles_group)
    legend_elements.append(mlines.Line2D([], [], color='none', label="Estudo"))
    legend_elements.extend(handles_study)

    fig.legend(
        handles=legend_elements,
        loc='center left',
        bbox_to_anchor=(0.92, 0.5),
        frameon=True,
        title=None,
        handletextpad=0.6
    )

    plt.show()