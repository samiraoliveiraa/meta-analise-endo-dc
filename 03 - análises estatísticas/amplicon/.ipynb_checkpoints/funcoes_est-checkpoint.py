import numpy as np
import pandas as pd
from itertools import combinations

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines 
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
import os

from scipy.stats import kruskal, ttest_ind, mannwhitneyu, shapiro
from skbio.diversity import beta_diversity
from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa
import scikit_posthocs as sp
from statsmodels.stats.multitest import multipletests
import statsmodels.stats.multitest as smm
from statannotations.Annotator import Annotator


vaginal = ['endo_hernandes2020', 'endo_jimenez2024', 'endo_perrotta2020', 'meta-analise-vag']
intestinal = ['cfs_giloteaux2016', 'endo_ata2019', 'endo_wei2023', 'fbm_garcia2019', 'fbm_minerbi2019', 'ibs_jacobs2023', 'ibs_vork2021', 'meta-analise-int', 'meta-analise-int-agrupada']

map_labels={"controle": "Controle", 
            "endo": "Endometriose",
            "fbm": "Fibromialgia",
            "cfs": "SFC",
            "ibs": "SII",
            "cpp": "DPC",
            "cpp_endo": "Endometriose + DPC"}


def remove_prefix(taxa_name):
    if taxa_name == 'Outros':
        return taxa_name
    return taxa_name.split("__")[-1]


def taxonomia(artigo, nivel="genero"):

    if artigo in vaginal:
        microbiota = 'vaginal'
    elif artigo in intestinal:
        microbiota = 'intestinal'
    else:
        raise ValueError(f"Artigo '{artigo}' não identificado em vaginal ou intestinal.")

    nivel_map = {
        "filo": "p__",
        "classe": "c__",
        "ordem": "o__",
        "familia": "f__",
        "genero": "g__"
    }
    
    if nivel not in nivel_map:
        raise ValueError(f"Nível '{nivel}' inválido. Use: {list(nivel_map.keys())}")

    prefixo = nivel_map[nivel]

    def extrair_nivel(taxon):
        if pd.isna(taxon):
            return "NA"
        partes = taxon.split(";")
        for p in partes[::-1]:
            if p.startswith(prefixo):
                return p
        return "NA"

    data = pd.read_csv(f'../../02 - processamento/amplicon/{microbiota}/{artigo}/taxonomia_{artigo}.csv')
    data = data.set_index("index")
    data = data.apply(pd.to_numeric, errors="coerce")

    data = data.rename(columns=extrair_nivel)
    data = data.T.groupby(level=0).sum().T
    data = data.dropna(axis=1, how="all")

    data_rel = data.div(data.sum(axis=1), axis=0)
    
    prevalencia = (data_rel >= 0.001).sum(axis=0) / data_rel.shape[0]
    data_rel = data_rel.loc[:, prevalencia >= 0.10]

    total_abundance = data_rel.sum().sort_values(ascending=False)
    top_10_taxa = total_abundance.head(10).index
    data_plot_prep = data_rel[top_10_taxa].copy()

    other_taxa = total_abundance.index.difference(top_10_taxa)
    if not other_taxa.empty:
        data_plot_prep['Outros'] = data_rel[other_taxa].sum(axis=1)

    metadata_path = f'../../02 - processamento/amplicon/{microbiota}/{artigo}/{artigo}_metadata.tsv'
    metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
    metadata = metadata.rename(columns={'group': 'Group'})
    
    metadata["Group"] = metadata["Group"].replace(map_labels)

    data_to_sort = data_plot_prep.join(metadata)
    data_sorted = data_to_sort.sort_values(by=['Group'])
    groups_sorted = data_sorted['Group']
    data_for_plotting = data_sorted.drop(columns=['Group'])
    
    n_taxa = len(data_for_plotting.columns)
    paleta = sns.color_palette("tab20", n_colors=n_taxa)
    
    def remove_prefix(name):
        if pd.isna(name) or name in ["NA", ""]:
            return "Indefinido"
        for prefix in ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]:
            if isinstance(name, str) and name.startswith(prefix):
                return name.replace(prefix, "")
        return name
    
    data_for_plotting_legend = data_for_plotting.rename(columns=remove_prefix)

    fig, ax = plt.subplots(figsize=(13, 7))
    data_for_plotting_legend.plot(kind="bar", stacked=True, width=0.8, ax=ax, color=paleta)
    
    grupos_unicos = groups_sorted.unique()
    group_colors = sns.color_palette("Dark2", n_colors=len(grupos_unicos))
    group_colors = dict(zip(grupos_unicos, group_colors))

    line_y = -0.3
    text_y = -0.35

    start = -0.5
    for grupo in grupos_unicos:
        n = (groups_sorted == grupo).sum()
        end = start + n
        ax.hlines(y=line_y, xmin=start, xmax=end - 0.5,
                  color=group_colors[grupo], linewidth=5,
                  transform=ax.get_xaxis_transform(), clip_on=False)
        ax.text(x=start + (n - 1) / 2, y=text_y, s=grupo,
                ha='center', va='center', transform=ax.get_xaxis_transform())
        start = end

    ax.set_ylabel("Abundância Relativa")
    ax.set_xlabel("")
    ax.set_ylim(0, 1)
    ax.tick_params(axis='x', rotation=90)
    #ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title=f"Taxa ({nivel})")
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3)
    plt.show()
    
    return fig, ax
    
    
def teste_estatistico(df, teste, coluna_grupo, coluna_valor):

    grupos = df[coluna_grupo].unique()
    resultados = []

    if teste == 'kruskal' and len(grupos) > 2:
        dados_por_grupo = [df[df[coluna_grupo] == g][coluna_valor] for g in grupos]

        stat, p = kruskal(*dados_por_grupo)
        resultados.append({
            'teste': 'Kruskal-Wallis',
            'stat': stat,
            'p_value': p
        })

        if p < 0.05:
            dunn = sp.posthoc_dunn(df, val_col=coluna_valor, group_col=coluna_grupo, p_adjust='fdr_bh')

            for g1, g2 in combinations(dunn.columns, 2):
                resultados.append({
                    'teste': 'Dunn',
                    'grupo1': g1,
                    'grupo2': g2,
                    'p_value': dunn.loc[g1, g2],
                    'stat': None
                })

        return resultados

    for grupo1, grupo2 in combinations(grupos, 2):
        dados_grupo1 = df[df[coluna_grupo] == grupo1][coluna_valor]
        dados_grupo2 = df[df[coluna_grupo] == grupo2][coluna_valor]

        if teste == 't':
            stat, p = ttest_ind(dados_grupo1, dados_grupo2, equal_var=False)
        elif teste == 'u':
            stat, p = mannwhitneyu(dados_grupo1, dados_grupo2, alternative='two-sided')
        else:
            raise ValueError(f"Teste '{teste}' não disponível.")

        resultados.append({
            'grupo1': grupo1,
            'grupo2': grupo2,
            'stat': stat,
            'p_value': p
        })

    return resultados


def alfa_diversidade(artigo, teste='u', coluna_grupo='group', coluna_valor='shannon_entropy', grupo_controle='controle', map_labels=map_labels):
    
    if artigo in vaginal:
        microbiota = 'vaginal'
    elif artigo in intestinal:
        microbiota = 'intestinal'
    else:
        raise ValueError(f"Artigo '{artigo}' não identificado em vaginal ou intestinal.")
    
    shannon = pd.read_csv(f'../../02 - processamento/amplicon/{microbiota}/{artigo}/shannon_{artigo}.tsv', sep='\t')
    shannon = shannon.drop(0)
    shannon["shannon_entropy"] = pd.to_numeric(shannon["shannon_entropy"], errors="coerce")

    test_map = {'u': 'Mann-Whitney U', 't': 'Teste t de Welch', 'kruskal': 'Kruskal-Wallis'}
    nome_teste = test_map.get(teste, 'Teste Desconhecido')

    resultados_lista = teste_estatistico(shannon, teste, coluna_grupo, coluna_valor)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    grupos = sorted(shannon[coluna_grupo].unique())
    ordem_grupos = []

    if grupo_controle in grupos:
        ordem_grupos.append(grupo_controle)
    if "endo" in grupos:
        ordem_grupos.append("endo")

    ordem_grupos += [g for g in grupos if g not in ordem_grupos]
    
    sns.boxplot(data=shannon, x=coluna_grupo, y=coluna_valor, order=ordem_grupos,
                hue=coluna_grupo, palette="Set2", legend=False, ax=ax)
    sns.stripplot(data=shannon, x=coluna_grupo, y=coluna_valor, order=ordem_grupos,
                  color="black", alpha=0.6, ax=ax)
    
    ax.set_ylabel("Shannon Index")
    ax.set_xlabel("")

    if map_labels is not None:
        new_labels = [map_labels.get(g, g) for g in ordem_grupos]
        ax.set_xticks(range(len(ordem_grupos)))
        ax.set_xticklabels(new_labels)

    y_max = shannon[coluna_valor].max()
    altura_barra = 0.02 * y_max  
    incremento_altura = 0.08 * y_max 
    y_atual = y_max + altura_barra + 0.05
    
    resultados_significativos = [
        res for res in resultados_lista 
        if 'p_value' in res and res['p_value'] < 0.05 and 'grupo1' in res and 'grupo2' in res
    ]

    for res in resultados_significativos:
        grupo1 = res['grupo1']
        grupo2 = res['grupo2']
        p_val = res['p_value']
        
        x1 = ordem_grupos.index(grupo1)
        x2 = ordem_grupos.index(grupo2)

        ax.plot([x1, x1, x2, x2], 
                [y_atual, y_atual + altura_barra, y_atual + altura_barra, y_atual], 
                lw=0.8, c="black")
        ax.text((x1 + x2) * 0.5, y_atual + altura_barra - 0.02, '*', 
                ha='center', va='bottom', color='black', fontsize=12)
        
        y_atual += incremento_altura 

    ax.set_ylim(0, y_atual + 0.3)
    plt.show()

    df_resultados = pd.DataFrame(resultados_lista)
    print(f'Teste aplicado: {nome_teste}')
    
    return fig, ax, df_resultados


def distribuicao(artigo):
    
    if artigo in vaginal:
        microbiota = 'vaginal'
    elif artigo in intestinal:
        microbiota = 'intestinal'
    else:
        raise ValueError(f"Artigo '{artigo}' não identificado em vaginal ou intestinal.")
    
    df = pd.read_csv(f'../../02 - processamento/amplicon/{microbiota}/{artigo}/shannon_{artigo}.tsv', sep='\t')
    df = df.drop(0)
    df["shannon_entropy"] = pd.to_numeric(df["shannon_entropy"], errors="coerce")

    estatistica = df['shannon_entropy'].describe()

    plt.figure(figsize=(6, 4))
    sns.histplot(df['shannon_entropy'], bins=20, kde=True)
    plt.xlabel("Shannon Entropy")
    plt.ylabel("Frequência")
    plt.title("Distribuição geral dos dados")
    plt.show()
    print(estatistica)
    
    for g, subset in df.groupby('group'):
        plt.figure(figsize=(6,4))
        sns.histplot(subset['shannon_entropy'], bins=15, kde=True)
        plt.xlabel("Shannon Entropy")
        plt.ylabel("Frequência")
        plt.title(f"Distribuição - Grupo {g}")
        plt.show()

        stat, p = shapiro(subset['shannon_entropy'].dropna())
        print(f"\nGrupo: {g}")
        print(f"Tamanho: {len(subset)}")
        print(f"Estatística W = {stat:.4f}, p-valor = {p:.4f}")
        if p > 0.05:
            print("Não rejeitamos H0: Segue uma distribuição normal")
        else:
            print("Rejeitamos H0: Não segue uma distribuição normal")
        print()


def beta_diversidade(artigo, coluna_grupo='group', grupo_controle='controle', map_labels=map_labels, n_permutations=999):
    
    if artigo in vaginal:
        microbiota = 'vaginal'
    elif artigo in intestinal:
        microbiota = 'intestinal'
    else:
        raise ValueError(f"Artigo '{artigo}' não identificado em vaginal ou intestinal.")

    braycurtis = pd.read_csv(f'../../02 - processamento/amplicon/{microbiota}/{artigo}/bc_{artigo}.tsv', sep='\t', index_col=0)
    braycurtis.index = braycurtis.columns
    
    metadata = pd.read_csv(f"../../02 - processamento/amplicon/{microbiota}/{artigo}/{artigo}_metadata.tsv", sep='\t')
    metadata = metadata.rename(columns={"sample-id": "id"})
    metadata["id"] = metadata["id"].astype(str)
    metadata = metadata.set_index("id")

    common_ids = braycurtis.index.intersection(metadata.index)
    braycurtis = braycurtis.loc[common_ids, common_ids]
    metadata = metadata.loc[common_ids]

    contiguous_array = np.ascontiguousarray(braycurtis.values)

    dm = DistanceMatrix(contiguous_array, ids=braycurtis.index)

    permanova_result = permanova(distance_matrix=dm, grouping=metadata[coluna_grupo], permutations=n_permutations)

    ordination = pcoa(dm)

    coords_df = ordination.samples[['PC1', 'PC2']]
    coords_df.index = braycurtis.index.astype(str)
    
    coords_df = coords_df.merge(metadata, left_index=True, right_index=True)

    grupos = sorted(coords_df[coluna_grupo].unique())
    ordem_grupos = [grupo_controle] + [g for g in grupos if g != grupo_controle] if grupo_controle in grupos else grupos

    fig, ax = plt.subplots(figsize=(8, 6))
    
    if map_labels is not None:
        coords_df['grupo_plot'] = coords_df[coluna_grupo].map(lambda x: map_labels.get(x, x))
        ordem_grupos_plot = [map_labels.get(g, g) for g in ordem_grupos]
    else:
        coords_df['grupo_plot'] = coords_df[coluna_grupo]
        ordem_grupos_plot = ordem_grupos

    sns.scatterplot(
        data=coords_df,
        x="PC1",
        y="PC2",
        hue='grupo_plot',
        style='grupo_plot',
        hue_order=ordem_grupos_plot,
        style_order=ordem_grupos_plot,
        s=100,
        palette="Set2",
        ax=ax
    )
    
    pc1_var = ordination.proportion_explained['PC1'] * 100
    pc2_var = ordination.proportion_explained['PC2'] * 100
    ax.set_xlabel(f"PCoA 1 ({pc1_var:.2f}%)")
    ax.set_ylabel(f"PCoA 2 ({pc2_var:.2f}%)")

    leg = ax.legend()
    leg.set_title("Grupo")

    plt.show()
    
    return fig, ax, permanova_result
    

def alfa_multiplot(
    artigos: list,
    teste: str = 'u',
    coluna_grupo: str = 'group',
    coluna_valor: str = 'shannon_entropy',
    grupo_controle: str = 'controle',
    ncols: int = 2,
    salvar: bool = False,
    nome_arquivo: str = "alfa_diversidade"
):
    
    original_params = plt.rcParams.copy()
    try:
        sns.set_theme(context="paper", style="white")
        plt.rcParams.update({
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 1.2,
            "axes.edgecolor": "black",
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
        })

        test_map = {'u': 'Mann-Whitney U', 't': 'Teste t de Welch', 'kruskal': 'Kruskal-Wallis'}
        nome_teste = test_map.get(teste, 'Teste Desconhecido')

        nrows = int(np.ceil(len(artigos) / ncols))
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5*ncols, 4.5*nrows)) 
        axes = axes.flatten() if len(artigos) > 1 else [axes]

        resultados_total = []
        
        mapa_cores_conhecido = {
            grupo_controle: "#5D6D7E",   
            'endo': "#E74C3C",           
            'cpp': "#9B59B6",            
            'cfs': "#ffc759",            
            'fbm': "#E84393",            
            'ibs': "#80ed99",            
            'cpp_endo': "#F39C12",      
            'enxaqueca': "#16A085"       
        }
        
        cor_default = "#AAAAAA" 

        for i, artigo in enumerate(artigos):
            
            map_titles = {
                'endo_hernandes2020': "Endometriose\n(HERNANDES et al., 2020)",
                'endo_jimenez2024': "Endometriose/DPC\n(JIMENEZ et al., 2024)",
                'cfs_giloteaux2016': "SFC\n(GILOTEAUX et al., 2016)",
                'endo_ata2019': "Endometriose\n(ATA et al., 2019)",
                'endo_wei2023': "Endometriose\n(WEI et al., 2023)",
                'fbm_minerbi2019': "Fibromialgia\n(MINERBI et al., 2019)",
                'fbm_garcia2019': "Fibromialgia\n(CLOS-GARCIA et al., 2019)",
                'ibs_jacobs2023': "SII\n(JACOBS et al., 2023)"
            }
            
            if artigo in vaginal:
                microbiota = 'vaginal'
            elif artigo in intestinal:
                microbiota = 'intestinal'
            else:
                raise ValueError(f"Artigo '{artigo}' não identificado em vaginal ou intestinal.")
            
            shannon = pd.read_csv(f'../../02 - processamento/amplicon/{microbiota}/{artigo}/shannon_{artigo}.tsv', sep='\t')
            shannon = shannon.drop(0)
            shannon[coluna_valor] = pd.to_numeric(shannon[coluna_valor], errors='coerce').dropna()

            resultados_lista = teste_estatistico(shannon, teste, coluna_grupo, coluna_valor)
            df_res = pd.DataFrame(resultados_lista)
            df_res["artigo"] = artigo
            resultados_total.append(df_res)

            ax = axes[i]
            grupos = sorted(shannon[coluna_grupo].unique())
            ordem_grupos = []
            if grupo_controle in grupos:
                ordem_grupos.append(grupo_controle)
            if "endo" in grupos:
                ordem_grupos.append("endo")
            ordem_grupos += [g for g in grupos if g not in ordem_grupos]
            paleta_map = {g: mapa_cores_conhecido.get(g, cor_default) for g in ordem_grupos}

            sns.boxplot(
                data=shannon, x=coluna_grupo, y=coluna_valor, order=ordem_grupos,
                hue=coluna_grupo, palette=paleta_map, legend=False, ax=ax,
                fliersize=0, 
                boxprops={'edgecolor':'black', 'linewidth': 1.5},
                whiskerprops={'color':'black', 'linewidth': 1.5},
                capprops={'color':'black', 'linewidth': 1.5},
                medianprops={'color':'black', 'linewidth': 1.5} 
            )
            
            sns.stripplot(
                data=shannon, x=coluna_grupo, y=coluna_valor, order=ordem_grupos,
                hue=coluna_grupo, 
                palette=paleta_map, 
                legend=False, 
                edgecolor="black", 
                linewidth=1.0, 
                alpha=1.0, 
                jitter=0.15, 
                size=5, 
                ax=ax
            )

            titulo = map_titles.get(artigo, artigo)
            ax.set_title(titulo, fontsize=12, fontweight='bold', linespacing=1.6, pad=15)
            ax.set_xlabel("")
            if i % ncols == 0:
                ax.set_ylabel(r"Índice de Shannon (H$'$)", fontsize=11, fontweight='bold')
            else:
                ax.set_ylabel("")

            if map_labels is not None:
                new_labels = [map_labels.get(g, g) for g in ordem_grupos]
                ax.set_xticks(range(len(ordem_grupos)))
                ax.set_xticklabels(new_labels, fontsize=11)

            ax.set_ylim(bottom=0) 
            y_min_ax, y_max_ax = ax.get_ylim()
            y_range = y_max_ax - y_min_ax 
            
            data_max = shannon[coluna_valor].max()
            
            altura_barra = 0.02 * y_range
            incremento_altura = 0.1 * y_range 
            y_atual = data_max + 0.05 * y_range 

            resultados_significativos = [
                res for res in resultados_lista
                if 'p_value' in res and res['p_value'] < 0.05 and 'grupo1' in res and 'grupo2' in res
            ]
            
            num_barras = 0
            for res in resultados_significativos:
                grupo1, grupo2 = res['grupo1'], res['grupo2']
                if grupo1 not in ordem_grupos or grupo2 not in ordem_grupos:
                    continue 
                
                x1 = ordem_grupos.index(grupo1)
                x2 = ordem_grupos.index(grupo2)
                
                ax.plot([x1, x1, x2, x2],
                        [y_atual, y_atual + altura_barra, y_atual + altura_barra, y_atual],
                        lw=1.5, c="black") 
               
                p_valor = res.get('p_value', None)
                if p_valor is not None:
                    p_formatado = f"p = {p_valor:.3f}" if p_valor >= 0.001 else "p < 0.001"
                else:
                    p_formatado = ""

                ax.text(
                    (x1 + x2) * 0.5, 
                    y_atual + altura_barra, 
                    '*',
                    ha='center', 
                    va='bottom', 
                    color='black', 
                    fontsize=14,
                    fontweight='bold'
                )

                y_atual += incremento_altura 
                num_barras += 1

            if num_barras > 0:
                ax.set_ylim(bottom=0, top=y_atual + 0.1 * y_range)
            
            sns.despine(ax=ax)
            
        n = len(artigos)
        for j in range(i + 1, len(axes)):
            fig.delaxes(axes[j])

        fig.tight_layout(rect=[0, 0, 1, 0.95], h_pad=2.5, w_pad=1.2)
        
        if (n % ncols != 0) and (n > ncols): 
            ax_last = axes[n-1]
            pos = ax_last.get_position()
            new_x0 = 0.5 - pos.width / 2
            ax_last.set_position([
                new_x0, pos.y0, pos.width, pos.height
            ])

        df_resultados_total = pd.concat(resultados_total, ignore_index=True)
        if salvar:
            fig.savefig(f"{nome_arquivo}.png", dpi=300, bbox_inches='tight', transparent=True)
            print(f"✅ Figura salva como '{nome_arquivo}.png'")
        return fig, axes, df_resultados_total
    
    finally:
        plt.rcParams.update(original_params)
        

def beta_multiplot(
    artigos: list,
    coluna_grupo: str = 'group',
    grupo_controle: str = 'controle',
    n_permutations: int = 999,
    ncols: int = 2,
    salvar: bool = False,
    nome_arquivo: str = "beta_diversidade"
):

    original_params = plt.rcParams.copy()
    try:
        sns.set_theme(context="paper", style="white")
        plt.rcParams.update({
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 1.2,
            "axes.edgecolor": "black",
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
        })

        mapa_cores_conhecido = {
            grupo_controle: "#5D6D7E",   
            'endo': "#E74C3C",           
            'cpp': "#9B59B6",            
            'cfs': "#ffc759",            
            'fbm': "#E84393",            
            'ibs': "#80ed99",            
            'cpp_endo': "#F39C12",    
            'enxaqueca': "#16A085"     
        }
        cor_default = "#999999"

        n = len(artigos)
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(
            nrows=nrows, ncols=ncols,
            figsize=(5.8 * ncols, 5.2 * nrows)
        )
        axes = axes.flatten() if n > 1 else [axes]

        permanova_results_list = []

        for i, artigo in enumerate(artigos):
            ax = axes[i]
            map_titles = {
                'endo_hernandes2020': "Endometriose\n(HERNANDES et al., 2020)",
                'endo_jimenez2024': "Endometriose/DPC\n(JIMENEZ et al., 2024)",
                'cfs_giloteaux2016': "SFC\n(GILOTEAUX et al., 2016)",
                'endo_ata2019': "Endometriose\n(ATA et al., 2019)",
                'endo_wei2023': "Endometriose\n(WEI et al., 2023)",
                'fbm_minerbi2019': "Fibromialgia\n(MINERBI et al., 2019)",
                'fbm_garcia2019': "Fibromialgia\n(CLOS-GARCIA et al., 2019)",
                'ibs_jacobs2023': "SII\n(JACOBS et al., 2023)"
            }

            if artigo in vaginal:
                microbiota = 'vaginal'
            elif artigo in intestinal:
                microbiota = 'intestinal'
            else:
                raise ValueError(f"Artigo '{artigo}' não identificado em vaginal ou intestinal.")
                
            braycurtis = pd.read_csv(f'../../02 - processamento/amplicon/{microbiota}/{artigo}/bc_{artigo}.tsv', sep='\t', index_col=0)
            braycurtis.index = braycurtis.columns
            metadata = pd.read_csv(f"{microbiota}/{artigo}/{artigo}_metadata.tsv", sep='\t')
            metadata = metadata.rename(columns={"sample-id": "id"}).set_index("id")

            common_ids = braycurtis.index.intersection(metadata.index)
            braycurtis = braycurtis.loc[common_ids, common_ids]
            metadata = metadata.loc[common_ids]

            dm = DistanceMatrix(np.ascontiguousarray(braycurtis.values), ids=braycurtis.index)
            permanova_result = permanova(distance_matrix=dm, grouping=metadata[coluna_grupo], permutations=n_permutations)

            permanova_results_list.append(pd.Series({
                'artigo': artigo,
                'F': permanova_result['test statistic'],
                'p': permanova_result['p-value']
            }))

            ordination = pcoa(dm)
            coords_df = ordination.samples[['PC1', 'PC2']]
            coords_df = coords_df.merge(metadata, left_index=True, right_index=True)

            grupos = coords_df[coluna_grupo].unique()
            paleta_local = {g: mapa_cores_conhecido.get(g, cor_default) for g in grupos}

            sns.scatterplot(
                data=coords_df,
                x="PC1",
                y="PC2",
                hue=coluna_grupo,
                palette=paleta_local,
                s=65,
                alpha=1,
                edgecolor="black",
                linewidth=0.8,
                ax=ax
            )

            ax.axhline(0, ls='--', color='gray', lw=0.8, zorder=0)
            ax.axvline(0, ls='--', color='gray', lw=0.8, zorder=0)

            pc1_var = ordination.proportion_explained['PC1'] * 100
            pc2_var = ordination.proportion_explained['PC2'] * 100

            ax.set_xlabel(f"PCoA 1 ({pc1_var:.2f}%)", fontsize=11, fontweight='bold')
            ax.set_ylabel(f"PCoA 2 ({pc2_var:.2f}%)", fontsize=11, fontweight='bold')
            ax.set_title(map_titles.get(artigo, artigo), fontweight='bold', pad=12)

            f_val = permanova_result['test statistic']
            p_val = permanova_result['p-value']
            texto = f"PERMANOVA\nF = {f_val:.3f}\np = {p_val:.4f}"
            ax.text(
                0.02, 0.98, texto,
                transform=ax.transAxes,
                fontsize=9.5,
                fontweight='bold',
                va='top', ha='left',
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.35')
            )

            leg = ax.legend(
                title=None,
                loc='upper right',
                frameon=True,
                facecolor='white',
                edgecolor='black',
                fontsize=9
            )
            leg.get_frame().set_linewidth(0.8)
            
            for text in leg.get_texts():
                label_original = text.get_text()
                text.set_text(map_labels.get(label_original, label_original))

            sns.despine(ax=ax)

        for j in range(n, len(axes)):
            fig.delaxes(axes[j])

        fig.tight_layout(h_pad=2.5, w_pad=1.4)

        if salvar:
            fig.savefig(f"{nome_arquivo}.png", dpi=600, bbox_inches='tight')
            print(f"✅ Figura salva como '{nome_arquivo}.png'")

        return fig, axes, pd.DataFrame(permanova_results_list).set_index('artigo')

    finally:
        plt.rcParams.update(original_params)
        
        
def taxonomia_multiplot(
    estudos=None,
    grupos_doencas=None,
    niveis=('filo', 'genero'),
    top_n=10,
    width=18,
    height=6,
    salvar_arquivo=None
):

    doenca_map = {
        'cfs_giloteaux2016': 'SFC',
        'endo_ata2019': 'Endometriose',
        'endo_wei2023': 'Endometriose',
        'fbm_garcia2019': 'Fibromialgia',
        'fbm_minerbi2019': 'Fibromialgia',
        'ibs_jacobs2023': 'SII',
        'endo_hernandes': 'Endometriose',
    }

    nivel_map = {"filo": "p__", "genero": "g__"}
    nomes_niveis = {"filo": "Filo", "genero": "Gênero"}
    ordem_doencas = ["Controle", "SFC", "Endometriose", "Fibromialgia", "SII"]

    cores = [
        '#6D597A', '#355C7D', '#E84A5F', '#FF847C', '#6A0572',
        '#1A936F', '#114B5F', '#45062E', '#3A015C', '#1B512D',
        '#5D2A42', '#3C1642', '#1E3F20', '#8F2D56', '#218380'
    ]
    cor_outros = "#d9d9d9"

    def configurar_estilo():
        plt.rcParams.update({
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
        })

    def extrair_nivel(taxon, nivel):
        prefixo = nivel_map[nivel]
        if pd.isna(taxon):
            return "NA"
        for p in taxon.split(";"):
            if p.strip().startswith(prefixo):
                return p.replace(prefixo, "").strip()
        return "NA"

    def agrupar_por_nivel(df, nivel):
        grupos = {}
        for col in df.columns:
            nome = extrair_nivel(col, nivel)
            grupos[nome] = grupos.get(nome, 0) + df[col]
        return pd.DataFrame(grupos)

    if estudos is None:
        estudos = ['cfs_giloteaux2016', 'endo_ata2019', 'endo_wei2023',
                   'fbm_garcia2019', 'fbm_minerbi2019', 'ibs_jacobs2023']

    dados = {}
    for art in estudos:
        caminho_data = f'../../02 - processamento/amplicon/intestinal/{art}/taxonomia_{art}.csv'
        caminho_meta = f'../../02 - processamento/amplicon/intestinal/{art}/{art}_metadata.tsv'

        if not os.path.exists(caminho_data) or not os.path.exists(caminho_meta):
            continue

        df = pd.read_csv(caminho_data).set_index("index")
        df = df.apply(pd.to_numeric, errors="coerce")

        meta = pd.read_csv(caminho_meta, sep="\t", index_col=0)
        meta["Group"] = meta["group"].replace(map_labels)
        meta["Doenca"] = doenca_map.get(art, "Indefinido")

        dados[art] = {"relativa": df, "meta": meta}

    if grupos_doencas is None:
        grupos_doencas = {
            "Controle": [
                "cfs_giloteaux2016", "endo_ata2019", "endo_wei2023",
                "fbm_garcia2019", "fbm_minerbi2019", "ibs_jacobs2023"
            ],
            "SFC": ["cfs_giloteaux2016"],
            "Endometriose": ["endo_ata2019", "endo_wei2023"],
            "Fibromialgia": ["fbm_garcia2019", "fbm_minerbi2019"],
            "SII": ["ibs_jacobs2023"]
        }

    configurar_estilo()
    fig, axes = plt.subplots(1, len(niveis), figsize=(width, height))

    for ax, nivel in zip(axes, niveis):
        dados_por_doenca = {}
        n_amostras = {}

        for doenca, lista_estudos in grupos_doencas.items():
            combinados = []

            for est in lista_estudos:
                if est not in dados:
                    continue

                df = dados[est]["relativa"]
                meta = dados[est]["meta"]

                if doenca == "Controle":
                    amostras = meta[meta["Group"] == "Controle"].index
                else:
                    amostras = meta[meta["Group"] != "Controle"].index

                amostras = amostras.intersection(df.index)
                if len(amostras) == 0:
                    continue

                df_sel = df.loc[amostras]
                df_agr = agrupar_por_nivel(df_sel, nivel)
                df_rel = df_agr.div(df_agr.sum(axis=1), axis=0) * 100
                combinados.append(df_rel)

            if combinados:
                df_final = pd.concat(combinados)
                dados_por_doenca[doenca] = df_final.mean()
                n_amostras[doenca] = len(df_final)

        if not dados_por_doenca:
            ax.text(0.5, 0.5, "Sem dados", ha="center")
            continue

        df_medias = pd.DataFrame(dados_por_doenca)

        top_taxa = df_medias.mean(axis=1).nlargest(top_n).index
        df_plot = df_medias.loc[top_taxa].copy()
        df_plot.loc["Outros"] = df_medias.loc[~df_medias.index.isin(top_taxa)].sum()

        ordem_taxons = df_plot.sum(axis=1).sort_values(ascending=False).index
        colunas = [c for c in ordem_doencas if c in df_plot.columns]
        df_plot = df_plot.loc[ordem_taxons, colunas].fillna(0)

        bottom = np.zeros(len(colunas))
        patches = []

        for i, taxon in enumerate(df_plot.index):
            vals = df_plot.loc[taxon]
            cor = cor_outros if taxon == "Outros" else cores[i % len(cores)]

            ax.bar(colunas, vals, bottom=bottom, color=cor, edgecolor="black")
            bottom += vals

            patches.append(Patch(facecolor=cor, edgecolor="black", label=taxon))

        ax.set_ylim(0, 100)
        ax.set_ylabel("Abundância relativa média (%)")
        ax.set_xticklabels([f"{c}\n(n={n_amostras.get(c,0)})" for c in colunas])

        ax.legend(handles=patches, bbox_to_anchor=(1.05,1), loc="upper left")

        ax.set_title(nomes_niveis[nivel])

    plt.tight_layout()

    if salvar_arquivo:
        plt.savefig(f"{salvar_arquivo}.png", bbox_inches="tight")
        print(f"Figura salva como {salvar_arquivo}.png")

    plt.show()
    return fig, axes
