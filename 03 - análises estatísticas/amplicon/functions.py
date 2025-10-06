import numpy as np
import pandas as pd
from itertools import combinations

import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import kruskal, ttest_ind, mannwhitneyu, shapiro
from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa
import scikit_posthocs as sp
from statsmodels.stats.multitest import multipletests


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
    

def taxonomia_extendida(artigo):
    
    if artigo in vaginal:
        microbiota == 'vaginal'
    elif artigo in intestinal:
        microbiota == 'intestinal'
    else:
        print('Artigo não identificado.')
    
    data = pd.read_csv(f'{microbiota}/{artigo}/bar-plot_{artigo}.csv')
    data = data.set_index("index")
    data = data.apply(pd.to_numeric, errors="coerce")
    data = data.rename(columns=extrair_ultimo_nome)
    data = data.dropna(axis=1, how="all")

    data_rel = data.div(data.sum(axis=1), axis=0)

    metadata_path = f'{microbiota}/{artigo}/{artigo}_metadata.tsv'
    metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
    metadata.columns = ['Group']

    data_to_sort = data_rel.join(metadata)

    data_sorted = data_to_sort.sort_values(by=['Group', 'index'])

    groups_sorted = data_sorted['Group']
    data_for_plotting = data_sorted.drop(columns=['Group'])
    
    fig, ax = plt.subplots(figsize=(18, 7))
    data_for_plotting_legend = data_for_plotting.rename(columns=remove_prefix)
    data_for_plotting_legend.plot(kind="bar", stacked=True, width=0.8, ax=ax)

    group_colors = {'controle': 'darkcyan', 'endo': 'goldenrod'}

    n_samples = len(data_for_plotting)
    n_control = (groups_sorted == 'controle').sum()

    line_y = -0.3
    text_y = -0.35

    ax.hlines(y=line_y, xmin=-0.5, xmax=n_control - 0.5, 
              color=group_colors['controle'], linewidth=5, 
              transform=ax.get_xaxis_transform(), clip_on=False)
    ax.text(x=(n_control - 1) / 2, y=text_y, s='Controle', 
            ha='center', va='center', transform=ax.get_xaxis_transform())

    if n_control < n_samples:
        ax.hlines(y=line_y, xmin=n_control - 0.5, xmax=n_samples - 0.5, 
                  color=group_colors['endo'], linewidth=5, 
                  transform=ax.get_xaxis_transform(), clip_on=False)
        ax.text(x=n_control + (n_samples - n_control - 1) / 2, y=text_y, s='Endometriose', 
                ha='center', va='center', transform=ax.get_xaxis_transform())

    ax.set_title(f"Composição Taxonômica - {artigo}")
    ax.set_ylabel("Proporção Relativa das Leituras")
    ax.set_xlabel("")
    ax.set_ylim(0, 1)
    ax.tick_params(axis='x', rotation=90)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Taxa (Gênero)")

    plt.subplots_adjust(bottom=0.25)

    plt.show()
    
    return fig, ax



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
        
    def extrair_nivel(taxon):
        if pd.isna(taxon):
            return "NA"
        partes = taxon.split(";")
        for p in partes[::-1]:
            if p.startswith(prefixo):
                return p
        return "NA"
    
    prefixo = nivel_map[nivel]

    data = pd.read_csv(f'{microbiota}/{artigo}/taxonomia_{artigo}.csv')
    data = data.set_index("index")
    data = data.apply(pd.to_numeric, errors="coerce")

    data = data.rename(columns=extrair_nivel)
    data = data.T.groupby(level=0).sum().T
    data = data.dropna(axis=1, how="all")

    data_rel = data.div(data.sum(axis=1), axis=0)

    total_abundance = data_rel.sum().sort_values(ascending=False)
    top_10_taxa = total_abundance.head(10).index
    data_plot_prep = data_rel[top_10_taxa].copy()

    other_taxa = total_abundance.index.difference(top_10_taxa)
    if not other_taxa.empty:
        data_plot_prep['Outros'] = data_rel[other_taxa].sum(axis=1)

    metadata_path = f'{microbiota}/{artigo}/{artigo}_metadata.tsv'
    metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
    metadata.columns = ['Group']

    data_to_sort = data_plot_prep.join(metadata)
    data_sorted = data_to_sort.sort_values(by=['Group'])
    groups_sorted = data_sorted['Group']
    data_for_plotting = data_sorted.drop(columns=['Group'])
    
    n_taxa = len(data_for_plotting.columns)
    paleta = sns.color_palette("tab20", n_colors=n_taxa)

    fig, ax = plt.subplots(figsize=(18, 7))
    data_for_plotting_legend = data_for_plotting.rename(columns=remove_prefix)
    data_for_plotting_legend.plot(kind="bar", stacked=True, width=0.8, ax=ax)
    
    group_colors = {'controle': '#262626', 'endo': '#FF7F11'}
    n_samples = len(data_for_plotting)
    n_control = (groups_sorted == 'controle').sum()

    line_y = -0.3
    text_y = -0.35

    ax.hlines(y=line_y, xmin=-0.5, xmax=n_control - 0.5,
              color=group_colors['controle'], linewidth=5,
              transform=ax.get_xaxis_transform(), clip_on=False)
    ax.text(x=(n_control - 1) / 2, y=text_y, s='Controle',
            ha='center', va='center', transform=ax.get_xaxis_transform())

    if n_control < n_samples:
        ax.hlines(y=line_y, xmin=n_control - 0.5, xmax=n_samples - 0.5,
                  color=group_colors['endo'], linewidth=5,
                  transform=ax.get_xaxis_transform(), clip_on=False)
        ax.text(x=n_control + (n_samples - n_control - 1) / 2, y=text_y, s='Endometriose',
                ha='center', va='center', transform=ax.get_xaxis_transform())

    #ax.set_title(f"Composição Taxonômica ({nivel}) - {artigo}")
    ax.set_ylabel("Abundância Relativa")
    ax.set_xlabel("")
    ax.set_ylim(0, 1)
    ax.tick_params(axis='x', rotation=90)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title=f"Taxa ({nivel})")

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
            dunn = sp.posthoc_dunn(df, val_col=coluna_valor, group_col=coluna_grupo, p_adjust='bonferroni')

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
    
    shannon = pd.read_csv(f'{microbiota}/{artigo}/shannon_{artigo}.tsv', sep='\t')
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
    
    df = pd.read_csv(f'{microbiota}/{artigo}/shannon_{artigo}.tsv', sep='\t')
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

    braycurtis = pd.read_csv(f'{microbiota}/{artigo}/bc_{artigo}.tsv', sep='\t', index_col=0)
    braycurtis.index = braycurtis.columns
    
    metadata = pd.read_csv(f"{microbiota}/{artigo}/{artigo}_metadata.tsv", sep='\t')
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

    
def vias_metabolicas(artigo):
    
    if artigo in vaginal:
        microbiota = 'vaginal'
    elif artigo in intestinal:
        microbiota = 'intestinal'
    else:
        raise ValueError(f"Artigo '{artigo}' não identificado em vaginal ou intestinal.")

    
    metadata_path = f'{microbiota}/{artigo}/{artigo}_metadata.tsv'
    metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
    
    df = pd.read_csv(f'{microbiota}/{artigo}/pathway_{artigo}.tsv', sep='\t', skiprows=1, index_col=0)

    df_log = np.log1p(df)

    g = sns.clustermap(df_log,
                       z_score=0,
                       cmap='rocket',
                       figsize=(12, 18),
                       xticklabels=True,
                       cbar_kws={'label': 'Abundância Relativa (Z-score)'})

    g.fig.suptitle(f'Mapa de Calor das Vias Metabólicas - {artigo}', fontsize=14, y=1.02)

    reordered_samples = df_log.columns[g.dendrogram_col.reordered_ind]

    new_labels = [f"{metadata.loc[sample, coluna_grupo]} | {sample}"
                  for sample in reordered_samples]

    g.ax_heatmap.set_xticklabels(new_labels, rotation=90, ha='center')
    g.ax_heatmap.set_ylabel("Via Metabólica") 

    plt.show()
    
    return g

def clr_transform(df):
    
    df = df + 1 
    log_df = np.log(df)
    gm = log_df.mean(axis=0)
    
    return log_df.sub(gm, axis=1)


def abundancia_funcional(artigo, coluna_grupo='group'):
    
    if artigo in vaginal:
        microbiota = 'vaginal'
    elif artigo in intestinal:
        microbiota = 'intestinal'
    else:
        raise ValueError(f"Artigo '{artigo}' não identificado em vaginal ou intestinal.")

    table = pd.read_csv(f"{microbiota}/{artigo}/taxonomia_{artigo}.csv", index_col=0)
    meta = pd.read_csv(f"{microbiota}/{artigo}/{artigo}_metadata.tsv", sep='\t', index_col=0)
    table = table.drop(columns=['group'])
    table = table.T

    clr_table = clr_transform(table)
    
    grupos = meta[coluna_grupo].unique()

    results = []
    for genus in clr_table.index:
        group1 = clr_table.loc[genus, meta[meta[coluna_grupo]=="endo"].index]
        group2 = clr_table.loc[genus, meta[meta[coluna_grupo]=="controle"].index]
        stat, p = mannwhitneyu(group1, group2, alternative="two-sided")
        results.append([genus, p])

    res_df = pd.DataFrame(results, columns=["Gênero", "p-valor"])
    res_df["FDR"] = res_df["p-valor"] * len(res_df) 
    res_df = res_df.sort_values("p-valor")

    return res_df

def analise_diferencial_abundancia(artigo, coluna_grupo='group', p_adjust_method='fdr_bh', feature_level="Táxon", nivel='genero'):
    
    if artigo in vaginal:
        microbiota = 'vaginal'
    elif artigo in intestinal:
        microbiota = 'intestinal'
    else:
        raise ValueError(f"Artigo '{artigo}' não identificado em listas de microbiota.")
        
    table = pd.read_csv(f"{microbiota}/{artigo}/taxonomia_{artigo}.csv", index_col=0)
    meta = pd.read_csv(f"{microbiota}/{artigo}/{artigo}_metadata.tsv", sep='\t', index_col=0)
    
    nivel_map = {
        "filo": "p__", "classe": "c__", "ordem": "o__",
        "familia": "f__", "genero": "g__"
    }
    
    prefixo = nivel_map.get(nivel)
    if not prefixo:
        raise ValueError(f"Nível taxonômico '{nivel}' não reconhecido.")

    def extrair_nivel(taxon):
        if pd.isna(taxon): return "NA"
        partes = str(taxon).split(";")
        for p in partes[::-1]:
            p = p.strip() # Limpa espaços em branco
            if p.startswith(prefixo):
                return p
        return "NA"
    
    table.columns = table.columns.map(extrair_nivel)

    table = table.groupby(by=table.columns, axis=1).sum()

    if 'NA' in table.columns:
        table = table.drop(columns=['NA'])
    if 'group' in table.columns:
        table = table.drop(columns=['group'])

    table = table.T
    
    amostras_comuns = table.columns.intersection(meta.index)
    table = table[amostras_comuns]
    meta = meta.loc[amostras_comuns]

    clr_table = clr_transform(table)
    
    grupos = meta[coluna_grupo].unique()

    if len(grupos) > 2:
        kruskal_results = []
        for feature in clr_table.index:
            dados_por_grupo = [clr_table.loc[feature, meta[meta[coluna_grupo] == g].index].dropna() for g in grupos]
            if any(len(d) < 1 for d in dados_por_grupo): continue
            stat, p = kruskal(*dados_por_grupo)
            kruskal_results.append({feature_level: feature, 'statistic_kruskal': stat, 'p_kruskal': p})
        
        res_kruskal_df = pd.DataFrame(kruskal_results)
        
        clr_table_transposed = clr_table.T
        clr_table_transposed[coluna_grupo] = meta[coluna_grupo]
        
        p_values_dunn = sp.posthoc_dunn(clr_table_transposed, val_col=clr_table.index.tolist(), group_col=coluna_grupo, p_adjust=p_adjust_method)
        
        p_values_dunn.columns.name = feature_level
        p_values_dunn.index.name = "grupo1"
        p_values_dunn = p_values_dunn.stack().reset_index(name='p_valor_ajustado')
        p_values_dunn.rename(columns={'level_1': 'grupo2'}, inplace=True)
        
        return res_kruskal_df, p_values_dunn

    elif len(grupos) == 2:
        results = []
        grupo1_nome, grupo2_nome = grupos[0], grupos[1]
        
        for feature in clr_table.index:
            grupo1_data = clr_table.loc[feature, meta[meta[coluna_grupo] == grupo1_nome].index]
            grupo2_data = clr_table.loc[feature, meta[meta[coluna_grupo] == grupo2_nome].index]
            
            if grupo1_data.empty or grupo2_data.empty: continue
            
            stat, p = mannwhitneyu(grupo1_data, grupo2_data, alternative="two-sided")
            results.append([feature, stat, p])

        res_df = pd.DataFrame(results, columns=[feature_level, "statistic", "p-valor"])
        
        if not res_df.empty:
            res_df.dropna(subset=['p-valor'], inplace=True)
            reject, pvals_corrected, _, _ = multipletests(res_df["p-valor"], alpha=0.05, method=p_adjust_method)
            res_df["p-valor_ajustado"] = pvals_corrected
            res_df["significativo"] = reject
        
        res_df = res_df.sort_values("p-valor")
        
        return res_df