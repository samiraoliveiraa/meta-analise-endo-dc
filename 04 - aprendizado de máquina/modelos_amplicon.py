import numpy as np
import pandas as pd
import joblib
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, confusion_matrix
from sklearn.linear_model import LogisticRegression
from interpret.glassbox import ExplainableBoostingClassifier
from catboost import CatBoostClassifier
from sklearn.ensemble import RandomForestClassifier
from skbio.stats.composition import clr
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import optuna
import os

# Constantes
SEMENTE_ALEATORIA = 24
TAMANHO_TESTE = 0.15
PSEUDOCOUNT = 1
N_TRIALS = 500

os.makedirs("modelos-amplicon", exist_ok=True)

def k_valor(y, min_per_fold=3, k_min=3, k_max=10):
    counts = np.bincount(y)
    n_min = counts.min()
    k = max(k_min, min(k_max, n_min // min_per_fold))
    return k


def otimizar_modelo(modelo_nome, X, y):
    N_SPLITS = k_valor(y)

    def objetivo(trial):
        if modelo_nome == "LogReg":
            C = trial.suggest_float("C", 1e-3, 1.0, log=True)
            penalty = trial.suggest_categorical("penalty", ["l1", "l2"])
            modelo = LogisticRegression(
                solver='saga',
                max_iter=5000,
                C=C,
                penalty=penalty,
                random_state=SEMENTE_ALEATORIA
            )

        elif modelo_nome == "RandomForest":
            n_estimators = trial.suggest_int("n_estimators", 200, 600)
            max_depth = trial.suggest_int("max_depth", 3, 10)
            min_samples_split = trial.suggest_int("min_samples_split", 4, 12)
            min_samples_leaf = trial.suggest_int("min_samples_leaf", 2, 6)
            modelo = RandomForestClassifier(
                n_estimators=n_estimators,
                max_depth=max_depth,
                min_samples_split=min_samples_split,
                min_samples_leaf=min_samples_leaf,
                random_state=SEMENTE_ALEATORIA
            )

        elif modelo_nome == "EBM":
            max_bins = trial.suggest_int("max_bins", 64, 256)
            max_rounds = trial.suggest_int("max_rounds", 50, 200)
            learning_rate = trial.suggest_float("learning_rate", 0.001, 0.1)
            max_leaves = trial.suggest_int("max_leaves", 2, 5)
            modelo = ExplainableBoostingClassifier(
                max_rounds=max_rounds,
                max_leaves=max_leaves,
                learning_rate=learning_rate,
                max_bins=max_bins,
                random_state=SEMENTE_ALEATORIA,
                interactions=0
            )

        elif modelo_nome == "CatBoost":
            iterations = trial.suggest_int("iterations", 100, 300)
            depth = trial.suggest_int("depth", 3, 8)
            learning_rate = trial.suggest_float("learning_rate", 0.001, 0.1)
            l2_leaf_reg = trial.suggest_float("l2_leaf_reg", 1, 10)
            min_data_in_leaf = trial.suggest_int("min_data_in_leaf", 5, 30)
            modelo = CatBoostClassifier(
                iterations=iterations,
                depth=depth,
                min_data_in_leaf=min_data_in_leaf,
                learning_rate=learning_rate,
                l2_leaf_reg=l2_leaf_reg,
                verbose=0,
                random_state=SEMENTE_ALEATORIA
            )
        else:
            raise ValueError(f"Modelo '{modelo_nome}' não suportado.")
        
        cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEMENTE_ALEATORIA)
        scores = cross_val_score(modelo, X, y, cv=cv, scoring='balanced_accuracy', n_jobs=-1)
        return np.mean(scores)

    study = optuna.create_study(direction="maximize")
    study.optimize(objetivo, n_trials=N_TRIALS)

    print(f"\nMelhores parâmetros para {modelo_nome}: {study.best_params}")
    print(f"Número de folds: {N_SPLITS}")
    print(f"Melhor acurácia balanceada média (CV): {study.best_value:.3f}")
    return study.best_params


def avaliar_modelos(nome, dados):
    print(f"\n=== {nome} ===")

    y = (dados["Condição"] != "controle").astype(int).values
    dados = dados.drop(columns=["Condição"])

    print(f'Número de táxons inicial: {len(dados.columns)}')

    dados_rel = dados.div(dados.sum(axis=1), axis=0)
    keep_taxa = (dados_rel >= 0.0001).sum(axis=0) >= (0.10 * dados_rel.shape[0])
    dados = dados.loc[:, keep_taxa]

    print(f'Número de táxons pós-filtro: {len(dados.columns)}')

    X = dados.values

    X_treino, X_teste, y_treino, y_teste = train_test_split(
        X, y, test_size=TAMANHO_TESTE, random_state=SEMENTE_ALEATORIA, stratify=y
    )

    X_treino_clr = clr(X_treino + PSEUDOCOUNT)
    X_teste_clr = clr(X_teste + PSEUDOCOUNT)

    grupos = np.unique(y_treino)
    g1 = X_treino_clr[y_treino == grupos[0]]
    g2 = X_treino_clr[y_treino == grupos[1]]

    pvals = []

    for i in range(X_treino_clr.shape[1]):
        try:
            stat, p = mannwhitneyu(g1[:, i], g2[:, i], alternative="two-sided")
        except ValueError:
            p = 1.0
        pvals.append(p)

    pvals = np.array(pvals)

    _, pvals_corr, _, _ = multipletests(pvals, method="fdr_bh")

    n_features_desejado = int(0.3 * X_treino_clr.shape[0])
    idx_top = np.argsort(pvals_corr)[:n_features_desejado]

    X_treino_sel = X_treino_clr[:, idx_top]
    X_teste_sel  = X_teste_clr[:, idx_top]

    dados = dados.iloc[:, idx_top]

    print(f'Número de táxons pós-seleção: {len(dados.columns)}')

    normalizador = StandardScaler()
    X_treino_norm = normalizador.fit_transform(X_treino_sel)
    X_teste_norm = normalizador.transform(X_teste_sel)

    pesos = compute_sample_weight(class_weight='balanced', y=y_treino)

    modelos = ["LogReg", "RandomForest", "EBM", "CatBoost"]

    for nome_modelo in modelos:
        print(f"\n--- Otimizando {nome_modelo} ---")
        best_params = otimizar_modelo(nome_modelo, X_treino_norm, y_treino)

        if nome_modelo == "LogReg":
            modelo = LogisticRegression(max_iter=5000, solver='saga', class_weight='balanced', random_state=SEMENTE_ALEATORIA, **best_params)
        elif nome_modelo == "RandomForest":
            modelo = RandomForestClassifier(class_weight='balanced', random_state=SEMENTE_ALEATORIA, **best_params)
        elif nome_modelo == "EBM":
            modelo = ExplainableBoostingClassifier(random_state=SEMENTE_ALEATORIA, interactions=0, **best_params)
        elif nome_modelo == "CatBoost":
            modelo = CatBoostClassifier(verbose=0, random_state=SEMENTE_ALEATORIA, **best_params)

        modelo.fit(X_treino_norm, y_treino, sample_weight=pesos if nome_modelo in ["LogReg", "EBM", "CatBoost"] else None)
        
        y_pred = modelo.predict(X_teste_norm)
        y_proba = modelo.predict_proba(X_teste_norm)[:, 1]
        
        auc = roc_auc_score(y_teste, y_proba)
        bacc = balanced_accuracy_score(y_teste, y_pred)
        cm = confusion_matrix(y_teste, y_pred)
        
        print(f"ROC AUC no conjunto de teste: {auc:.3f}")
        print(f"Acurácia balanceada no conjunto de teste: {bacc:.3f}")
        print(f"Matriz de confusão:\n{cm}")

        nome_limpo = nome.replace(" ", "_").replace(".", "")
        caminho_modelo = f"modelos-amplicon/{nome_limpo}_{nome_modelo}.joblib"
        joblib.dump(modelo, caminho_modelo)
        print(f"Modelo salvo em: {caminho_modelo}")

    print("-" * 60)


data = pd.read_csv('taxonomia_amplicon_corrigida.csv').rename({'Unnamed: 0': 'index'}, axis=1).set_index('index')
metadata = pd.read_csv('../02 - processamento/amplicon/intestinal/meta-analise-int/meta-analise-int_metadata.tsv', sep='\t').set_index("sample-id")
data["Condição"] = metadata.loc[data.index, "group"]

pares_separados = [
    ("Controle vs. FBM", data.loc[~data["Condição"].isin(['endo', 'ibs', 'cfs'])]),
    ("Controle vs. CFS", data.loc[~data["Condição"].isin(['endo', 'ibs', 'fbm'])]),
    ("Controle vs. IBS", data.loc[~data["Condição"].isin(['endo', 'fbm', 'cfs'])]),
    ("Controle vs. Endometriose", data.loc[~data["Condição"].isin(['ibs', 'fbm', 'cfs'])]),
]

print("\n=== Análises por doença ===")
for nome, dados in pares_separados:
    avaliar_modelos(nome, dados)
