#!/usr/bin/env python3
"""
GSE1456 — 06 Machine Learning: ER+ vs ER- classification
Setup-1: biological gene list (DGE + WGCNA hubs)
Setup-2: VarianceThreshold + SelectKBest (k=100)
"""

from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    auc,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

warnings.filterwarnings("ignore")
RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

# Project root
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
PROC_DIR = PROJECT_ROOT / "data" / "processed"
TAB_DIR = PROJECT_ROOT / "results" / "tables"
FIG_DIR = PROJECT_ROOT / "results" / "figures" / "ML"

FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)


def load_data():
    """Load expression matrix and metadata from processed CSV/RDS."""
    expr_csv = PROC_DIR / "expression_matrix.csv"
    meta_csv = PROC_DIR / "metadata.csv"

    if not expr_csv.exists() or not meta_csv.exists():
        print(
            "Hata: data/processed/ altinda expression_matrix.csv ve metadata.csv bulunamadi.\n"
            "Once R/01_preprocessing.R scriptini calistirin."
        )
        sys.exit(1)

    print("Ifade matrisi yukleniyor...")
    expr = pd.read_csv(expr_csv, index_col=0)
    meta = pd.read_csv(meta_csv)

    if "sample_id" in meta.columns:
        meta = meta.set_index("sample_id")

    common = expr.columns.intersection(meta.index)
    expr = expr[common]
    meta = meta.loc[common]

    meta["er_status"] = meta["er_status"].astype(str).str.strip()
    valid = meta["er_status"].isin(["ER+", "ER-"])
    expr = expr.loc[:, valid]
    meta = meta.loc[valid]

    y = (meta["er_status"] == "ER+").astype(int).values  # ER+ = 1, ER- = 0
    X = expr.T.values
    feature_names = expr.index.tolist()
    sample_ids = expr.columns.tolist()

    print(f"Ornek: {X.shape[0]}, Gen: {X.shape[1]}")
    print(f"ER+ (1): {y.sum()}, ER- (0): {(1 - y).sum()}")
    return X, y, feature_names, sample_ids


def setup1_features(feature_names):
    """Union of top DGE genes and WGCNA hub genes."""
    dge_path = TAB_DIR / "DGE_results.csv"
    hub_path = TAB_DIR / "WGCNA_hub_genes.csv"
    genes = set()

    if dge_path.exists():
        dge = pd.read_csv(dge_path)
        sig = (dge["adj.P.Val"] < 0.05) & (dge["logFC"].abs() > 1)
        up = dge.loc[sig & (dge["logFC"] > 1), "gene"].tolist()
        down = dge.loc[sig & (dge["logFC"] < -1), "gene"].tolist()
        genes.update(up)
        genes.update(down)
        print(f"Setup-1: DGE genleri — up {len(up)}, down {len(down)}")
    else:
        print("Uyari: DGE_results.csv yok; once 02_DGE.R calistirin.")

    if hub_path.exists():
        hub = pd.read_csv(hub_path)
        gene_col = "gene" if "gene" in hub.columns else hub.columns[0] if len(hub.columns) else None
        if gene_col and len(hub) > 0:
            hub_genes = hub[gene_col].dropna().astype(str).tolist()
            genes.update(hub_genes)
            print(f"Setup-1: WGCNA hub genleri — {len(set(hub_genes))} benzersiz")
        else:
            print("Uyari: WGCNA_hub_genes.csv bos; sadece DGE genleri kullanilacak.")
    else:
        print("Uyari: WGCNA_hub_genes.csv yok; once 04_WGCNA.R calistirin.")

    available = [g for g in genes if g in feature_names]
    print(f"Setup-1: toplam {len(available)} feature (matriste mevcut)")
    return available


def subset_features(X, feature_names, selected_genes):
    idx = [feature_names.index(g) for g in selected_genes]
    return X[:, idx], selected_genes


def get_models():
    return {
        "LogisticRegression": LogisticRegression(
            max_iter=1000, solver="lbfgs", random_state=RANDOM_STATE
        ),
        "RandomForest": RandomForestClassifier(
            n_estimators=100, random_state=RANDOM_STATE
        ),
        "SVC": SVC(kernel="rbf", probability=True, random_state=RANDOM_STATE),
    }


def evaluate_setup(setup_name, X, y, feature_names):
    """5-fold stratified CV for all models."""
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    models = get_models()
    rows = []

    for model_name, clf in models.items():
        print(f"\n--- {setup_name} | {model_name} ---")
        pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("clf", clf),
        ])

        accs, precs, recs, f1s = [], [], [], []

        for fold, (train_idx, test_idx) in enumerate(skf.split(X, y), 1):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            pipe.fit(X_train, y_train)
            y_pred = pipe.predict(X_test)
            accs.append(accuracy_score(y_test, y_pred))
            precs.append(precision_score(y_test, y_pred, zero_division=0))
            recs.append(recall_score(y_test, y_pred, zero_division=0))
            f1s.append(f1_score(y_test, y_pred, zero_division=0))
            print(f"  Fold {fold}: acc={accs[-1]:.3f} f1={f1s[-1]:.3f}")

        # OOF probabilities for ROC (cross_val_predict)
        y_prob_all = cross_val_predict(pipe, X, y, cv=skf, method="predict_proba")[:, 1]

        # Confusion matrix: last fold
        train_idx, test_idx = list(skf.split(X, y))[-1]
        pipe.fit(X[train_idx], y[train_idx])
        y_pred_last = pipe.predict(X[test_idx])
        y_test_last = y[test_idx]

        for metric, vals in [
            ("accuracy", accs),
            ("precision", precs),
            ("recall", recs),
            ("f1", f1s),
        ]:
            rows.append({
                "setup": setup_name,
                "model": model_name,
                "metric": metric,
                "mean": np.mean(vals),
                "std": np.std(vals),
            })

        print(
            f"  Ortalama accuracy: {np.mean(accs):.3f} ± {np.std(accs):.3f} | "
            f"F1: {np.mean(f1s):.3f} ± {np.std(f1s):.3f}"
        )

        # Confusion matrix (last fold)
        cm = confusion_matrix(y_test_last, y_pred_last)
        fig, ax = plt.subplots(figsize=(5, 4))
        sns.heatmap(
            cm,
            annot=True,
            fmt="d",
            cmap="Blues",
            ax=ax,
            xticklabels=["ER-", "ER+"],
            yticklabels=["ER-", "ER+"],
        )
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        ax.set_title(f"{setup_name} — {model_name}\nConfusion matrix (last fold)")
        fig.savefig(
            FIG_DIR / f"06_cm_{setup_name}_{model_name}.png",
            dpi=300, bbox_inches="tight",
        )
        plt.close(fig)

        # ROC (OOF probabilities from CV folds)
        fpr, tpr, _ = roc_curve(y, y_prob_all)
        roc_auc = auc(fpr, tpr)
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.plot(fpr, tpr, label=f"AUC = {roc_auc:.3f}")
        ax.plot([0, 1], [0, 1], "k--", alpha=0.5)
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.set_title(f"{setup_name} — {model_name}\nROC (5-fold OOF)")
        ax.legend(loc="lower right")
        fig.savefig(
            FIG_DIR / f"06_roc_{setup_name}_{model_name}.png",
            dpi=300, bbox_inches="tight",
        )
        plt.close(fig)

    return pd.DataFrame(rows)


def main():
    print("=== GSE1456 ML Classification ===\n")
    X_full, y, feature_names, sample_ids = load_data()

    # Setup-0: HVG-3000 baseline
    print("\n=== Setup-0: HVG-3000 (Baseline) ===")
    hvg_csv = PROC_DIR / "expression_matrix_hvg3000.csv"
    if not hvg_csv.exists():
        print("Hata: expression_matrix_hvg3000.csv bulunamadi; once 01_preprocessing.R calistirin.")
        sys.exit(1)
    hvg_df = pd.read_csv(hvg_csv, index_col=0)
    available_hvg = [s for s in sample_ids if s in hvg_df.columns]
    if len(available_hvg) < len(sample_ids):
        print(f"Uyari: {len(sample_ids) - len(available_hvg)} ornek HVG matrisinde yok, atlandı.")
    idx_hvg = [sample_ids.index(s) for s in available_hvg]
    X0 = hvg_df[available_hvg].T.values
    y0 = y[idx_hvg]
    print(f"Setup-0: {X0.shape[0]} ornek, {X0.shape[1]} gen (HVG-3000)")
    res0 = evaluate_setup("Setup0_HVG3000", X0, y0, hvg_df.index.tolist())

    # Setup-1
    bio_genes = setup1_features(feature_names)
    if len(bio_genes) < 2:
        print("Setup-1 icin yeterli gen yok; en az 2 gen gerekli.")
        sys.exit(1)
    X1, _ = subset_features(X_full, feature_names, bio_genes)
    res1 = evaluate_setup("Setup1_Biological", X1, y, bio_genes)

    # Setup-2
    print("\n=== Setup-2: VarianceThreshold + SelectKBest ===")
    selector_pipe = Pipeline([
        ("var", VarianceThreshold(threshold=0.1)),
        ("kbest", SelectKBest(f_classif, k=100)),
    ])
    X2 = selector_pipe.fit_transform(X_full, y)
    print(f"Setup-2: secilen feature sayisi: {X2.shape[1]}")
    res2 = evaluate_setup("Setup2_Statistical", X2, y, feature_names)

    results = pd.concat([res0, res1, res2], ignore_index=True)
    out_path = TAB_DIR / "ML_results.csv"
    results.to_csv(out_path, index=False)
    print(f"\nSonuclar kaydedildi: {out_path}")
    print(results.pivot_table(index=["setup", "model"], columns="metric", values="mean"))
    print("\n06_ML_classification.py tamamlandi.")


if __name__ == "__main__":
    main()
