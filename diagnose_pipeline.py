import os
import sys
import traceback
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import LinearRegression, Lasso
from sklearn.ensemble import RandomForestRegressor
from lightgbm import LGBMRegressor
from sklearn.preprocessing import StandardScaler
from scipy.stats import ks_2samp

from src.data.io import load_train_test
from src.features.logsum_selector import select_features as select_logsum
from src.features.shap_selector import select_shap_features
from src.features.corr_filter import drop_highly_correlated
from src.features.ols_selector import backward_elimination_ols
from src.utils.eval import regression_metrics


OUTDIR = "./src/outputs/diagnostics"
os.makedirs(OUTDIR, exist_ok=True)


def save_fig(fig, name):
    path = os.path.join(OUTDIR, name)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(train_path="./src/data/train.csv", test_path="./src/data/test.csv", smoke=False):
    report_lines = []
    try:
        X_train, X_test, y_train, y_test = load_train_test(train_path, test_path)
    except Exception:
        report_lines.append("Failed to load train/test: \n" + traceback.format_exc())
        open(os.path.join(OUTDIR, "report.txt"), "w").write("\n".join(report_lines))
        print("Failed to load data—see report")
        return

    report_lines.append(f"Shapes: X_train {X_train.shape}, X_test {X_test.shape}, y_train {y_train.shape}, y_test {y_test.shape}")
    report_lines.append("y_train describe:\n" + str(y_train.describe()))
    report_lines.append("y_test describe:\n" + str(y_test.describe()))

    # basic checks
    report_lines.append(f"NaNs in train (per col count>0): {int((X_train.isna().sum()>0).sum())}, test: {int((X_test.isna().sum()>0).sum())}")
    const_cols = [c for c in X_train.columns if X_train[c].nunique()<=1]
    report_lines.append(f"Constant columns in train: {len(const_cols)}; example: {const_cols[:10]}")
    dupes = int(X_train.duplicated().sum())
    report_lines.append(f"Duplicate rows in train: {dupes}")

    # feature range
    rng = (X_train.max() - X_train.min()).abs()
    top_ranges = rng.sort_values(ascending=False).head(10)
    report_lines.append("Top feature ranges (train):\n" + str(top_ranges))

    # mean diffs
    mean_diff = (X_train.mean() - X_test.mean()).abs().sort_values(ascending=False)
    report_lines.append("Top mean diffs between train/test:\n" + str(mean_diff.head(20)))

    # KS tests for top differing features
    top_feats = mean_diff.head(30).index.tolist()
    ks_results = {}
    for col in top_feats:
        try:
            a = X_train[col].dropna()
            b = X_test[col].dropna()
            if len(a) < 10 or len(b) < 10:
                ks_results[col] = np.nan
                continue
            ks_results[col] = ks_2samp(a, b).pvalue
        except Exception:
            ks_results[col] = np.nan
    ks_series = pd.Series(ks_results).sort_values()
    report_lines.append("KS-test p-values for top differing features:\n" + str(ks_series.head(20)))

    # PCA scatter of train vs test
    try:
        pca = PCA(2)
        Z = pca.fit_transform(pd.concat([X_train.fillna(0), X_test.fillna(0)], axis=0))
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(Z[:len(X_train),0], Z[:len(X_train),1], alpha=0.4, label='train')
        ax.scatter(Z[len(X_train):,0], Z[len(X_train):,1], alpha=0.4, label='test')
        ax.legend()
        ax.set_title('PCA train vs test')
        save_fig(fig, 'pca_train_test.png')
        report_lines.append('Saved PCA plot to pca_train_test.png')
    except Exception:
        report_lines.append('PCA failed:\n' + traceback.format_exc())

    # Log-sum and SHAP selection introspection
    try:
        sel_logsum, coef_series = select_logsum(X_train, y_train, top_k=None)
        report_lines.append(f"Log-sum selected {len(sel_logsum)} features (example 10): {sel_logsum[:10]}")
    except Exception:
        report_lines.append('Log-sum selection failed:\n' + traceback.format_exc())

    try:
        nf = 10 if smoke else 50
        sel_shap, shap_series = select_shap_features(X_train, y_train, num_features=nf)
        report_lines.append(f"SHAP selected {len(sel_shap)} features (top {nf}): {sel_shap[:10]}")
    except Exception:
        report_lines.append('SHAP selection failed:\n' + traceback.format_exc())

    # OLS selection on logsum features
    try:
        if len(sel_logsum) == 0:
            report_lines.append('Log-sum selected zero features, skipping OLS')
        else:
            X_train_sel = X_train[sel_logsum]
            X_test_sel = X_test[sel_logsum]
            X_train_corr, dropped = drop_highly_correlated(X_train_sel, threshold=0.8)
            X_test_corr = X_test_sel[X_train_corr.columns]
            selected_ols, ols_model = backward_elimination_ols(X_train_corr, y_train)
            report_lines.append(f"After corr filter: {X_train_corr.shape[1]} features, dropped {len(dropped)}\nOLS selected {len(selected_ols)} features: {selected_ols[:20]}")
    except Exception:
        report_lines.append('OLS selection failed:\n' + traceback.format_exc())

    # Quick model fits on full data (no selection)
    try:
        report_lines.append('Training quick baseline models on full feature set...')
        # LinearRegression
        lr = LinearRegression()
        lr.fit(X_train.values, y_train)
        ylr = lr.predict(X_test.values)
        metrics_lr = regression_metrics(y_test, ylr, len(y_test), X_test.values.shape[1])
        report_lines.append('LinearRegression full -> ' + str(metrics_lr))

        # Lasso scaled
        scaler = StandardScaler()
        Xtr_s = scaler.fit_transform(X_train)
        Xte_s = scaler.transform(X_test)
        lasso = Lasso(alpha=0.01, max_iter=5000)
        lasso.fit(Xtr_s, y_train)
        yl = lasso.predict(Xte_s)
        metrics_lasso = regression_metrics(y_test, yl, len(y_test), Xte_s.shape[1])
        report_lines.append('Lasso_scaled full -> ' + str(metrics_lasso))

        # RandomForest
        rf = RandomForestRegressor(random_state=42, n_estimators=100)
        rf.fit(X_train.values, y_train)
        yrf = rf.predict(X_test.values)
        metrics_rf = regression_metrics(y_test, yrf, len(y_test), X_test.values.shape[1])
        report_lines.append('RandomForest full -> ' + str(metrics_rf))

        # LGBM
        lgb = LGBMRegressor(random_state=42)
        lgb.fit(X_train.values, y_train)
        ylg = lgb.predict(X_test.values)
        metrics_lgb = regression_metrics(y_test, ylg, len(y_test), X_test.values.shape[1])
        report_lines.append('LightGBM full -> ' + str(metrics_lgb))

        # Dummy baseline
        dum = DummyRegressor(strategy='mean')
        dum.fit(X_train.values, y_train)
        yd = dum.predict(X_test.values)
        metrics_d = regression_metrics(y_test, yd, len(y_test), 1)
        report_lines.append('Dummy mean -> ' + str(metrics_d))

        # Save a comparison scatter RF vs true
        try:
            fig, ax = plt.subplots(figsize=(6,5))
            ax.scatter(y_test, yrf, s=8, alpha=0.5)
            ax.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--')
            ax.set_xlabel('y_test'); ax.set_ylabel('y_pred_rf'); ax.set_title('RF preds vs true')
            save_fig(fig, 'rf_preds_vs_true.png')
            report_lines.append('Saved rf_preds_vs_true.png')
        except Exception:
            report_lines.append('RF scatter save failed:\n' + traceback.format_exc())

    except Exception:
        report_lines.append('Model training block failed:\n' + traceback.format_exc())

    # write report
    report_path = os.path.join(OUTDIR, 'report.txt')
    open(report_path, 'w').write("\n".join(report_lines))
    print(f"Diagnosis complete. Report written to {report_path}")


if __name__ == '__main__':
    smoke = True
    if len(sys.argv) > 1 and sys.argv[1] == '--full':
        smoke = False
    main(smoke=smoke)
