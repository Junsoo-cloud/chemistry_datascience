import os
import json
from pathlib import Path

import pandas as pd
import time
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler

from src.data.io import load_train_test
from src.features.logsum_selector import select_features as select_logsum
from src.features.shap_selector import select_shap_features
from src.features.corr_filter import drop_highly_correlated
from src.models.runner import run_baseline, run_optuna_tuning, train_and_eval
from src.features.ols_selector import backward_elimination_ols
from src.utils.eval import results_to_df, plot_metric_bar, plot_feature_importance


def ensure_outdir(outdir):
    Path(outdir).mkdir(parents=True, exist_ok=True)


def run_all(train_path="./src/data/train_new.csv", test_path="./src/data/test_new.csv", target_col="LogP", outdir="outputs", smoke=False):
    ensure_outdir(outdir)

    X_train, X_test, y_train, y_test = load_train_test(train_path, test_path, target_col)

    # Standardize if needed externally; here we rely on models that handle scaling.

    # 1) Log-sum selection
    selected_logsum, coef_series = select_logsum(X_train, y_train, lam=0.01, eps=1e-4, top_k=None, threshold=0.0)
    # save coefficients
    coef_series.to_csv(os.path.join(outdir, "logsum_coefficients_12.csv"))
    plot_feature_importance(coef_series.abs(), top_n=50, out_path=os.path.join(outdir, "logsum_top30_12.png"), title="Log-sum abs coefficients")

    # 2) SHAP selection (top N)
    num_features = 50
    if smoke:
        num_features = 20
    selected_shap, shap_series = select_shap_features(X_train, y_train, num_features=num_features)
    shap_series.to_csv(os.path.join(outdir, "shap_importance_12.csv"))
    plot_feature_importance(shap_series, top_n=30, out_path=os.path.join(outdir, "shap_top30_12.png"), title="SHAP mean abs")

    experiments = []
    feature_selection_records = []

    # For each selector, construct datasets and run baseline + tuned experiments
    selectors = [
        ("logsum", selected_logsum),
        ("shap", selected_shap),
    ]

    # Pre-calc number of tasks to provide ETA
    # For each selector:
    # - Linear models (LinearRegression + Lasso) -> each with scaled/unscaled variants = 2 models * 2 scaling = 4 baseline runs
    # - Tree models baseline: RF + LGBM = 2 runs
    # - Optuna tuning: Lasso (scaled/unscaled) -> 2 runs, RF ->1, LGBM->1 => 4
    per_selector_tasks = 4 + 2 + 4
    total_selectors = len(selectors)
    total_tasks = per_selector_tasks * total_selectors

    # crude avg seconds per task; reduce during smoke
    avg_sec_per_task = 3 if smoke else 20
    estimated_seconds = total_tasks * avg_sec_per_task
    print(f"Planned experiments: {total_tasks} runs. Estimated total time: {estimated_seconds/3600:.2f} hours")

    pbar = tqdm(total=total_tasks, desc="Pipeline progress")

    # Add a "no selection" baseline (use full feature set) as its own selector
    selectors_with_none = [("none", list(X_train.columns))] + selectors

    for selector_name, features in selectors_with_none:
        if len(features) == 0:
            print(f"No features selected for {selector_name}, skipping")
            # advance pbar by per_selector_tasks to keep ETA reasonable
            pbar.update(per_selector_tasks)
            continue

        X_train_sel = X_train[features]
        X_test_sel = X_test[features]

        # record which features are being used for this selector
        feature_selection_records.append({"selector": selector_name, "n_features": len(features), "features": features})

        # --- Linear models: apply correlation filter ---
        corr_thresh = 0.8
        X_train_corr, dropped = drop_highly_correlated(X_train_sel, threshold=corr_thresh)
        X_test_corr = X_test_sel[X_train_corr.columns]

        # 1) Linear Regression (OLS backward elimination)
        # For linear models we apply scaling at the top-level only (per user request)
        selected_ols, ols_model = backward_elimination_ols(X_train_corr, y_train, significance_level=0.05)
        if len(selected_ols) == 0:
            # If OLS dropped everything, fallback to using all corr-filtered features
            selected_ols = list(X_train_corr.columns)

        lr = __import__("sklearn.linear_model", fromlist=["LinearRegression"]).LinearRegression()
        # scale for linear regression
        scaler = StandardScaler()
        X_train_lr = scaler.fit_transform(X_train_corr[selected_ols].values)
        X_test_lr = scaler.transform(X_test_corr[selected_ols].values)
        metrics_lr = train_and_eval(lr, X_train_lr, X_test_lr, y_train, y_test)
        experiments.append({"selector": selector_name, "experiment": "baseline", "model": "LinearRegression_OLS", "selected_features": selected_ols, **metrics_lr})
        pbar.update(1)

        # 2) Lasso baseline
        # Lasso does its own feature selection; user requested to keep Lasso's behavior separate.
        # Apply scaling for Lasso at the top-level only.
        scaler = StandardScaler()
        X_train_lasso = scaler.fit_transform(X_train_corr)
        X_test_lasso = scaler.transform(X_test_corr)
        lasso = __import__("sklearn.linear_model", fromlist=["Lasso"]).Lasso(alpha=0.01, max_iter=5000)
        metrics_lasso = train_and_eval(lasso, X_train_lasso, X_test_lasso, y_train, y_test)
        experiments.append({"selector": selector_name, "experiment": "baseline", "model": "Lasso_scaled", "selected_features": list(X_train_corr.columns), **metrics_lasso})
        pbar.update(1)

        # --- Tree model baselines (no correlation drop or scaling necessary) ---
        from sklearn.ensemble import RandomForestRegressor
        from lightgbm import LGBMRegressor

        rf = RandomForestRegressor(random_state=42)
        metrics_rf = train_and_eval(rf, X_train_sel.values, X_test_sel.values, y_train, y_test)
        experiments.append({"selector": selector_name, "experiment": "baseline", "model": "RandomForest", "selected_features": features, **metrics_rf})
        pbar.update(1)

        lgb = LGBMRegressor(random_state=42)
        metrics_lgb = train_and_eval(lgb, X_train_sel.values, X_test_sel.values, y_train, y_test)
        experiments.append({"selector": selector_name, "experiment": "baseline", "model": "LightGBM", "selected_features": features, **metrics_lgb})
        pbar.update(1)

        # --- Optuna tuning: RF and LGBM only ---
        trials = 4 if smoke else 30

        # RF
        tuned_rf = run_optuna_tuning(X_train_sel.values, y_train, X_test_sel.values, y_test, model_name="RandomForest", n_trials=trials)
        tuned_rf_record = {"selector": selector_name, "experiment": "optuna", "model": "RandomForest", "selected_features": features, **tuned_rf}
        experiments.append(tuned_rf_record)
        pbar.update(1)

        # LGBM
        tuned_lgb = run_optuna_tuning(X_train_sel.values, y_train, X_test_sel.values, y_test, model_name="LightGBM", n_trials=trials)
        tuned_lgb_record = {"selector": selector_name, "experiment": "optuna", "model": "LightGBM", "selected_features": features, **tuned_lgb}
        experiments.append(tuned_lgb_record)
        pbar.update(1)

    pbar.close()

    df = results_to_df(experiments)
    df.to_csv(os.path.join(outdir, "experiments_results_251212.csv"), index=False)

    # Plot metrics
    plot_metric_bar(df, "r2", out_path=os.path.join(outdir, "r2_by_model_selector_12.png"), title="R2 by model and selector")
    plot_metric_bar(df, "rmse", out_path=os.path.join(outdir, "rmse_by_model_selector_12.png"), title="RMSE by model and selector")

    # Save summary JSON
    with open(os.path.join(outdir, "experiments_summary_12.json"), "w") as f:
        json.dump({"n_experiments": len(experiments)}, f, indent=2)

    # Save feature selection summary (which features were used per selector)
    feat_df = pd.DataFrame([{"selector": r["selector"], "n_features": r["n_features"], "features": json.dumps(r["features"])} for r in feature_selection_records])
    feat_df.to_csv(os.path.join(outdir, "feature_selection_summary.csv"), index=False)
    with open(os.path.join(outdir, "feature_selection_summary.json"), "w") as f:
        json.dump(feature_selection_records, f, indent=2)

    print(f"Pipeline completed. Results in: {outdir}")


if __name__ == "__main__":
    # quick smoke run if executed directly
    run_all()
