import os
import json
from pathlib import Path

import pandas as pd

from src.data.io import load_train_test
from src.features.logsum_selector import select_features as select_logsum
from src.features.shap_selector import select_shap_features
from src.models.runner import run_baseline, run_optuna_tuning
from src.utils.eval import results_to_df, plot_metric_bar, plot_feature_importance


def ensure_outdir(outdir):
    Path(outdir).mkdir(parents=True, exist_ok=True)


def run_all(train_path="./src/data/train.csv", test_path="./src/data/test.csv", target_col="LogP", outdir="outputs", smoke=False):
    ensure_outdir(outdir)

    X_train, X_test, y_train, y_test = load_train_test(train_path, test_path, target_col)

    # Standardize if needed externally; here we rely on models that handle scaling.

    # 1) Log-sum selection
    selected_logsum, coef_series = select_logsum(X_train, y_train, lam=0.01, eps=1e-4, top_k=None, threshold=0.0)
    # save coefficients
    coef_series.to_csv(os.path.join(outdir, "logsum_coefficients.csv"))
    plot_feature_importance(coef_series.abs(), top_n=30, out_path=os.path.join(outdir, "logsum_top30.png"), title="Log-sum abs coefficients")

    # 2) SHAP selection (top N)
    num_features = 50
    if smoke:
        num_features = 10
    selected_shap, shap_series = select_shap_features(X_train, y_train, num_features=num_features)
    shap_series.to_csv(os.path.join(outdir, "shap_importance.csv"))
    plot_feature_importance(shap_series, top_n=30, out_path=os.path.join(outdir, "shap_top30.png"), title="SHAP mean abs")

    experiments = []

    # For each selector, construct datasets and run baseline + tuned experiments
    selectors = [
        ("logsum", selected_logsum),
        ("shap", selected_shap),
    ]

    for selector_name, features in selectors:
        if len(features) == 0:
            print(f"No features selected for {selector_name}, skipping")
            continue

        X_train_sel = X_train[features]
        X_test_sel = X_test[features]

        # Baseline
        baseline_results = run_baseline(X_train_sel, X_test_sel, y_train, y_test)
        for r in baseline_results:
            experiments.append({"selector": selector_name, "experiment": "baseline", **r})

        # Optuna tuning
        trials = 8 if smoke else 30
        # Tune Lasso, RF, LightGBM (we don't tune plain LinearRegression)
        for model_name in ["Lasso", "RandomForest", "LightGBM"]:
            tuned = run_optuna_tuning(X_train_sel.values, y_train, X_test_sel.values, y_test, model_name=model_name, n_trials=trials)
            tuned_record = {"selector": selector_name, "experiment": "optuna", **tuned}
            experiments.append(tuned_record)

    df = results_to_df(experiments)
    df.to_csv(os.path.join(outdir, "experiments_results.csv"), index=False)

    # Plot metrics
    plot_metric_bar(df, "r2", out_path=os.path.join(outdir, "r2_by_model_selector.png"), title="R2 by model and selector")
    plot_metric_bar(df, "rmse", out_path=os.path.join(outdir, "rmse_by_model_selector.png"), title="RMSE by model and selector")

    # Save summary JSON
    with open(os.path.join(outdir, "experiments_summary.json"), "w") as f:
        json.dump({"n_experiments": len(experiments)}, f, indent=2)

    print(f"Pipeline completed. Results in: {outdir}")


if __name__ == "__main__":
    # quick smoke run if executed directly
    run_all(smoke=True)
