import os
import numpy as np
import pandas as pd

from sklearn.linear_model import LinearRegression, Lasso
from sklearn.ensemble import RandomForestRegressor
from lightgbm import LGBMRegressor
from sklearn.metrics import root_mean_squared_error

import optuna

from ..utils.eval import regression_metrics


def _train_and_eval(estimator, X_train, X_test, y_train, y_test):
    estimator.fit(X_train, y_train)
    y_pred = estimator.predict(X_test)
    # number of data points (n) and features (p) for adjusted R^2 and RSE
    try:
        num_data = int(len(y_test))
    except Exception:
        num_data = None

    try:
        if hasattr(X_test, "shape") and len(X_test.shape) >= 2:
            num_features = int(X_test.shape[1])
        else:
            num_features = 1
    except Exception:
        num_features = None

    return regression_metrics(y_test, y_pred, num_data, num_features)


def train_and_eval(estimator, X_train, X_test, y_train, y_test):
    """Public wrapper to train an estimator and return metrics."""
    return _train_and_eval(estimator, X_train, X_test, y_train, y_test)


def run_baseline(X_train, X_test, y_train, y_test):
    """Run baseline models without hyperparameter tuning.
    Returns list of result dicts.
    """
    results = []

    models = {
        "LinearRegression": LinearRegression(),
        "Lasso": Lasso(alpha=0.01, max_iter=5000),
        "RandomForest": RandomForestRegressor(random_state=42),
        "LightGBM": LGBMRegressor(random_state=42),
    }

    for name, model in models.items():
        metrics = _train_and_eval(model, X_train, X_test, y_train, y_test)
        results.append({"model": name, **metrics})

    return results


def run_optuna_tuning(X, y, X_test, y_test, model_name="RandomForest", n_trials=30, random_state=42):
    """Tune hyperparameters with Optuna and evaluate on the provided test set.

    This function does NOT use cross-validation; it fits on (X, y) and evaluates on (X_test, y_test).
    model_name: one of 'RandomForest', 'LightGBM', 'Lasso'
    """
    def objective_rf(trial):
        n_estimators = trial.suggest_int("n_estimators", 50, 300)
        max_depth = trial.suggest_int("max_depth", 3, 30)
        min_samples_split = trial.suggest_int("min_samples_split", 2, 10)
        model = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth,
                                      min_samples_split=min_samples_split, random_state=random_state)
        model.fit(X, y)
        y_pred = model.predict(X_test)
        rmse = root_mean_squared_error(y_test, y_pred)
        return float(rmse)

    def objective_lgb(trial):
        params = {
            "num_leaves": trial.suggest_int("num_leaves", 8, 64),
            "learning_rate": trial.suggest_float("learning_rate", 1e-3, 0.3, log=True),
            "n_estimators": trial.suggest_int("n_estimators", 50, 300),
            "max_depth": trial.suggest_int("max_depth", -1, 30),
        }
        model = LGBMRegressor(random_state=random_state, **params)
        model.fit(X, y)
        y_pred = model.predict(X_test)
        rmse = root_mean_squared_error(y_test, y_pred)
        return float(rmse)

    if model_name == "RandomForest":
        objective = objective_rf
    elif model_name == "LightGBM":
        objective = objective_lgb
    else:
        raise ValueError("Unsupported model_name for tuning")

    study = optuna.create_study(direction="minimize")
    study.optimize(objective, n_trials=n_trials)

    # Train final model with best params
    if model_name == "RandomForest":
        best = study.best_params
        final_model = RandomForestRegressor(n_estimators=best.get("n_estimators", 100),
                                            max_depth=best.get("max_depth", None),
                                            min_samples_split=best.get("min_samples_split", 2),
                                            random_state=random_state)
    elif model_name == "LightGBM":
        best = study.best_params
        final_model = LGBMRegressor(random_state=random_state, **best)
    else:
        # Should not happen due to earlier check, but keep fallback
        best = study.best_params
        final_model = LGBMRegressor(random_state=random_state, **best)

    metrics = _train_and_eval(final_model, X, X_test, y, y_test)
    return {"model": model_name, "best_params": study.best_params, **metrics}
