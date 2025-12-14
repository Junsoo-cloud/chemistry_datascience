import os
import numpy as np
import pandas as pd

from sklearn.linear_model import LinearRegression, Lasso
from sklearn.ensemble import RandomForestRegressor
from lightgbm import LGBMRegressor
from sklearn.metrics import r2_score, root_mean_squared_error

import optuna

from ..utils.eval import regression_metrics


def _train_and_eval(estimator, X_train, X_test, y_train, y_test):
    # Fit and evaluate on both train and test to produce training metrics as well
    estimator.fit(X_train, y_train)
    y_pred_test = estimator.predict(X_test)
    y_pred_train = estimator.predict(X_train)

    # number of data points (n) and features (p) for adjusted R^2 and RSE
    try:
        num_data_test = int(len(y_test))
    except Exception:
        num_data_test = None

    try:
        if hasattr(X_test, "shape") and len(X_test.shape) >= 2:
            num_features_test = int(X_test.shape[1])
        else:
            num_features_test = 1
    except Exception:
        num_features_test = None

    try:
        num_data_train = int(len(y_train))
    except Exception:
        num_data_train = None

    try:
        if hasattr(X_train, "shape") and len(X_train.shape) >= 2:
            num_features_train = int(X_train.shape[1])
        else:
            num_features_train = 1
    except Exception:
        num_features_train = None

    test_metrics = regression_metrics(y_test, y_pred_test, num_data_test, num_features_test)
    train_metrics = regression_metrics(y_train, y_pred_train, num_data_train, num_features_train)

    # prefix train metrics
    prefixed_train = {f"train_{k}": v for k, v in train_metrics.items()}
    # test metrics keep original keys (r2, rmse, ...)
    return {**prefixed_train, **test_metrics}


def train_and_eval(estimator, X_train, X_test, y_train, y_test):
    """Public wrapper to train an estimator and return metrics."""
    return _train_and_eval(estimator, X_train, X_test, y_train, y_test)


def run_baseline(X_train, X_test, y_train, y_test):
    """Run baseline models without hyperparameter tuning.
    Returns list of result dicts.
    """
    results = []

    # Baseline models: do not embed scaling here. Scaling should be controlled by the caller.
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


def run_optuna_tuning(X, y, X_test, y_test, model_name="RandomForest", n_trials=30, random_state=42, val_size=0.2):
    """Tune hyperparameters with Optuna using an internal train/validation split.

    Optuna trials are evaluated on the internal validation set.
    After tuning, the final model is trained on the entire (X, y) and evaluated on
    the external test set (X_test, y_test) passed by the caller.

    Parameters
    - X, y: arrays for tuning (will be split into train/val)
    - X_test, y_test: external test set used only for final evaluation
    - model_name: one of 'RandomForest', 'LightGBM', or 'Lasso'
    - n_trials: number of Optuna trials
    - val_size: fraction of (X,y) to use as validation during tuning
    """
    from sklearn.model_selection import train_test_split

    X_arr = np.asarray(X)
    y_arr = np.asarray(y)

    # split training data into tuning train / validation
    X_tune_train, X_val, y_tune_train, y_val = train_test_split(
        X_arr, y_arr, test_size=val_size, random_state=random_state
    )

    def objective_rf(trial):
        n_estimators = trial.suggest_int("n_estimators", 50, 300)
        max_depth = trial.suggest_int("max_depth", 3, 30)
        min_samples_split = trial.suggest_int("min_samples_split", 2, 10)
        model = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth,
                                      min_samples_split=min_samples_split, random_state=random_state)
        model.fit(X_tune_train, y_tune_train)
        y_pred = model.predict(X_val)
        r2 = r2_score(y_val, y_pred)
        rmse = root_mean_squared_error(y_val, y_pred)
        return float(r2)
        # rmse = root_mean_squared_error(y_val, y_pred)
        # return float(rmse)

    def objective_lgb(trial):
        params = {
            "num_leaves": trial.suggest_int("num_leaves", 8, 64),
            "learning_rate": trial.suggest_float("learning_rate", 1e-3, 0.3, log=True),
            "n_estimators": trial.suggest_int("n_estimators", 50, 300),
            "max_depth": trial.suggest_int("max_depth", -1, 30),
        }
        model = LGBMRegressor(random_state=random_state, **params)
        model.fit(X_tune_train, y_tune_train)
        y_pred = model.predict(X_val)
        r2 = r2_score(y_val, y_pred)
        rmse = root_mean_squared_error(y_val, y_pred)
        return float(r2)
        # rmse = root_mean_squared_error(y_val, y_pred)
        # return float(rmse)

    def objective_lasso(trial):
        alpha = trial.suggest_float("alpha", 1e-6, 10.0, log=True)
        model = Lasso(alpha=alpha, max_iter=5000)
        model.fit(X_tune_train, y_tune_train)
        y_pred = model.predict(X_val)
        r2 = r2_score(y_val, y_pred)
        return float(r2)
        # rmse = root_mean_squared_error(y_val, y_pred)
        # return float(rmse)

    if model_name == "RandomForest":
        objective = objective_rf
    elif model_name == "LightGBM":
        objective = objective_lgb
    elif model_name == "Lasso":
        objective = objective_lasso
    else:
        raise ValueError("Unsupported model_name for tuning")

    # objectives return R2 (higher is better)
    study = optuna.create_study(direction="maximize")
    study.optimize(objective, n_trials=n_trials)

    # Train final model with best params on the entire training set
    best = study.best_params
    if model_name == "RandomForest":
        final_model = RandomForestRegressor(n_estimators=best.get("n_estimators", 100),
                                            max_depth=best.get("max_depth", None),
                                            min_samples_split=best.get("min_samples_split", 2),
                                            random_state=random_state)
    elif model_name == "LightGBM":
        final_model = LGBMRegressor(random_state=random_state, **best)
    elif model_name == "Lasso":
        # Do not embed scaling here; caller should scale if desired
        final_model = Lasso(alpha=best.get("alpha", 1.0), max_iter=5000)
    else:
        final_model = LGBMRegressor(random_state=random_state, **best)

    # fit on full training data and evaluate on external test set
    final_model.fit(X_arr, y_arr)

    # compute train + test metrics
    metrics = _train_and_eval(final_model, X_arr, X_test, y_arr, y_test)
    return {"model": model_name, "best_params": study.best_params, **metrics}
