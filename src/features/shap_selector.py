import numpy as np
import pandas as pd
import re

def select_shap_features(X_train, y_train, model=None, num_features=50, random_state=42):
    """Select top features by mean absolute SHAP values.

    X_train: DataFrame
    y_train: array-like
    model: optional estimator; if None, a LightGBMRegressor will be trained.
    Returns selected_feature_names, shap_values_series (abs-mean per feature)
    """
    try:
        import lightgbm as lgb
        import shap
    except Exception as e:
        raise ImportError("shap and lightgbm are required for SHAP-based selection. Install them: pip install shap lightgbm")

    if model is None:
        model = lgb.LGBMRegressor(random_state=random_state)

    # LightGBM errors if feature names contain special JSON characters.
    # Create a sanitized copy of the training DataFrame for fitting, but keep
    # original column names to map SHAP importance back to the user-facing names.
    orig_cols = list(X_train.columns)
    safe_cols = [re.sub(r"[^0-9A-Za-z_]", "_", str(c)) for c in orig_cols]
    X_train_safe = X_train.copy()
    X_train_safe.columns = safe_cols

    model.fit(X_train_safe, y_train)

    explainer = shap.TreeExplainer(model)
    shap_vals = explainer.shap_values(X_train_safe)

    # shap_values for regression returns array (n_samples, n_features)
    mean_abs = np.mean(np.abs(shap_vals), axis=0)
    series = pd.Series(mean_abs, index=safe_cols)

    # map index back to original feature names
    series.index = orig_cols
    series = series.sort_values(ascending=False)

    selected = series.head(num_features).index.tolist()

    return selected, series
