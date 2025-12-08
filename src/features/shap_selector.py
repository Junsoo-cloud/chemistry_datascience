import numpy as np
import pandas as pd

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

    model.fit(X_train, y_train)

    explainer = shap.TreeExplainer(model)
    shap_vals = explainer.shap_values(X_train)

    # shap_values for regression returns array (n_samples, n_features)
    mean_abs = np.mean(np.abs(shap_vals), axis=0)
    series = pd.Series(mean_abs, index=X_train.columns).sort_values(ascending=False)

    selected = series.head(num_features).index.tolist()

    return selected, series
