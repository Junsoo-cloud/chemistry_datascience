import numpy as np
import pandas as pd
from sklearn.metrics import r2_score, root_mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns


def regression_metrics(y_true, y_pred, num_data, num_features):
    r2 = r2_score(y_true, y_pred)
    rmse = root_mean_squared_error(y_true, y_pred)
    adj_r2 = 1 - (1 - r2) * (num_data - 1) / (num_data - num_features - 1)
    rss = np.sum((y_true - y_pred)**2)
    rse = np.sqrt(rss / (num_data - num_features - 1))
    return {"r2": float(r2), "rmse": float(rmse), "adj_r2" : float(adj_r2), "rse": float(rse)}


def results_to_df(results):
    """Convert results list/dict to DataFrame for easier saving/plotting.
    results: list of dicts with keys ['selector','experiment','model','r2','rmse']
    """
    return pd.DataFrame(results)


def plot_metric_bar(df, metric, out_path=None, title=None):
    """Plot grouped bar chart of metric by selector and model."""
    plt.figure(figsize=(8, 5))
    sns.barplot(data=df, x="model", y=metric, hue="selector")
    plt.title(title or f"{metric} by model and selector")
    plt.tight_layout()
    if out_path:
        plt.savefig(out_path, dpi=150)
    plt.close()


def plot_feature_importance(series, top_n=30, out_path=None, title=None):
    plt.figure(figsize=(6, min(0.25 * top_n + 2, 12)))
    series.head(top_n).sort_values().plot(kind="barh")
    plt.title(title or "Feature importance (top {})".format(top_n))
    plt.tight_layout()
    if out_path:
        plt.savefig(out_path, dpi=150)
    plt.close()
