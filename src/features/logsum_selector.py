import numpy as np
import pandas as pd


def logsum_regression(X, y, lam=0.01, eps=1e-6, max_iter=1000, tol=1e-6):
    """Coordinate-wise update for log-sum penalized regression.

    X: numpy array (n_samples, n_features)
    y: array-like (n_samples,)
    returns beta (numpy array of coefficients)
    """
    n, p = X.shape
    beta = np.zeros(p)

    for iteration in range(max_iter):
        beta_old = beta.copy()

        for j in range(p):
            residual = y - X.dot(beta) + X[:, j] * beta[j]
            wj = np.dot(X[:, j], residual) / n

            c1 = np.clip(wj - eps, -1e3, 1e3)
            c2 = np.clip(c1 ** 2 - 4 * (lam - wj * eps), -1e3, 1e3)

            if c2 > 0:
                beta[j] = np.sign(wj) * (c1 + np.sqrt(c2)) / 2
            else:
                beta[j] = 0

        beta = np.nan_to_num(beta, nan=0.0, posinf=1e6, neginf=-1e6)

        if np.linalg.norm(beta - beta_old) < tol:
            break

    return beta


def select_features(X_df, y, lam=0.01, eps=1e-4, top_k=None, threshold=0.0):
    """Run log-sum selection on X_df (DataFrame) and y (array-like).

    Returns: selected_feature_names (list), coefficients (pd.Series indexed by feature)
    """
    X = X_df.values
    beta = logsum_regression(X, y, lam=lam, eps=eps)

    coef_series = pd.Series(beta, index=X_df.columns)

    if top_k is not None:
        selected = coef_series.abs().sort_values(ascending=False).head(top_k).index.tolist()
    else:
        selected = coef_series[coef_series.abs() > threshold].index.tolist()

    return selected, coef_series
