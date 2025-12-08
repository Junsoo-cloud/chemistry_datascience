import pandas as pd
import numpy as np


def drop_highly_correlated(X: pd.DataFrame, threshold: float = 0.8):
    """Drop features from X that are highly correlated with others.

    Strategy:
    - Compute absolute correlation matrix.
    - For pairs with correlation > threshold, greedily drop the column with higher mean absolute correlation (more redundant).

    Returns: X_reduced (DataFrame), dropped_features (list)
    """
    if X.shape[1] == 0:
        return X.copy(), []

    corr = X.corr().abs()
    # Mask diagonal
    np.fill_diagonal(corr.values, 0)

    to_drop = set()

    # Compute mean absolute correlation per column
    mean_abs = corr.mean()

    # Iterate over pairs sorted by correlation descending
    # We'll get upper triangle indices
    # Stack and sort
    pairs = (
        corr.stack()
        .reset_index()
        .rename(columns={0: 'corr', 'level_0': 'f1', 'level_1': 'f2'})
        .query('corr > @threshold')
        .sort_values('corr', ascending=False)
    )

    for _, row in pairs.iterrows():
        f1 = row['f1']
        f2 = row['f2']
        if f1 in to_drop or f2 in to_drop:
            continue
        # drop the one with higher mean abs correlation
        if mean_abs[f1] >= mean_abs[f2]:
            to_drop.add(f1)
        else:
            to_drop.add(f2)

    X_reduced = X.drop(columns=list(to_drop)) if to_drop else X.copy()
    return X_reduced, list(to_drop)
