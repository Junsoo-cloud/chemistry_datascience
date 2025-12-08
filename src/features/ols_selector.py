import pandas as pd
import statsmodels.api as sm


def backward_elimination_ols(X: pd.DataFrame, y, significance_level: float = 0.05):
    """Perform backward elimination based on p-values using OLS from statsmodels.

    Returns: selected_features (list), summary (statsmodels summary)
    """
    X_work = X.copy()
    # add constant
    X_with_const = sm.add_constant(X_work)
    selected = list(X_work.columns)

    while True:
        model = sm.OLS(y, sm.add_constant(X_work[selected])).fit()
        pvalues = model.pvalues.drop('const', errors='ignore')
        if pvalues.empty:
            break
        max_pval = pvalues.max()
        if max_pval > significance_level:
            # drop the feature with highest p-value
            drop_feature = pvalues.idxmax()
            selected.remove(drop_feature)
        else:
            break

    final_model = sm.OLS(y, sm.add_constant(X_work[selected])).fit()
    return selected, final_model
