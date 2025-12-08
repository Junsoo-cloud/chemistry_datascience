import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.linear_model import LinearRegression, Lasso
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from lightgbm import LGBMRegressor

train = pd.read_csv("train.csv")
test = pd.read_csv("test.csv")

X_train = train.drop(columns='LogP')
y_train = train['LogP']
X_test = test.drop(columns='LogP')
y_test = test['LogP']

X_train = X_train.replace([np.inf, -np.inf], np.nan).fillna(0)
X_test = X_test.replace([np.inf, -np.inf], np.nan).fillna(0)
y_train = np.nan_to_num(y_train)
y_test = np.nan_to_num(y_test)

X_train_np = X_train.values
X_test_np = X_test.values

scaler = StandardScaler()
X_train_np = scaler.fit_transform(X_train_np)
X_test_np = scaler.transform(X_test_np)

def logsum_regression(X, y, lam=0.01, eps=1e-4, max_iter=2000, tol=1e-6):
    n, p = X.shape
    beta = np.zeros(p)

    for iteration in range(max_iter):
        beta_old = beta.copy()

        for j in range(p):
            residual = y - X.dot(beta) + X[:, j] * beta[j]
            wj = np.dot(X[:, j], residual) / n

            c1 = np.clip(wj - eps, -1e3, 1e3)
            c2 = np.clip(c1**2 - 4 * (lam - wj * eps), -1e3, 1e3)

            if c2 > 0:
                beta[j] = np.sign(wj) * (c1 + np.sqrt(c2)) / 2
            else:
                beta[j] = 0

        if np.linalg.norm(beta - beta_old) < tol:
            break

    return beta

beta = logsum_regression(X_train_np, y_train, lam=0.01, eps=1e-4)

selected_idx = np.where(beta != 0)[0]
print(f"선택된 변수 개수: {len(selected_idx)} / {len(beta)}\n")

coef_df = pd.DataFrame({
    "Descriptor": X_train.columns[selected_idx],
    "LogSum Coef": beta[selected_idx]
}).sort_values(by="LogSum Coef", ascending=False)

print("🔹 Top 10 Positive Descriptors")
print(coef_df.head(10))

X_train_sel = X_train_np[:, selected_idx]
X_test_sel = X_test_np[:, selected_idx]

lr = LinearRegression()
lr.fit(X_train_sel, y_train)
y_pred_lr = lr.predict(X_test_sel)

lasso = Lasso(alpha=0.01, max_iter=5000)
lasso.fit(X_train_sel, y_train)
y_pred_lasso = lasso.predict(X_test_sel)

rf = RandomForestRegressor(random_state=42)
rf.fit(X_train_sel, y_train)
y_pred_rf = rf.predict(X_test_sel)

xgb = XGBRegressor(random_state=42, n_jobs=-1)
xgb.fit(X_train_sel, y_train)
y_pred_xgb = xgb.predict(X_test_sel)

lgb = LGBMRegressor(random_state=42)
lgb.fit(X_train_sel, y_train)
y_pred_lgb = lgb.predict(X_test_sel)

def evaluate_model(name, y_true, y_pred):
    print(f"\n {name}")
    print(f"R²: {r2_score(y_true, y_pred):.4f}")
#    print(f"RMSE: {mean_squared_error(y_true, y_pred, squared=False):.4f}")

evaluate_model("Linear Regression", y_test, y_pred_lr)
evaluate_model("LASSO", y_test, y_pred_lasso)
evaluate_model("Random Forest", y_test, y_pred_rf)
evaluate_model("XGBoost", y_test, y_pred_xgb)
evaluate_model("LightGBM", y_test, y_pred_lgb)
