import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error

df = pd.read_csv("train.csv")
X = df.drop(columns='LogP')
y = df['LogP']

feature_names = X.columns
X = X.replace([np.inf, -np.inf], np.nan).fillna(0)
y = np.nan_to_num(y)
X = X.values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

def logsum_regression(X, y, lam=0.01, eps=1e-6, max_iter=1000, tol=1e-6):
    n, p = X.shape
    beta = np.zeros(p)

    for iteration in range(max_iter):
        beta_old = beta.copy()

        for j in range(p):
            residual = y - X.dot(beta) + X[:, j] * beta[j]
            wj = np.dot(X[:, j], residual) / n

            c1 = np.clip(wj - eps, -1e3, 1e3)
            c2 = np.clip(c1**2 - 4*(lam - wj*eps), -1e3, 1e3)
            if c2 > 0:
                beta[j] = np.sign(wj) * (c1 + np.sqrt(c2)) / 2
            else:
                beta[j] = 0

        # 수치 안정화
        beta = np.nan_to_num(beta, nan=0.0, posinf=1e6, neginf=-1e6)

        # 수렴 조건
        if np.linalg.norm(beta - beta_old) < tol:
            break

    return beta

df = pd.read_csv("test.csv")
X_test = df.drop(columns='LogP')
y_test = df['LogP']

beta = logsum_regression(X, y, lam=0.01, eps=1e-4)
y_pred = X_test.dot(beta)

# 수치 안정화 후 평가
y_pred = np.nan_to_num(y_pred, nan=0.0, posinf=np.max(y), neginf=np.min(y))

r2 = r2_score(y_test, y_pred)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))

print(f"R² = {r2:.4f}, RMSE = {rmse:.4f}")

#유의미한 descriptor 확인
coef_df = pd.DataFrame({
    "Descriptor": feature_names,
    "Coefficient": beta
}).sort_values(by="Coefficient", ascending=False)

print("\n top 10 Positive Descriptors:")
print(coef_df.head(10))