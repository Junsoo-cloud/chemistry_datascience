import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Lasso
from sklearn.metrics import r2_score, mean_squared_error

df_train = pd.read_csv("data/train.csv")
X = df_train.drop(columns="LogP")
y = df_train["LogP"]

feature_names = X.columns
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

lasso = Lasso(alpha=0.01, max_iter=10000, fit_intercept=True)
lasso.fit(X_scaled, y)

df_test = pd.read_csv("data/test.csv")
X_test = df_test.drop(columns="LogP")
y_test = df_test["LogP"]

X_test_scaled = scaler.transform(X_test)

y_pred = lasso.predict(X_test_scaled)

r2 = r2_score(y_test, y_pred)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))

print(f"\n LASSO Results")
print(f"R² = {r2:.4f}")
print(f"RMSE = {rmse:.4f}")

# 8) Coefficient DataFrame 생성
coef_df = pd.DataFrame({
    "Descriptor": feature_names,
    "Coefficient": lasso.coef_
}).sort_values(by="Coefficient", ascending=False)

# 9) Top 10 Positive Descriptors
print("\n Top 10 Positive Descriptors (LASSO):")
print(coef_df.head(10))