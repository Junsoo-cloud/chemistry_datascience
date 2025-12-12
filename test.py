import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from lightgbm import LGBMRegressor
from sklearn.metrics import root_mean_squared_error, r2_score


df_train = pd.read_csv("./src/data/train.csv")
df_test = pd.read_csv("./src/data/test.csv")


X = df_train.drop(columns="LogP")
y = df_train["LogP"]


model = RandomForestRegressor(n_estimators=300, max_depth=20,
                                      min_samples_split=3, random_state=42)
model.fit(X, y)
y_pred = model.predict(X)
r_squared = r2_score(y, y_pred)
print("training r_squared score : ")
print(r_squared)

print("===================================")

X = df_test.drop(columns="LogP")
y = df_test["LogP"]

y_pred_test = model.predict(X)
r_squared = r2_score(y, y_pred)
print("training r_squared score : ")
print(r_squared)


