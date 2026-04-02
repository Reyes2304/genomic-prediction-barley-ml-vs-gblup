from barley_gblup import GBLUPRegressor, simulate_barley_article_dataset
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import numpy as np

sim = simulate_barley_article_dataset(n_snps=3000, random_state=7, return_numpy=True)
X = sim["X"]
mask = sim["phenotype"]["integrated_subset"].to_numpy()
y = sim["phenotype"].loc[mask, "PLH"].to_numpy()
X = X[mask]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=7)
model = GBLUPRegressor(grm="vanraden", optimize_delta=True)
model.fit(X_train, y_train)
pred = model.predict(X_test)
rmse = mean_squared_error(y_test, pred, squared=False)
corr = np.corrcoef(y_test, pred)[0, 1]
print({"rmse": round(float(rmse), 4), "corr": round(float(corr), 4), "delta": round(float(model.delta_), 4)})
