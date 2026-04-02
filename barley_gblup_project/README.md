# barley-gblup

Mini library for genomic BLUP in Python with a scikit-learn style API, plus an article-inspired barley simulator.

## Main objects

- `barley_gblup.GBLUPRegressor`
- `barley_gblup.simulate_barley_article_dataset`

## Quick start

```python
from barley_gblup import GBLUPRegressor, simulate_barley_article_dataset

sim = simulate_barley_article_dataset(n_snps=4000, return_numpy=True)
X = sim["X"]
y = sim["phenotype"]["PLH"].to_numpy()

model = GBLUPRegressor(grm="vanraden", optimize_delta=True)
model.fit(X, y)
pred = model.predict(X[:10])
print(pred)
```

## Notes

The simulator is inspired by the barley article discussed in the conversation:

- 1,110 genotypes
- four groups: `Elite_Spring`, `Elite_Winter`, `PGR_Spring`, `PGR_Winter`
- seven traits: `HD`, `PLH`, `LOD`, `BLU`, `PUC`, `RHY`, `RAM`
- integrated subset flag inspired by the genomic–phenotypic overlap reported in the article

It is synthetic data, not a reconstruction of the original dataset.
