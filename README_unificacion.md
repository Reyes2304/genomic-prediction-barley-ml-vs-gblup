# Notebooks Unificados (Nested CV 5x3)

Carpeta de trabajo con notebooks estandarizados para comparar representaciones PCA20, PCA50, SNP direct y referencia GBLUP, manteniendo trazabilidad de resultados.

## Notebooks principales
- `10_pca20_svr_unificado_nested5x3.ipynb`
- `11_pca50_svr_unificado_nested5x3.ipynb`
- `20_pca20_tree_unificado_nested5x3.ipynb`
- `21_pca50_tree_unificado_nested5x3.ipynb`
- `40_pca20_ensemble_unificado_nested5x3.ipynb`
- `40_pca50_ensemble_unificado_nested5x3.ipynb`
- `45_pca20_blocks_unificado_nested5x3.ipynb`
- `45_pca50_blocks_unificado_nested5x3.ipynb`
- `53_snp_direct_blocks_unificado_nested5x3.ipynb`
- `60_snp_direct_gblup_unificado_nested5x3.ipynb`
- `70_snp_direct_nn_keras_pilot.ipynb`
- `99_comparador_resultados_nested5x3.ipynb`

## Notebooks de resumen comparables (estilo comun)
- `pca20_summary_nested5x3.ipynb`
- `pca50_summary_nested5x3.ipynb`
- `snp_summary_nested5x3.ipynb`

## Salidas clave
- `outputs_nested5x3/`:
  - salidas por run unificado (`summary_nested5x3.csv`, `fold_metrics_nested5x3.csv`, etc.)
  - `_comparativas/` con exportes consolidados del notebook 99.
- `results/pca20_summary_nested5x3/`:
  - master, model summary, heatmap y barplot para PCA20.
- `results/pca50_summary_nested5x3/`:
  - master, model summary, heatmap y barplot para PCA50.
- `results/snp_summary_nested5x3/`:
  - master, model summary, heatmap y barplot para SNP direct.
- `results/snp_summary/`:
  - export de compatibilidad (`snp_model_summary.xlsx/csv`).
- `results/snp_blocks_notebook_nested5x3/` y `results/snp_blocks_notebook_nested5x3_gblup/`:
  - resultados legacy SNP direct y GBLUP usados por el comparador 99.

## Actualizacion reciente en notebook 99
- Se anadio bloque `GBLUP vs Mejor ML Por Trait` que calcula:
  - mejor ML por trait (pool global de PCA20 + PCA50 + SNP direct),
  - top-3 ML por trait,
  - scatter `GBLUP vs mejor ML` con etiquetas por trait y diagonal `y=x`.
- Exportes nuevos en `_comparativas/`:
  - `best_ml_by_trait_all_sources.csv`
  - `top3_ml_by_trait_all_sources.csv`
  - `gblup_vs_best_ml_by_trait.csv`
  - `scatter_gblup_vs_best_ml_pearson.png`
  - `gblup_vs_best_ml_report.xlsx`
