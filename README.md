# Repo reproducible - Nested CV 5x3 (TFM)

Este repositorio contiene los notebooks y artefactos necesarios para reproducir el flujo de trabajo del TFM en cebada con validacion cruzada anidada `5x3`.

## Contenido principal

- Notebooks de estructura y merge:
  - `01_pca20_estructura_y_merge.ipynb`
  - `02_pca50_estructura_y_merge.ipynb`
- Notebooks de modelos unificados:
  - `10_pca20_svr_unificado_nested5x3.ipynb`
  - `11_pca50_svr_unificado_nested5x3.ipynb`
  - `20_pca20_tree_unificado_nested5x3.ipynb`
  - `21_pca50_tree_unificado_nested5x3.ipynb`
  - `40_pca20_ensemble_unificado_nested5x3.ipynb`
  - `40_pca50_ensemble_unificado_nested5x3.ipynb`
  - `45_pca20_blocks_unificado_nested5x3.ipynb`
  - `45_pca50_blocks_unificado_nested5x3.ipynb`
  - `50_snp_direct_xgboost_unificado_nested5x3.ipynb`
  - `53_snp_direct_blocks_unificado_nested5x3.ipynb`
  - `60_snp_direct_gblup_unificado_nested5x3.ipynb`
  - `70_snp_direct_nn_keras_pilot.ipynb`
  - `71_snp_direct_nn_keras_ramularia_mlp_a.ipynb`
  - `72_snp_direct_nn_keras_nested5x3_mlp_a_all_traits.ipynb`
  - `73_snp_direct_nn_keras_pca_nested5x3.ipynb`
- Notebooks de resumen/comparativa:
  - `pca20_summary_nested5x3.ipynb`
  - `pca50_summary_nested5x3.ipynb`
  - `snp_summary_nested5x3.ipynb`
  - `99_comparador_resultados_nested5x3.ipynb`
- Datos y resultados:
  - `pheno_data/`
  - `outputs_nested5x3/`
  - `results/`
- Baseline GBLUP local:
  - `barley_gblup_project/`

## Requisito de dato externo (archivo grande)

El archivo `merged_all_mind0.20_g0.10_mac1_nomaf_pruned_raw.raw` no se incluye en este repo para evitar problemas de subida por tamano.

Descarga (OneDrive): [archivo .raw](https://1drv.ms/i/c/23a3a5fd50a4290f/IQDael9XlGW-SqilTJUJRN40AdqGIngDlOT81TRdVvzHA8c?e=vcNBN2)


Pasos:
1. Descarga el fichero desde el enlace.
2. Verifica que el nombre final sea exactamente `merged_all_mind0.20_g0.10_mac1_nomaf_pruned_raw.raw`.
3. Colocalo en la raiz del repo (`./`) antes de ejecutar notebooks SNP/GBLUP/NN.

Alternativa:
1. Si decides versionarlo en GitHub, usa Git LFS.

## Entorno recomendado

- Python `3.11.x`
- Dependencias minimas usadas en notebooks:
  - `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `seaborn`, `plotly`, `xgboost`, `tensorflow`
- Para GBLUP local:

```bash
pip install -e ./barley_gblup_project
```

## Orden sugerido de ejecucion

1. `01_pca20_estructura_y_merge.ipynb`
2. `02_pca50_estructura_y_merge.ipynb`
3. Notebooks de modelos (`10/11/20/21/40/45/50/53/60/70/71/72/73`)
4. Resumenes (`pca20_summary`, `pca50_summary`, `snp_summary`)
5. Comparador final `99_comparador_resultados_nested5x3.ipynb`

