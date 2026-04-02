# PCA20 nested5x3 consolidation report

## Scope
- Solo fuentes dentro de `nested5x3_unificados/outputs_nested5x3`.
- Orden de modelos fijo por bloques: lineales -> SVR -> arbol/ensembles.
- Heatmap en una sola figura con traits agrupados: enfermedades | agronomicos.

## Additional metrics added
- `median_pearson_r`
- `mean_rank_pearson` (menor es mejor)
- `best_pearson_trait_count`