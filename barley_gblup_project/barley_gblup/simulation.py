from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd


TRAIT_CONFIG = {
    "HD": {"h2": 0.72, "n_qtl": 180, "kind": "continuous"},
    "PLH": {"h2": 0.68, "n_qtl": 220, "kind": "continuous"},
    "LOD": {"h2": 0.78, "n_qtl": 120, "kind": "ordinal"},
    "BLU": {"h2": 0.60, "n_qtl": 140, "kind": "ordinal"},
    "PUC": {"h2": 0.56, "n_qtl": 140, "kind": "ordinal"},
    "RHY": {"h2": 0.18, "n_qtl": 90, "kind": "ordinal"},
    "RAM": {"h2": 0.08, "n_qtl": 80, "kind": "ordinal"},
}

GROUPS = [
    ("Elite_Spring", 128),
    ("Elite_Winter", 170),
    ("PGR_Spring", 288),
    ("PGR_Winter", 524),
]


def _sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


def _make_ld_block(
    rng: np.random.Generator,
    n: int,
    size: int,
    af_low: float,
    af_high: float,
    rho: float,
) -> np.ndarray:
    """Generate one SNP block with approximate LD using a shared latent factor.

    This is intentionally approximate: it aims for realistic-looking marker correlation
    and allele-frequency ranges without reproducing the original data.
    """
    shared = rng.normal(size=(n, 1)).astype(np.float32)
    marker_noise = rng.normal(size=(1, size)).astype(np.float32)
    individual_noise = rng.normal(size=(n, size)).astype(np.float32)

    base_p = rng.uniform(af_low, af_high, size=(1, size)).astype(np.float32)
    base_logit = np.log(base_p / (1.0 - base_p))

    logits = base_logit + 0.8 * math.sqrt(rho) * shared + 0.35 * math.sqrt(rho) * marker_noise + 0.25 * math.sqrt(max(1e-8, 1.0 - rho)) * individual_noise
    probs = np.clip(_sigmoid(logits), 1e-4, 1 - 1e-4)
    return rng.binomial(2, probs).astype(np.float32)


def _simulate_group(
    rng: np.random.Generator,
    n: int,
    n_snps: int,
    ld_block_size: int,
    af_low: float,
    af_high: float,
    rho: float,
) -> np.ndarray:
    blocks = []
    full_blocks = n_snps // ld_block_size
    remainder = n_snps % ld_block_size
    for _ in range(full_blocks):
        blocks.append(_make_ld_block(rng, n, ld_block_size, af_low, af_high, rho))
    if remainder:
        blocks.append(_make_ld_block(rng, n, remainder, af_low, af_high, rho))
    return np.concatenate(blocks, axis=1)


def _ordinalize(x: np.ndarray, n_classes: int = 9) -> np.ndarray:
    qs = np.quantile(x, np.linspace(0, 1, n_classes + 1)[1:-1])
    return np.digitize(x, qs) + 1


def simulate_barley_article_dataset(
    n_snps: int = 8000,
    random_state: int | None = 42,
    return_numpy: bool = False,
) -> dict[str, Any]:
    """Simulate an article-inspired barley genomic prediction dataset.

    The defaults are chosen to be fairly large but still comfortable on a laptop.
    The output is synthetic and should only be used for software testing or demonstrations.
    """
    rng = np.random.default_rng(random_state)

    n_total = sum(n for _, n in GROUPS)
    ld_block_size = 40

    group_params = {
        "Elite_Spring": {"af_low": 0.12, "af_high": 0.45, "rho": 0.88},
        "Elite_Winter": {"af_low": 0.10, "af_high": 0.48, "rho": 0.84},
        "PGR_Spring": {"af_low": 0.05, "af_high": 0.50, "rho": 0.66},
        "PGR_Winter": {"af_low": 0.05, "af_high": 0.50, "rho": 0.70},
    }

    X_parts = []
    rows = []
    start = 0
    for group_name, n in GROUPS:
        params = group_params[group_name]
        X_g = _simulate_group(
            rng=rng,
            n=n,
            n_snps=n_snps,
            ld_block_size=ld_block_size,
            af_low=params["af_low"],
            af_high=params["af_high"],
            rho=params["rho"],
        )
        X_parts.append(X_g)
        season = "Spring" if "Spring" in group_name else "Winter"
        panel = "Elite" if "Elite" in group_name else "PGR"
        end = start + n
        for i in range(start, end):
            rows.append({"sample_id": f"G{i+1:04d}", "group": group_name, "panel": panel, "season": season})
        start = end

    X = np.vstack(X_parts).astype(np.float32)
    metadata = pd.DataFrame(rows)

    # Low-rank summary covariates, analogous to genotype PCs.
    X_slice = X[:, : min(800, X.shape[1])]
    X_centered = X_slice - X_slice.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(X_centered, full_matrices=False)
    pcs = u[:, :5] * s[:5]
    for j in range(5):
        metadata[f"PC{j+1}"] = pcs[:, j]

    panel_effect = np.where(metadata["panel"].to_numpy() == "Elite", -0.35, 0.25)
    season_effect = np.where(metadata["season"].to_numpy() == "Winter", 0.20, -0.05)

    phenotypes = metadata[["sample_id", "group", "panel", "season"]].copy()
    qtl_library = {}

    for trait_name, cfg in TRAIT_CONFIG.items():
        qtl_idx = rng.choice(n_snps, size=min(cfg["n_qtl"], n_snps), replace=False)
        effects = rng.normal(0.0, 1.0 / math.sqrt(len(qtl_idx)), size=len(qtl_idx))
        signal = X[:, qtl_idx] @ effects

        trait_shift = np.zeros(n_total, dtype=float)
        if trait_name in {"BLU", "PUC", "RHY", "RAM"}:
            trait_shift += 0.6 * panel_effect
            if trait_name in {"RHY", "RAM"}:
                trait_shift += 0.5 * season_effect
        else:
            if trait_name == "HD":
                trait_shift += 0.8 * season_effect
            if trait_name == "PLH":
                trait_shift += 0.3 * panel_effect
            if trait_name == "LOD":
                trait_shift += 0.4 * panel_effect + 0.2 * season_effect

        latent = signal + trait_shift + 0.3 * pcs[:, 0] - 0.15 * pcs[:, 1]
        var_signal = np.var(latent)
        noise_var = var_signal * max(1e-6, (1.0 - cfg["h2"]) / cfg["h2"])
        y_cont = latent + rng.normal(0.0, np.sqrt(noise_var), size=n_total)

        y = _ordinalize(y_cont, n_classes=9).astype(int) if cfg["kind"] == "ordinal" else y_cont.astype(float)
        phenotypes[trait_name] = y
        qtl_library[trait_name] = {"indices": qtl_idx, "effects": effects, "target_h2": cfg["h2"]}

    availability = np.zeros(n_total, dtype=bool)
    spring_idx = metadata.index[metadata["season"] == "Spring"].to_numpy()
    winter_idx = metadata.index[metadata["season"] == "Winter"].to_numpy()
    availability[rng.choice(spring_idx, size=min(652, len(spring_idx)), replace=False)] = True
    availability[rng.choice(winter_idx, size=min(458, len(winter_idx)), replace=False)] = True
    phenotypes["integrated_subset"] = availability

    snp_columns = [f"SNP_{j+1:05d}" for j in range(n_snps)]
    genotype_df = pd.DataFrame(X, columns=snp_columns)
    genotype_df.insert(0, "sample_id", metadata["sample_id"])

    out = {
        "genotype": genotype_df,
        "phenotype": phenotypes,
        "metadata": metadata,
        "qtl_library": qtl_library,
        "trait_config": TRAIT_CONFIG,
        "article_notes": {
            "total_genotypes": n_total,
            "group_counts": dict(GROUPS),
            "integrated_subset_sizes": {"Spring": 652, "Winter": 458},
            "source_inspiration": "Yuan et al., GigaScience 2025, PMCID: PMC11811526",
        },
    }
    if return_numpy:
        out["X"] = X
        out["sample_ids"] = metadata["sample_id"].to_numpy()
    return out
