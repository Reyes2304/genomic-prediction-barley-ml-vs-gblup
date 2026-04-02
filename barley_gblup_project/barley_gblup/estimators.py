from __future__ import annotations

import numpy as np
from scipy.linalg import cho_factor, cho_solve, eigh
from scipy.optimize import minimize_scalar
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.metrics import r2_score
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from threadpoolctl import threadpool_limits


class GBLUPRegressor(BaseEstimator, RegressorMixin):
    """Genomic BLUP regressor with a scikit-learn style API.

    Parameters
    ----------
    grm : {'vanraden', 'standardized'}, default='vanraden'
        Strategy used to construct the genomic relationship matrix.
    fit_intercept : bool, default=True
        Whether to include an intercept in the fixed-effects design.
    optimize_delta : bool, default=True
        If True, estimate delta = sigma_e^2 / sigma_g^2 by REML.
        Otherwise, use ``delta_fixed``.
    delta_fixed : float, default=1.0
        Fixed delta value used when ``optimize_delta=False``.
    alpha_bounds : tuple[float, float], default=(-10.0, 10.0)
        Bounds on log(delta) for the scalar REML optimisation.
    jitter : float, default=1e-8
        Small diagonal stabiliser added to the genomic relationship matrix.
    copy_X : bool, default=True
        Whether to copy X during validation.
    """

    def __init__(
        self,
        grm: str = "vanraden",
        fit_intercept: bool = True,
        optimize_delta: bool = True,
        delta_fixed: float = 1.0,
        alpha_bounds: tuple[float, float] = (-10.0, 10.0),
        jitter: float = 1e-8,
        copy_X: bool = True,
    ):
        self.grm = grm
        self.fit_intercept = fit_intercept
        self.optimize_delta = optimize_delta
        self.delta_fixed = delta_fixed
        self.alpha_bounds = alpha_bounds
        self.jitter = jitter
        self.copy_X = copy_X

    def _compute_grm_train(self, X: np.ndarray):
        if self.grm == "vanraden":
            p = X.mean(axis=0) / 2.0
            W = X - 2.0 * p
            denom = 2.0 * np.sum(p * (1.0 - p))
            if denom <= 0:
                raise ValueError("Non-positive denominator while building the genomic relationship matrix.")
            G = (W @ W.T) / denom
            aux = {"p": p, "denom": denom}
            return G, aux

        if self.grm == "standardized":
            mean_ = X.mean(axis=0)
            scale_ = X.std(axis=0, ddof=0)
            scale_[scale_ == 0.0] = 1.0
            Z = (X - mean_) / scale_
            G = (Z @ Z.T) / X.shape[1]
            aux = {"mean_": mean_, "scale_": scale_}
            return G, aux

        raise ValueError("grm must be either 'vanraden' or 'standardized'.")

    def _compute_grm_cross(self, X_new: np.ndarray, X_train: np.ndarray) -> np.ndarray:
        if self.grm == "vanraden":
            p = self._grm_aux_["p"]
            denom = self._grm_aux_["denom"]
            W_new = X_new - 2.0 * p
            W_train = X_train - 2.0 * p
            return (W_new @ W_train.T) / denom

        mean_ = self._grm_aux_["mean_"]
        scale_ = self._grm_aux_["scale_"]
        Z_new = (X_new - mean_) / scale_
        Z_train = (X_train - mean_) / scale_
        return (Z_new @ Z_train.T) / X_train.shape[1]

    @staticmethod
    def _reml_objective(log_delta: float, Uy: np.ndarray, UX: np.ndarray, eigvals: np.ndarray) -> float:
        delta = np.exp(log_delta)
        s = eigvals + delta

        Uy_w = Uy / np.sqrt(s)
        UX_w = UX / np.sqrt(s[:, None]) if UX.shape[1] else UX

        if UX.shape[1]:
            XtVinvX = UX_w.T @ UX_w
            try:
                c, low = cho_factor(XtVinvX, lower=True, check_finite=False)
                beta = cho_solve((c, low), UX_w.T @ Uy_w, check_finite=False)
            except np.linalg.LinAlgError:
                return np.inf
            resid = Uy_w - UX_w @ beta
            sign, logdet_xtvinvx = np.linalg.slogdet(XtVinvX)
            if sign <= 0:
                return np.inf
        else:
            resid = Uy_w
            logdet_xtvinvx = 0.0

        rss = float(np.sum(resid * resid))
        n = Uy.shape[0]
        p = UX.shape[1]
        if rss <= 0 or n <= p:
            return np.inf

        logdet_v = np.sum(np.log(s))
        return 0.5 * (logdet_v + logdet_xtvinvx + (n - p) * np.log(rss / (n - p)))

    def fit(self, X: np.ndarray, y: np.ndarray):
        X, y = check_X_y(
            X,
            y,
            accept_sparse=False,
            ensure_2d=True,
            y_numeric=True,
            copy=self.copy_X,
        )
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float).reshape(-1)

        self.n_features_in_ = X.shape[1]
        self.X_train_ = X

        G, self._grm_aux_ = self._compute_grm_train(X)
        G = 0.5 * (G + G.T)
        G += self.jitter * np.eye(G.shape[0])

        n = X.shape[0]
        X_fixed = np.ones((n, 1), dtype=float) if self.fit_intercept else np.empty((n, 0), dtype=float)

        with threadpool_limits(limits=1):
            eigvals, eigvecs = eigh(G, check_finite=False)
        eigvals = np.maximum(eigvals, 0.0)
        Uy = eigvecs.T @ y
        UX = eigvecs.T @ X_fixed

        if self.optimize_delta:
            opt = minimize_scalar(
                self._reml_objective,
                bounds=self.alpha_bounds,
                method="bounded",
                args=(Uy, UX, eigvals),
                options={"xatol": 1e-5},
            )
            if not opt.success:
                raise RuntimeError(f"REML optimisation failed: {opt.message}")
            self.delta_ = float(np.exp(opt.x))
        else:
            if self.delta_fixed <= 0:
                raise ValueError("delta_fixed must be > 0.")
            self.delta_ = float(self.delta_fixed)

        s = eigvals + self.delta_
        Uy_w = Uy / np.sqrt(s)
        UX_w = UX / np.sqrt(s[:, None]) if X_fixed.shape[1] else X_fixed

        if X_fixed.shape[1]:
            XtVinvX = UX_w.T @ UX_w
            c, low = cho_factor(XtVinvX, lower=True, check_finite=False)
            self.beta_ = cho_solve((c, low), UX_w.T @ Uy_w, check_finite=False)
        else:
            self.beta_ = np.zeros((0,), dtype=float)

        y_minus_Xb = y - X_fixed @ self.beta_
        Uy_minus = eigvecs.T @ y_minus_Xb
        Vinv_y_minus = eigvecs @ (Uy_minus / s)

        dof = n - X_fixed.shape[1]
        self.sigma_g2_ = float(np.sum(y_minus_Xb * Vinv_y_minus) / dof)
        self.sigma_e2_ = float(self.delta_ * self.sigma_g2_)

        self.G_train_ = G
        self.X_fixed_train_ = X_fixed
        self.y_train_ = y.copy()
        self.eigvals_ = eigvals
        self.eigvecs_ = eigvecs

        with threadpool_limits(limits=1):
            self.alpha_ = np.linalg.solve(G + self.delta_ * np.eye(n), y_minus_Xb)
        self.u_ = G @ self.alpha_
        self.fitted_ = X_fixed @ self.beta_ + self.u_
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        check_is_fitted(self, ["X_train_", "G_train_", "beta_", "alpha_", "delta_"])
        X = check_array(X, accept_sparse=False, ensure_2d=True)
        X = np.asarray(X, dtype=float)

        if X.shape[1] != self.n_features_in_:
            raise ValueError(
                f"X has {X.shape[1]} features, but the model expects {self.n_features_in_}."
            )

        G_cross = self._compute_grm_cross(X, self.X_train_)
        X_fixed_new = np.ones((X.shape[0], 1), dtype=float) if self.fit_intercept else np.empty((X.shape[0], 0), dtype=float)
        y_pred = X_fixed_new @ self.beta_ + G_cross @ self.alpha_
        return y_pred.ravel()

    def score(self, X: np.ndarray, y: np.ndarray) -> float:
        return r2_score(y, self.predict(X))
