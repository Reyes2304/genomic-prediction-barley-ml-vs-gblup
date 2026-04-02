"""barley_gblup: GBLUP tools and article-inspired simulation utilities."""

from .estimators import GBLUPRegressor
from .simulation import simulate_barley_article_dataset

__all__ = ["GBLUPRegressor", "simulate_barley_article_dataset"]
__version__ = "0.1.0"
