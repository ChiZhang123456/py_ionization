"""Photoionization cross sections and production rates."""

from .cross_section import photoionization_cross_section
from .production_rate import photoionization_rate

ph_cross_section = photoionization_cross_section
ph_prod_rate = photoionization_rate
cross_section = photoionization_cross_section
production_rate = photoionization_rate

__all__ = [
    "photoionization_cross_section",
    "photoionization_rate",
    "ph_cross_section",
    "ph_prod_rate",
    "cross_section",
    "production_rate",
]
