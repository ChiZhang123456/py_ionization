"""Electron impact ionization cross sections and production rates."""

from .cross_section import electron_impact_cross_section
from .production_rate import electron_impact_rate

ei_cross_section = electron_impact_cross_section
ei_prod_rate = electron_impact_rate
cross_section = electron_impact_cross_section
production_rate = electron_impact_rate

__all__ = [
    "electron_impact_cross_section",
    "electron_impact_rate",
    "ei_cross_section",
    "ei_prod_rate",
    "cross_section",
    "production_rate",
]
