"""Ionization cross sections and production-rate calculations.

The package provides three process-specific subpackages:

* :mod:`py_ionization.photoionization`
* :mod:`py_ionization.electron_impact`
* :mod:`py_ionization.charge_exchange`

All cross sections are returned in cm^2. Production rates are returned in
s^-1 for a neutral target particle when the input fluxes are differential
energy fluxes in eV cm^-2 s^-1 sr^-1 eV^-1.
"""

from .photoionization import photoionization_cross_section, photoionization_rate
from .electron_impact import electron_impact_cross_section, electron_impact_rate
from .charge_exchange import charge_exchange_cross_section, charge_exchange_rate

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "photoionization_cross_section",
    "photoionization_rate",
    "electron_impact_cross_section",
    "electron_impact_rate",
    "charge_exchange_cross_section",
    "charge_exchange_rate",
]
