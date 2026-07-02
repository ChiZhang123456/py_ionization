from __future__ import annotations

import numpy as np

from .cross_section import photoionization_cross_section


PLANCK_J_S = 6.626e-34
LIGHT_SPEED_M_S = 3.0e8


def photoionization_rate(wavelength_nm, irradiance_w_m2_nm, species: str = "H"):
    """Calculate a photoionization rate from spectral irradiance.

    Parameters
    ----------
    wavelength_nm
        1D wavelength grid in nm.
    irradiance_w_m2_nm
        Solar spectral irradiance in W m^-2 nm^-1. The last dimension must
        match ``wavelength_nm``.
    species
        Neutral target species passed to
        :func:`photoionization_cross_section`.

    Returns
    -------
    rate_s_1
        Photoionization rate in s^-1 for one neutral particle.

    Formula
    -------
    ``rate = integral sigma(lambda) Phi(lambda) d lambda``, where
    ``Phi(lambda) = I(lambda) / (h c / lambda)``.
    """
    wavelength = np.asarray(wavelength_nm, dtype=float)
    if wavelength.ndim != 1:
        raise ValueError("wavelength_nm must be a 1D array.")
    if wavelength.size < 2:
        raise ValueError("wavelength_nm must contain at least two samples.")

    irradiance = np.asarray(irradiance_w_m2_nm, dtype=float)
    if irradiance.shape[-1] != wavelength.size:
        raise ValueError("The last dimension of irradiance_w_m2_nm must match wavelength_nm.")

    wavelength_m = wavelength * 1e-9
    irradiance_w_m3 = irradiance * 1e9
    photon_energy_j = PLANCK_J_S * LIGHT_SPEED_M_S / wavelength_m
    photon_flux_m3_s = irradiance_w_m3 / photon_energy_j
    dlambda_m = np.nanmean(np.diff(wavelength_m))
    sigma_m2 = photoionization_cross_section(wavelength, species=species) * 1e-4
    rate = np.nansum(sigma_m2 * photon_flux_m3_s * dlambda_m, axis=-1)

    if irradiance.ndim == 1:
        return float(np.asarray(rate))
    return rate


ph_prod_rate = photoionization_rate
production_rate = photoionization_rate

__all__ = ["photoionization_rate", "ph_prod_rate", "production_rate"]
