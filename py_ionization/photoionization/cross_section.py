from __future__ import annotations

from functools import lru_cache
from importlib.resources import files

import numpy as np


_SPECIES_TO_FILE = {
    "h": "ph_H_cross_section.txt",
    "p": "ph_H_cross_section.txt",
    "h+": "ph_H_cross_section.txt",
    "o": "ph_O_cross_section.txt",
    "o+": "ph_O_cross_section.txt",
    "co": "ph_CO_cross_section.txt",
    "co+": "ph_CO_cross_section.txt",
    "co2": "ph_CO2_cross_section.txt",
    "co2+": "ph_CO2_cross_section.txt",
    "h2o": "ph_H2O_cross_section.txt",
    "h2o+": "ph_H2O_cross_section.txt",
    "n2": "ph_N2_cross_section.txt",
    "n2+": "ph_N2_cross_section.txt",
    "na": "ph_Na_cross_section.txt",
    "na+": "ph_Na_cross_section.txt",
    "mg": "ph_Mg_cross_section.txt",
    "mg+": "ph_Mg_cross_section.txt",
    "ca": "ph_Ca_cross_section.txt",
    "ca+": "ph_Ca_cross_section.txt",
}


def _cross_section_file(species: str) -> str:
    key = str(species).lower().strip().replace(" ", "")
    try:
        return _SPECIES_TO_FILE[key]
    except KeyError as exc:
        supported = ", ".join(sorted(set(_SPECIES_TO_FILE)))
        raise ValueError(f"Unsupported photoionization species {species!r}. Supported keys: {supported}") from exc


@lru_cache(maxsize=None)
def _load_cross_section(filename: str) -> tuple[np.ndarray, np.ndarray]:
    path = files("py_ionization.photoionization").joinpath("data", filename)
    data = np.loadtxt(path, ndmin=2)
    wavelength_nm = data[:, 0]
    sigma_cm2 = data[:, 1]
    order = np.argsort(wavelength_nm)
    return wavelength_nm[order], sigma_cm2[order]


def photoionization_cross_section(wavelength_nm, species: str = "H"):
    """Interpolate a photoionization cross section.

    Parameters
    ----------
    wavelength_nm
        Wavelength in nm. Scalar and array inputs are accepted.
    species
        Neutral target species. Supported species include ``"H"``, ``"O"``,
        ``"CO"``, ``"CO2"``, ``"H2O"``, ``"N2"``, ``"Na"``, ``"Mg"``,
        and ``"Ca"``.

    Returns
    -------
    sigma_cm2
        Photoionization cross section in cm^2. Values outside the tabulated
        wavelength range are set to zero.
    """
    wavelength_arr = np.asarray(wavelength_nm, dtype=float)
    scalar_input = wavelength_arr.ndim == 0
    grid_nm, sigma_cm2 = _load_cross_section(_cross_section_file(species))
    out = np.interp(wavelength_arr, grid_nm, sigma_cm2, left=0.0, right=0.0)
    if scalar_input:
        return float(np.asarray(out))
    return out


ph_cross_section = photoionization_cross_section
cross_section = photoionization_cross_section

__all__ = ["photoionization_cross_section", "ph_cross_section", "cross_section"]
