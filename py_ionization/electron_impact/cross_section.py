from __future__ import annotations

from functools import lru_cache
from importlib.resources import files

import numpy as np

from py_ionization._utils import read_x1e16_table


_SPECIES_TO_FILE = {
    "h": "electron_impact_H.txt",
    "h+": "electron_impact_H.txt",
    "c": "electron_impact_C.txt",
    "c+": "electron_impact_C.txt",
    "o": "electron_impact_O.txt",
    "o+": "electron_impact_O.txt",
    "o2": "electron_impact_O2.txt",
    "o2+": "electron_impact_O2.txt",
    "co": "electron_impact_CO.txt",
    "co+": "electron_impact_CO.txt",
    "n2": "electron_impact_N2.txt",
    "n2+": "electron_impact_N2.txt",
    "co2": "electron_impact_CO2.txt",
    "co2+": "electron_impact_CO2.txt",
    "ar": "electron_impact_Ar.txt",
    "ar+": "electron_impact_Ar.txt",
    "na": "electron_impact_Na.txt",
    "na+": "electron_impact_Na.txt",
    "mg": "electron_impact_Mg.txt",
    "mg+": "electron_impact_Mg.txt",
    "al": "electron_impact_Al.txt",
    "al+": "electron_impact_Al.txt",
    "si": "electron_impact_Si.txt",
    "si+": "electron_impact_Si.txt",
    "k": "electron_impact_K.txt",
    "k+": "electron_impact_K.txt",
    "ca": "electron_impact_Ca.txt",
    "ca+": "electron_impact_Ca.txt",
    "ti": "electron_impact_Ti.txt",
    "ti+": "electron_impact_Ti.txt",
    "fe": "electron_impact_Fe.txt",
    "fe+": "electron_impact_Fe.txt",
}


def _data_file(species: str) -> str:
    key = str(species).lower().strip().replace(" ", "")
    try:
        return _SPECIES_TO_FILE[key]
    except KeyError as exc:
        supported = ", ".join(sorted(set(_SPECIES_TO_FILE)))
        raise ValueError(f"Unsupported electron impact species {species!r}. Supported keys: {supported}") from exc


@lru_cache(maxsize=None)
def _load_cross_section(filename: str) -> tuple[np.ndarray, np.ndarray]:
    path = files("py_ionization.electron_impact").joinpath("data", filename)
    rows = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split(",") if "," in stripped else stripped.split()
            try:
                rows.append([float(part) for part in parts])
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"No numeric data lines found in {path}.")
    data = np.asarray(rows, dtype=float)
    column = 1 if data.shape[1] == 2 else 2
    return read_x1e16_table(path, column=column)


def electron_impact_cross_section(energy_eV, species: str = "H"):
    """Interpolate an electron impact ionization cross section.

    Parameters
    ----------
    energy_eV
        Electron energy in eV. Scalar and array inputs are accepted.
    species
        Neutral target species, for example ``"H"``, ``"O"``, ``"O2"``,
        ``"CO2"``, ``"N2"``, ``"Na"``, or ``"Mg"``.

    Returns
    -------
    sigma_cm2
        Ionization cross section in cm^2. Values outside the tabulated energy
        range are returned as NaN.
    """
    energy_arr = np.asarray(energy_eV, dtype=float)
    scalar_input = energy_arr.ndim == 0
    grid_eV, sigma_cm2 = _load_cross_section(_data_file(species))
    out = np.interp(energy_arr, grid_eV, sigma_cm2, left=np.nan, right=np.nan)
    if scalar_input:
        return float(np.asarray(out))
    return out


ei_cross_section = electron_impact_cross_section
cross_section = electron_impact_cross_section

__all__ = ["electron_impact_cross_section", "ei_cross_section", "cross_section"]
