from __future__ import annotations

from functools import lru_cache
from importlib.resources import files

import numpy as np

from py_ionization._utils import read_x1e16_table


_HP_METAL_TABLES = {
    "na": "charge_exchange_Hp_Na.txt",
    "mg": "charge_exchange_Hp_Mg.txt",
    "k": "charge_exchange_Hp_K.txt",
    "ca": "charge_exchange_Hp_Ca.txt",
    "al": "charge_exchange_Hp_Others.txt",
    "si": "charge_exchange_Hp_Others.txt",
    "ti": "charge_exchange_Hp_Others.txt",
    "other": "charge_exchange_Hp_Others.txt",
    "others": "charge_exchange_Hp_Others.txt",
}


def _clean_species(value: str) -> str:
    return str(value).lower().strip().replace(" ", "")


@lru_cache(maxsize=None)
def _load_table(filename: str) -> tuple[np.ndarray, np.ndarray]:
    path = files("py_ionization.charge_exchange").joinpath("data", filename)
    return read_x1e16_table(path, column=1)


def _table_cross_section(energy_eV, filename: str):
    energy = np.asarray(energy_eV, dtype=float)
    scalar_input = energy.ndim == 0
    grid, sigma = _load_table(filename)
    out = np.interp(energy, grid, sigma, left=np.nan, right=np.nan)
    if scalar_input:
        return float(np.asarray(out))
    return out


def charge_exchange_cross_section(energy_eV, *, ion: str = "H+", neutral: str = "H"):
    """Calculate a charge exchange cross section.

    Parameters
    ----------
    energy_eV
        Incident ion kinetic energy in eV.
    ion
        Incident ion species, for example ``"H+"`` or ``"O+"``.
    neutral
        Neutral target species, for example ``"H"``, ``"O"``, ``"N2"``,
        ``"O2"``, ``"Na"``, or ``"Mg"``.

    Returns
    -------
    sigma_cm2
        Charge exchange cross section in cm^2.
    """
    ion_key = _clean_species(ion)
    neutral_key = _clean_species(neutral)

    if ion_key == "h+" and neutral_key in _HP_METAL_TABLES:
        return _table_cross_section(energy_eV, _HP_METAL_TABLES[neutral_key])

    if ion_key == "h+" and neutral_key == "o":
        key = "HpO>OpH"
    elif ion_key == "h+" and neutral_key == "h":
        key = "HpH>HpH"
    elif ion_key == "o+" and neutral_key == "h":
        key = "OpH>HpO"
    elif ion_key == "o+" and neutral_key == "o":
        key = "OpO>OpO"
    elif ion_key == "h+" and neutral_key == "n2":
        key = "HpN2>N2PH"
    elif ion_key == "h+" and neutral_key == "o2":
        key = "HpO2>O2PH"
    else:
        raise ValueError(f"Unsupported charge exchange reaction for ion={ion!r}, neutral={neutral!r}.")

    energy = np.asarray(energy_eV, dtype=float)
    scalar_input = energy.ndim == 0
    log_ek = np.log(energy / 1e3)

    if key == "HpO>OpH":
        a1, a2, a3 = 2.91, 0.0886, 50.9
        a4, a5, a6 = 4.73, -0.862, 0.0306
        term1 = (a1 - a2 * log_ek) ** 2 * (1 - np.exp(-a3 / (energy / 1e3))) ** 2
        term2 = (a4 - a5 * log_ek) * (1 - np.exp(-a6 / (energy / 1e3))) ** 2
        result = term1 + term2
    elif key == "HpH>HpH":
        a1, a2, a3 = 4.15, 0.531, 67.3
        result = (a1 - a2 * log_ek) ** 2 * (1 - np.exp(-a3 / (energy / 1e3))) ** 4.5
    elif key == "OpH>HpO":
        a1, a2, a3 = 3.13, 0.17, 87.5
        result = (a1 - a2 * log_ek) ** 2 * (1 - np.exp(-a3 / (energy / 1e3))) ** 0.8
    elif key == "OpO>OpO":
        a1, a2, a3 = 4.07, 0.269, 415.0
        result = (a1 - a2 * log_ek) ** 2 * (1 - np.exp(-a3 / (energy / 1e3))) ** 0.8
    elif key == "HpN2>N2PH":
        a1, a2, a3 = 12.5, 1.52, 3.97
        a4, a5, a6, a7 = 0.36, -1.2, 0.208, 0.741
        term_a = a1 * np.exp(-(np.log(energy / 1e3) - a2) ** 2 / a3) * (1 - np.exp(-(energy / 1e3) / a4)) ** 2
        term_b = (a5 - a6 * log_ek) ** 2 * (1 - np.exp(-a7 * 1e3 / energy)) ** 2
        result = term_a + term_b
    elif key == "HpO2>O2PH":
        a1, a2, a3 = 1.83, -0.545, 15.8
        a4, a5, a6 = 6.35, -0.801, 0.24
        term1 = (a1 - a2 * log_ek) ** 2 * (1 - np.exp(-a3 * 1e3 / energy)) ** 1.5
        term2 = (a4 - a5 * log_ek) ** 2 * (1 - np.exp(-a6 * 1e3 / energy))
        result = term1 + term2
    else:
        raise RuntimeError("Unexpected charge exchange key.")

    sigma_cm2 = result * 1e-16
    if scalar_input:
        return float(np.asarray(sigma_cm2))
    return sigma_cm2


cex_cross_section = charge_exchange_cross_section
cross_section = charge_exchange_cross_section

__all__ = ["charge_exchange_cross_section", "cex_cross_section", "cross_section"]
