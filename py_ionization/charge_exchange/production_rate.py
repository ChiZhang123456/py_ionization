from __future__ import annotations

import numpy as np

from py_ionization._utils import energy_bin_width, solid_angle_from_theta_range

from .cross_section import charge_exchange_cross_section


_REACTION_MAP = {
    "hpo>oph": ("H+", "O"),
    "hph>hph": ("H+", "H"),
    "oph>hpo": ("O+", "H"),
    "opo>opo": ("O+", "O"),
    "hpn2>n2ph": ("H+", "N2"),
    "hpo2>o2ph": ("H+", "O2"),
}


def _reaction_to_species(reaction_type: str) -> tuple[str, str]:
    key = str(reaction_type).lower().strip().replace(" ", "")
    try:
        return _REACTION_MAP[key]
    except KeyError as exc:
        supported = ", ".join(sorted(_REACTION_MAP))
        raise ValueError(f"Unsupported reaction_type {reaction_type!r}. Supported keys: {supported}") from exc


def charge_exchange_rate(
    energy_eV,
    ion_def,
    *,
    reaction_type: str | None = None,
    ion: str | None = None,
    neutral: str | None = None,
    dphi: float = 2 * np.pi,
    theta_range=(-np.pi / 4, np.pi / 4),
    theta_type: str = "auto",
    max_energy_eV: float | None = 4000.0,
):
    """Compute a charge exchange production rate from ion differential energy flux.

    Parameters
    ----------
    energy_eV
        1D ion energy bins in eV.
    ion_def
        Ion differential energy flux in eV cm^-2 s^-1 sr^-1 eV^-1.
    reaction_type
        Optional compact reaction key, for example ``"HpH>HpH"``. If this is
        supplied, ``ion`` and ``neutral`` are inferred unless explicitly set.
    ion, neutral
        Incident ion and neutral target species.
    dphi, theta_range, theta_type
        Angular integration geometry.
    max_energy_eV
        Energies greater than or equal to this value are ignored. The default
        is 4000 eV.

    Returns
    -------
    rate_s_1
        Charge exchange rate in s^-1 for one neutral target particle.

    Formula
    -------
    ``rate = integral (DEF / E) sigma(E) dE dOmega``.
    """
    energy = np.asarray(energy_eV, dtype=float)
    if energy.ndim != 1:
        raise ValueError("energy_eV must be a 1D array.")
    if np.any(energy <= 0):
        raise ValueError("energy_eV must be positive.")

    flux = np.asarray(ion_def, dtype=float)
    scalar_time = flux.ndim == 1
    if scalar_time:
        flux = flux[None, :]
    elif flux.ndim != 2:
        raise ValueError("ion_def must be 1D or 2D.")
    if flux.shape[1] != energy.size:
        raise ValueError("ion_def must have the same number of energy bins as energy_eV.")

    if reaction_type is not None:
        inferred_ion, inferred_neutral = _reaction_to_species(reaction_type)
        ion = inferred_ion if ion is None else ion
        neutral = inferred_neutral if neutral is None else neutral
    if ion is None or neutral is None:
        raise ValueError("Specify either reaction_type or both ion and neutral.")

    dE = energy_bin_width(energy)[None, :]
    energy_mat = energy[None, :]
    domega, _ = solid_angle_from_theta_range(theta_range, dphi=dphi, theta_type=theta_type)
    differential_particle_flux = flux / energy_mat
    if max_energy_eV is not None:
        differential_particle_flux = np.where(energy_mat >= max_energy_eV, np.nan, differential_particle_flux)

    sigma = np.asarray(charge_exchange_cross_section(energy, ion=ion, neutral=neutral), dtype=float)
    rate = np.nansum(differential_particle_flux * sigma[None, :] * domega * dE, axis=1)
    if scalar_time:
        return float(rate[0])
    return rate


cex_prod_rate = charge_exchange_rate
charge_exchange_prod_rate = charge_exchange_rate
production_rate = charge_exchange_rate

__all__ = ["charge_exchange_rate", "cex_prod_rate", "charge_exchange_prod_rate", "production_rate"]
