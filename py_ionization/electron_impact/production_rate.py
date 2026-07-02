from __future__ import annotations

import numpy as np

from py_ionization._utils import energy_bin_width, solid_angle_from_theta_range

from .cross_section import electron_impact_cross_section


def electron_impact_rate(
    energy_eV,
    electron_def,
    *,
    species: str = "H",
    spacecraft_potential_eV=None,
    dphi: float = 2 * np.pi,
    theta_range=(-np.pi / 4, np.pi / 4),
    theta_type: str = "auto",
    min_energy_eV: float = 10.0,
):
    """Compute an electron impact ionization production rate.

    Parameters
    ----------
    energy_eV
        1D electron energy bins in eV.
    electron_def
        Electron differential energy flux in eV cm^-2 s^-1 sr^-1 eV^-1.
    species
        Neutral target species for the ionization cross section.
    spacecraft_potential_eV
        Optional spacecraft potential in eV. A 1D array is interpreted as one
        value per time sample and is subtracted from ``energy_eV``.
    dphi, theta_range, theta_type
        Angular integration geometry.
    min_energy_eV
        Corrected energies less than or equal to this value are ignored.

    Returns
    -------
    rate_s_1
        Ionization rate in s^-1 for one neutral particle.

    Formula
    -------
    ``rate = integral (DEF / E) sigma(E) dE dOmega``.
    """
    energy = np.asarray(energy_eV, dtype=float)
    if energy.ndim != 1:
        raise ValueError("energy_eV must be a 1D array.")
    if np.any(energy <= 0):
        raise ValueError("energy_eV must be positive.")

    flux = np.asarray(electron_def, dtype=float)
    scalar_time = flux.ndim == 1
    if scalar_time:
        flux = flux[None, :]
    elif flux.ndim != 2:
        raise ValueError("electron_def must be 1D or 2D.")
    if flux.shape[1] != energy.size:
        raise ValueError("electron_def must have the same number of energy bins as energy_eV.")

    n_time = flux.shape[0]
    energy_mat = np.tile(energy[None, :], (n_time, 1))
    if spacecraft_potential_eV is None:
        potential = np.zeros_like(energy_mat)
    else:
        scpot = np.asarray(spacecraft_potential_eV, dtype=float)
        potential = scpot[:, None] if scpot.ndim == 1 else scpot
    corrected_energy = energy_mat - potential
    corrected_energy = np.where(corrected_energy <= min_energy_eV, np.nan, corrected_energy)

    sigma = electron_impact_cross_section(energy, species=species)
    dE = energy_bin_width(energy)[None, :]
    domega, _ = solid_angle_from_theta_range(theta_range, dphi=dphi, theta_type=theta_type)
    differential_particle_flux = flux / corrected_energy
    differential_particle_flux = np.where(np.isfinite(corrected_energy), differential_particle_flux, np.nan)
    rate = np.nansum(differential_particle_flux * sigma[None, :] * domega * dE, axis=1)

    if scalar_time:
        return float(rate[0])
    return rate


ei_prod_rate = electron_impact_rate
production_rate = electron_impact_rate

__all__ = ["electron_impact_rate", "ei_prod_rate", "production_rate"]
