"""Create a Python reference file for the MATLAB comparison example."""

from pathlib import Path

import numpy as np
from scipy.io import savemat

from py_ionization.photoionization import photoionization_cross_section, photoionization_rate
from py_ionization.electron_impact import electron_impact_cross_section, electron_impact_rate
from py_ionization.charge_exchange import charge_exchange_cross_section, charge_exchange_rate


def main():
    here = Path(__file__).resolve().parent
    out_path = here / "python_reference_for_matlab_compare.mat"

    wavelength_nm = np.linspace(0.1, 120.0, 300)
    energy_eV = np.logspace(1.0, 3.5, 300)
    sample_index = np.arange(16, dtype=float)
    energy_rate_eV = np.logspace(np.log10(20.0), np.log10(3000.0), 16)

    irradiance = 1e-6 * (1.0 + 0.1 * np.sin(wavelength_nm / 15.0))
    irradiance_samples = (1.0 + 0.02 * sample_index[:, None]) * irradiance[None, :]

    electron_def = (
        1e6
        * (1.0 + 0.03 * sample_index[:, None])
        * (energy_rate_eV[None, :] / 100.0) ** -0.7
    )
    ion_def = (
        1e6
        * (1.0 + 0.025 * sample_index[:, None])
        * (energy_rate_eV[None, :] / 100.0) ** -0.5
    )

    savemat(
        out_path,
        {
            "wavelength_nm": wavelength_nm,
            "energy_eV": energy_eV,
            "energy_rate_eV": energy_rate_eV,
            "irradiance_samples": irradiance_samples,
            "electron_def": electron_def,
            "ion_def": ion_def,
            "ph_sigma_python": photoionization_cross_section(wavelength_nm, "H"),
            "ei_sigma_python": electron_impact_cross_section(energy_eV, "H"),
            "cex_sigma_python": charge_exchange_cross_section(energy_eV, ion="H+", neutral="H"),
            "ph_rate_python": photoionization_rate(wavelength_nm, irradiance_samples, "H"),
            "ei_rate_python": electron_impact_rate(energy_rate_eV, electron_def, species="H"),
            "cex_rate_python": charge_exchange_rate(energy_rate_eV, ion_def, reaction_type="HpH>HpH"),
        },
    )
    print(out_path)


if __name__ == "__main__":
    main()

