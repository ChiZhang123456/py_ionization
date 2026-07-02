"""Plot example cross sections for the three ionization processes."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from py_ionization.charge_exchange import charge_exchange_cross_section
from py_ionization.electron_impact import electron_impact_cross_section
from py_ionization.photoionization import photoionization_cross_section


out_dir = Path("example_outputs")
out_dir.mkdir(exist_ok=True)

energy = np.logspace(1, 4, 300)
wavelength = np.linspace(1.0, 120.0, 300)

fig, axes = plt.subplots(1, 3, figsize=(10, 3.2), constrained_layout=True)

axes[0].plot(wavelength, photoionization_cross_section(wavelength, "H"))
axes[0].set_title("Photoionization, H")
axes[0].set_xlabel("Wavelength (nm)")
axes[0].set_ylabel("Cross section (cm$^2$)")

axes[1].plot(energy, electron_impact_cross_section(energy, "H"))
axes[1].set_title("Electron impact, H")
axes[1].set_xlabel("Energy (eV)")
axes[1].set_xscale("log")
axes[1].set_yscale("log")

axes[2].plot(energy, charge_exchange_cross_section(energy, ion="H+", neutral="H"))
axes[2].set_title("Charge exchange, H+ + H")
axes[2].set_xlabel("Energy (eV)")
axes[2].set_xscale("log")
axes[2].set_yscale("log")

fig.savefig(out_dir / "cross_section_examples.png", dpi=200)
print(f"Saved {out_dir / 'cross_section_examples.png'}")
