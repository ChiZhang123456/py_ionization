"""Example: H electron impact ionization rate from MAVEN SWEA omni spectra."""

import numpy as np
from py_space_zc import maven

from py_ionization.electron_impact import electron_impact_rate


tint = ["2024-08-01T00:00:00", "2024-08-01T01:00:00"]
swea = maven.get_data(tint, "swea_omni")
rate = electron_impact_rate(
    swea.energy.data,
    swea.data,
    species="H",
    theta_range=(-np.pi / 3, np.pi / 3),
    theta_type="elevation",
)

print(f"Computed {rate.size} H electron impact ionization rates.")
print(f"Median rate: {float(np.nanmedian(rate)):.3e} s^-1")
