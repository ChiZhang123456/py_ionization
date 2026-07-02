"""Example: H production by H+ plus H charge exchange from MAVEN SWIA omni flux."""

import numpy as np
from py_space_zc import maven

from py_ionization.charge_exchange import charge_exchange_rate


tint = ["2024-08-01T00:00:00", "2024-08-01T01:00:00"]
swia = maven.get_data(tint, "swia_omni")
idef = swia["omni_flux"]

rate = charge_exchange_rate(
    idef.energy.data,
    idef.data,
    reaction_type="HpH>HpH",
)

print(f"Computed {rate.size} H charge exchange rates.")
print(f"Median rate: {float(np.nanmedian(rate)):.3e} s^-1")
