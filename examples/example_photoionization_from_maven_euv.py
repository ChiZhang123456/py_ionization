"""Example: H photoionization rate from MAVEN EUV L3 minute spectra.

Run from the parent project with the mars environment:

    C:\\Users\\Win\\.conda\\envs\\mars\\python.exe examples/example_photoionization_from_maven_euv.py
"""

from py_space_zc import maven

from py_ionization.photoionization import photoionization_rate


tint = ["2024-08-01T00:00:00", "2024-08-01T01:00:00"]
euv = maven.get_data(tint, "euv_l3")
rate = photoionization_rate(euv.wavelength.values, euv.values, species="H")

print(f"Computed {rate.size} H photoionization rates.")
print(f"Median rate: {float(__import__('numpy').nanmedian(rate)):.3e} s^-1")
