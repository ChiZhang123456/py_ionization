# photoionization

Functions:

- `photoionization_cross_section(wavelength_nm, species="H")`
- `photoionization_rate(wavelength_nm, irradiance_w_m2_nm, species="H")`

Files:

- `cross_section.py`: cross-section interpolation only.
- `production_rate.py`: rate calculation only.

Aliases:

- `ph_cross_section`
- `ph_prod_rate`
- `cross_section`
- `production_rate`

The rate calculation expects spectral irradiance in W m^-2 nm^-1 and returns
s^-1 for one neutral particle.
