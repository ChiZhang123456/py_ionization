# electron_impact

Functions:

- `electron_impact_cross_section(energy_eV, species="H")`
- `electron_impact_rate(energy_eV, electron_def, species="H")`

Files:

- `cross_section.py`: cross-section interpolation only.
- `production_rate.py`: rate calculation only.

Aliases:

- `ei_cross_section`
- `ei_prod_rate`
- `cross_section`
- `production_rate`

The rate calculation expects electron differential energy flux in
eV cm^-2 s^-1 sr^-1 eV^-1 and returns s^-1 for one neutral particle.
