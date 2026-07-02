# charge_exchange

Functions:

- `charge_exchange_cross_section(energy_eV, ion="H+", neutral="H")`
- `charge_exchange_rate(energy_eV, ion_def, reaction_type="HpH>HpH")`

Files:

- `cross_section.py`: cross-section calculation or interpolation only.
- `production_rate.py`: rate calculation only.

Aliases:

- `cex_cross_section`
- `cex_prod_rate`
- `charge_exchange_prod_rate`
- `cross_section`
- `production_rate`

The rate calculation expects ion differential energy flux in
eV cm^-2 s^-1 sr^-1 eV^-1 and returns s^-1 for one neutral particle.
