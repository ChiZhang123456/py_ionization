# py_ionization

`py_ionization` is a small, self-contained Python package for ionization cross
sections and production-rate calculations. It is not tied to Mars. MAVEN is
used only in some examples because the original application was MAVEN solar
wind ionization analysis.

The package separates three processes:

- `py_ionization.photoionization`
- `py_ionization.electron_impact`
- `py_ionization.charge_exchange`

Each process folder separates the two scientific steps:

- `cross_section.py`: get or interpolate the relevant cross section.
- `production_rate.py`: integrate the cross section with the measured flux or
  irradiance to obtain a production rate.

All cross-section functions return cm^2. The production-rate functions return
s^-1 for one neutral target particle when the input differential energy fluxes
are in eV cm^-2 s^-1 sr^-1 eV^-1.

## Installation for local use

From this directory:

```powershell
C:\Users\Win\.conda\envs\mars\python.exe -m pip install -e .
```

or import directly by adding the repository path to `PYTHONPATH`.

## Photoionization

```python
from py_ionization.photoionization import (
    photoionization_cross_section,
    photoionization_rate,
)

sigma = photoionization_cross_section(wavelength_nm, species="H")
rate = photoionization_rate(wavelength_nm, irradiance_w_m2_nm, species="H")
```

Inputs:

- `wavelength_nm`: EUV wavelength in nm.
- `irradiance_w_m2_nm`: spectral irradiance in W m^-2 nm^-1.

Formula:

```text
rate = integral sigma(lambda) Phi(lambda) d lambda
Phi(lambda) = I(lambda) / (h c / lambda)
```

The cross section is internally converted from cm^2 to m^2, and wavelength is
converted from nm to m.

## Electron impact ionization

```python
from py_ionization.electron_impact import (
    electron_impact_cross_section,
    electron_impact_rate,
)

sigma = electron_impact_cross_section(energy_eV, species="H")
rate = electron_impact_rate(energy_eV, electron_def, species="H")
```

Inputs:

- `energy_eV`: electron energy bins in eV.
- `electron_def`: electron differential energy flux in eV cm^-2 s^-1 sr^-1 eV^-1.

Formula:

```text
rate = integral (DEF / E) sigma(E) dE dOmega
```

The default angular integration is full azimuth and an elevation range of
`[-pi/4, pi/4]`. For the PCW project, Chi often used `theta_range=(-pi/3, pi/3)`
or an equivalent `dtheta = 2*pi/3` assumption, depending on the script.

## Charge exchange

```python
from py_ionization.charge_exchange import (
    charge_exchange_cross_section,
    charge_exchange_rate,
)

sigma = charge_exchange_cross_section(energy_eV, ion="H+", neutral="H")
rate = charge_exchange_rate(
    energy_eV,
    ion_def,
    reaction_type="HpH>HpH",
)
```

Inputs:

- `energy_eV`: incident ion energy bins in eV.
- `ion_def`: ion differential energy flux in eV cm^-2 s^-1 sr^-1 eV^-1.

Formula:

```text
rate = integral (DEF / E) sigma(E) dE dOmega
```

The default ignores ion energies `>= 4000 eV`, following the original
MAVEN/SWIA charge-exchange workflow.

Supported compact reaction keys include:

- `HpH>HpH` for H+ + H charge exchange
- `HpO>OpH`
- `OpH>HpO`
- `OpO>OpO`
- `HpN2>N2PH`
- `HpO2>O2PH`

## Examples

See the `examples/` directory:

- `example_photoionization_from_maven_euv.py`
- `example_electron_impact_from_swea.py`
- `example_charge_exchange_from_swia.py`
- `example_plot_cross_sections.py`

## Notes for collaborators

The production-rate functions calculate rates per neutral particle. They do
not multiply by a neutral density model. If a volumetric production rate is
needed, multiply the returned rate by the neutral density in cm^-3 or m^-3
with consistent units.
