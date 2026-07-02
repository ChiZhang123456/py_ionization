# MATLAB interface for py_ionization

This folder provides MATLAB functions for the same ionization calculations as
the Python package.

Add the MATLAB folder to your path:

```matlab
addpath(genpath('D:\Work_Work\Mars\MAVEN\Oxygen_ICW\py_ionization\matlab'))
```

## Photoionization

```matlab
wavelength_nm = linspace(1, 120, 300);
sigma = py_ionization.photoionization.cross_section(wavelength_nm, "H");
rate = py_ionization.photoionization.production_rate(wavelength_nm, irradiance_w_m2_nm, "H");
```

Aliases:

```matlab
sigma = py_ionization.photoionization.ph_cross_section(wavelength_nm, "H");
rate = py_ionization.photoionization.ph_prod_rate(wavelength_nm, irradiance_w_m2_nm, "H");
```

## Electron impact ionization

```matlab
energy_eV = logspace(1, 4, 200);
sigma = py_ionization.electron_impact.cross_section(energy_eV, "H");
rate = py_ionization.electron_impact.production_rate(energy_eV, electron_def, "Species", "H");
```

`electron_def` should be differential energy flux in
`eV cm^-2 s^-1 sr^-1 eV^-1`.

## Charge exchange

```matlab
energy_eV = logspace(1, 4, 200);
sigma = py_ionization.charge_exchange.cross_section(energy_eV, "Ion", "H+", "Neutral", "H");
rate = py_ionization.charge_exchange.production_rate(energy_eV, ion_def, "ReactionType", "HpH>HpH");
```

`ion_def` should be differential energy flux in
`eV cm^-2 s^-1 sr^-1 eV^-1`.

## Examples

- `examples/example_plot_cross_sections.m`
- `examples/compare_matlab_with_python_reference.m`

