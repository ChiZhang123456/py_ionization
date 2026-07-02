import numpy as np

from py_ionization.charge_exchange import charge_exchange_cross_section, charge_exchange_rate
from py_ionization.electron_impact import electron_impact_cross_section, electron_impact_rate
from py_ionization.photoionization import photoionization_cross_section, photoionization_rate


def test_basic_shapes():
    energy = np.array([20.0, 100.0, 1000.0])
    wavelength = np.array([10.0, 20.0, 30.0, 40.0])

    assert photoionization_cross_section(wavelength, "H").shape == wavelength.shape
    assert electron_impact_cross_section(energy, "H").shape == energy.shape
    assert charge_exchange_cross_section(energy, ion="H+", neutral="H").shape == energy.shape


def test_basic_rates_are_finite():
    energy = np.array([20.0, 100.0, 1000.0])
    wavelength = np.array([10.0, 20.0, 30.0, 40.0])
    electron_def = np.ones((2, energy.size)) * 1e6
    ion_def = np.ones((2, energy.size)) * 1e6
    irradiance = np.ones((2, wavelength.size)) * 1e-6

    ph = photoionization_rate(wavelength, irradiance, "H")
    ei = electron_impact_rate(energy, electron_def, species="H")
    cex = charge_exchange_rate(energy, ion_def, reaction_type="HpH>HpH")

    assert np.all(np.isfinite(ph))
    assert np.all(np.isfinite(ei))
    assert np.all(np.isfinite(cex))
