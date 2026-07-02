"""Compare py_ionization against py_space_zc.ionization.

This script checks that the standalone package reproduces the current
py_space_zc ionization calculations for representative synthetic inputs.
It writes a figure and a CSV summary under ``comparison_outputs``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from py_space_zc import ionization as old

from py_ionization.charge_exchange import charge_exchange_cross_section, charge_exchange_rate
from py_ionization.electron_impact import electron_impact_cross_section, electron_impact_rate
from py_ionization.photoionization import photoionization_cross_section, photoionization_rate


OUT_DIR = Path("comparison_outputs")
OUT_PNG = OUT_DIR / "py_ionization_vs_py_space_zc.png"
OUT_PDF = OUT_DIR / "py_ionization_vs_py_space_zc.pdf"
OUT_CSV = OUT_DIR / "py_ionization_vs_py_space_zc_summary.csv"


def _rel_err(new: np.ndarray, ref: np.ndarray) -> np.ndarray:
    new = np.asarray(new, dtype=float)
    ref = np.asarray(ref, dtype=float)
    return np.abs(new - ref) / np.maximum(np.abs(ref), 1e-300)


def _summary_row(name: str, new: np.ndarray, ref: np.ndarray) -> dict[str, float | str]:
    err = _rel_err(new, ref)
    finite = np.isfinite(err)
    return {
        "quantity": name,
        "n_finite": int(np.count_nonzero(finite)),
        "max_abs_diff": float(np.nanmax(np.abs(np.asarray(new) - np.asarray(ref)))),
        "max_relative_error": float(np.nanmax(err)),
        "median_relative_error": float(np.nanmedian(err)),
        "allclose": bool(np.allclose(new, ref, equal_nan=True)),
    }


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)

    energy = np.logspace(np.log10(12.0), np.log10(3500.0), 300)
    wavelength = np.linspace(1.0, 120.0, 300)

    ph_new = photoionization_cross_section(wavelength, "H")
    ph_ref = old.ph_cross_section(wavelength, "H")
    ei_new = electron_impact_cross_section(energy, "H")
    ei_ref = old.ei_cross_section(energy, "H")
    cex_new = charge_exchange_cross_section(energy, ion="H+", neutral="H")
    cex_ref = old.cex_cross_section(energy, ion="H+", neutral="H")

    n_time = 16
    ion_def = np.tile(np.exp(-energy / 900.0) * 1.5e7, (n_time, 1))
    electron_def = np.tile(np.exp(-energy / 180.0) * 2.0e7, (n_time, 1))
    irradiance = np.tile(np.exp(-wavelength / 35.0) * 1.0e-6, (n_time, 1))

    ph_rate_new = photoionization_rate(wavelength, irradiance, "H")
    ph_rate_ref = old.ph_prod_rate(wavelength, irradiance, "H")
    ei_rate_new = electron_impact_rate(energy, electron_def, species="H")
    ei_rate_ref = old.ei_prod_rate(energy, electron_def, species="H")
    cex_rate_new = charge_exchange_rate(energy, ion_def, reaction_type="HpH>HpH")
    cex_rate_ref = old.cex_prod_rate(energy, ion_def, reaction_type="HpH>HpH")

    rows = [
        _summary_row("photoionization_cross_section_H", ph_new, ph_ref),
        _summary_row("electron_impact_cross_section_H", ei_new, ei_ref),
        _summary_row("charge_exchange_cross_section_HpH", cex_new, cex_ref),
        _summary_row("photoionization_rate_H", ph_rate_new, ph_rate_ref),
        _summary_row("electron_impact_rate_H", ei_rate_new, ei_rate_ref),
        _summary_row("charge_exchange_rate_HpH", cex_rate_new, cex_rate_ref),
    ]
    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "dejavusans",
            "axes.unicode_minus": False,
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 8.5,
        }
    )
    fig, axes = plt.subplots(2, 3, figsize=(11.0, 6.2), constrained_layout=True)
    specs = [
        (axes[0, 0], wavelength, ph_new, ph_ref, "PH cross section", "Wavelength (nm)", "Cross section (cm$^2$)", "linear"),
        (axes[0, 1], energy, ei_new, ei_ref, "EI cross section", "Energy (eV)", "Cross section (cm$^2$)", "log"),
        (axes[0, 2], energy, cex_new, cex_ref, "CEX cross section", "Energy (eV)", "Cross section (cm$^2$)", "log"),
        (axes[1, 0], np.arange(n_time), ph_rate_new, ph_rate_ref, "PH production rate", "Synthetic sample", "Rate (s$^{-1}$)", "linear"),
        (axes[1, 1], np.arange(n_time), ei_rate_new, ei_rate_ref, "EI production rate", "Synthetic sample", "Rate (s$^{-1}$)", "linear"),
        (axes[1, 2], np.arange(n_time), cex_rate_new, cex_rate_ref, "CEX production rate", "Synthetic sample", "Rate (s$^{-1}$)", "linear"),
    ]

    for ax, x, new, ref, title, xlabel, ylabel, xscale in specs:
        ax.plot(x, ref, color="0.55", linewidth=2.5, label="py_space_zc")
        ax.plot(x, new, color="#1f77b4", linestyle="--", linewidth=1.5, label="py_ionization")
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.22)
        if xscale == "log":
            ax.set_xscale("log")
            ax.set_yscale("log")
        ax.legend(frameon=False)
        for spine in ax.spines.values():
            spine.set_visible(True)

    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    print(f"Saved figure: {OUT_PNG.resolve()}")
    print(f"Saved figure: {OUT_PDF.resolve()}")
    print(f"Saved summary: {OUT_CSV.resolve()}")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
