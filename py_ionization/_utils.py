from __future__ import annotations

from pathlib import Path

import numpy as np


def read_x1e16_table(path: Path, *, column: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Read a text table whose cross-section column is in units of 1e-16 cm^2."""
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split(",") if "," in stripped else stripped.split()
            try:
                rows.append([float(part) for part in parts])
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"No numeric data lines found in {path}.")

    data = np.asarray(rows, dtype=float)
    if data.shape[1] <= column:
        raise ValueError(f"Column {column} is not available in {path}.")
    order = np.argsort(data[:, 0])
    return data[order, 0], data[order, column] * 1e-16


def energy_bin_width(energy: np.ndarray) -> np.ndarray:
    """Return positive energy-bin widths in eV for ascending or descending bins."""
    energy = np.asarray(energy, dtype=float)
    if energy.ndim != 1:
        raise ValueError("energy must be a 1D array.")
    if energy.size < 2:
        raise ValueError("energy must contain at least two bins.")
    dE = np.abs(np.diff(energy))
    return np.concatenate([dE, dE[-1:]])


def solid_angle_from_theta_range(
    theta_range: tuple[float, float] | list[float],
    *,
    dphi: float = 2 * np.pi,
    theta_type: str = "auto",
) -> tuple[float, str]:
    """Compute a solid angle from a theta interval and azimuthal width.

    Parameters
    ----------
    theta_range
        Two theta limits in radians.
    dphi
        Azimuthal width in radians.
    theta_type
        ``"elevation"`` means theta is measured from the XY plane and uses
        cos(theta) weighting. ``"polar"`` means theta is measured from +Z and
        uses sin(theta) weighting. ``"auto"`` treats negative limits as
        elevation and non-negative limits as polar.

    Returns
    -------
    domega, theta_weight
        Solid angle in sr and a text label, ``"cos"`` or ``"sin"``.
    """
    theta = np.asarray(theta_range, dtype=float)
    if theta.shape != (2,):
        raise ValueError("theta_range must contain exactly two limits.")

    theta1, theta2 = theta
    theta_min = min(theta1, theta2)
    theta_max = max(theta1, theta2)

    kind = str(theta_type).lower()
    if kind in {"auto", ""}:
        kind = "elevation" if theta_min < 0.0 else "polar"

    if kind in {"elevation", "elev", "latitude", "lat"}:
        if theta_min < -np.pi / 2 or theta_max > np.pi / 2:
            raise ValueError("Elevation theta_range must be within [-pi/2, pi/2].")
        theta_factor = abs(np.sin(theta2) - np.sin(theta1))
        theta_weight = "cos"
    elif kind in {"polar", "colatitude", "colat"}:
        if theta_min < 0.0 or theta_max > np.pi:
            raise ValueError("Polar theta_range must be within [0, pi].")
        theta_factor = abs(np.cos(theta1) - np.cos(theta2))
        theta_weight = "sin"
    else:
        raise ValueError("theta_type must be 'auto', 'elevation', or 'polar'.")

    return float(dphi) * theta_factor, theta_weight
