% Plot cross sections from the MATLAB py_ionization interface.

clear; close all;

this_file = mfilename("fullpath");
examples_dir = fileparts(this_file);
matlab_dir = fileparts(examples_dir);
addpath(genpath(matlab_dir));

set(groot, "defaultAxesFontName", "Arial");
set(groot, "defaultTextFontName", "Arial");
set(groot, "defaultLegendFontName", "Arial");

wavelength_nm = linspace(0.1, 120, 300);
energy_eV = logspace(1, 3.5, 300);

ph_sigma = py_ionization.photoionization.cross_section(wavelength_nm, "H");
ei_sigma = py_ionization.electron_impact.cross_section(energy_eV, "H");
cex_sigma = py_ionization.charge_exchange.cross_section( ...
    energy_eV, "Ion", "H+", "Neutral", "H");

figure("Color", "w", "Position", [100, 100, 1200, 360]);

subplot(1, 3, 1);
plot(wavelength_nm, ph_sigma, "LineWidth", 2);
box on; grid on;
xlabel("Wavelength (nm)");
ylabel("Cross section (cm^2)");
title("Photoionization");

subplot(1, 3, 2);
loglog(energy_eV, ei_sigma, "LineWidth", 2);
box on; grid on;
xlabel("Energy (eV)");
ylabel("Cross section (cm^2)");
title("Electron impact");

subplot(1, 3, 3);
loglog(energy_eV, cex_sigma, "LineWidth", 2);
box on; grid on;
xlabel("Energy (eV)");
ylabel("Cross section (cm^2)");
title("Charge exchange");

out_dir = fullfile(fileparts(matlab_dir), "matlab", "comparison_outputs");
if ~isfolder(out_dir)
    mkdir(out_dir);
end
exportgraphics(gcf, fullfile(out_dir, "matlab_cross_sections.png"), "Resolution", 300);

