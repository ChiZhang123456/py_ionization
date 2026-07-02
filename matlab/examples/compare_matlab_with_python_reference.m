% Compare MATLAB py_ionization results against Python py_ionization reference.

clear; close all;

this_file = mfilename("fullpath");
examples_dir = fileparts(this_file);
matlab_dir = fileparts(examples_dir);
repo_dir = fileparts(matlab_dir);
addpath(genpath(matlab_dir));

ref_path = fullfile(examples_dir, "python_reference_for_matlab_compare.mat");
if ~isfile(ref_path)
    error("Reference file not found. Run make_python_reference_for_matlab_compare.py first.");
end

ref = load(ref_path);

wavelength_nm = ref.wavelength_nm(:).';
energy_eV = ref.energy_eV(:).';
energy_rate_eV = ref.energy_rate_eV(:).';

ph_sigma_matlab = py_ionization.photoionization.cross_section(wavelength_nm, "H");
ei_sigma_matlab = py_ionization.electron_impact.cross_section(energy_eV, "H");
cex_sigma_matlab = py_ionization.charge_exchange.cross_section(energy_eV, "Ion", "H+", "Neutral", "H");

ph_rate_matlab = py_ionization.photoionization.production_rate(wavelength_nm, ref.irradiance_samples, "H");
ei_rate_matlab = py_ionization.electron_impact.production_rate(energy_rate_eV, ref.electron_def, "Species", "H");
cex_rate_matlab = py_ionization.charge_exchange.production_rate( ...
    energy_rate_eV, ref.ion_def, "ReactionType", "HpH>HpH");

summary = table( ...
    ["PH cross section"; "EI cross section"; "CEX cross section"; "PH production rate"; "EI production rate"; "CEX production rate"], ...
    [max_abs_diff(ph_sigma_matlab, ref.ph_sigma_python); ...
     max_abs_diff(ei_sigma_matlab, ref.ei_sigma_python); ...
     max_abs_diff(cex_sigma_matlab, ref.cex_sigma_python); ...
     max_abs_diff(ph_rate_matlab, ref.ph_rate_python); ...
     max_abs_diff(ei_rate_matlab, ref.ei_rate_python); ...
     max_abs_diff(cex_rate_matlab, ref.cex_rate_python)], ...
    [max_rel_diff(ph_sigma_matlab, ref.ph_sigma_python); ...
     max_rel_diff(ei_sigma_matlab, ref.ei_sigma_python); ...
     max_rel_diff(cex_sigma_matlab, ref.cex_sigma_python); ...
     max_rel_diff(ph_rate_matlab, ref.ph_rate_python); ...
     max_rel_diff(ei_rate_matlab, ref.ei_rate_python); ...
     max_rel_diff(cex_rate_matlab, ref.cex_rate_python)], ...
    'VariableNames', {'Quantity', 'MaxAbsDiff', 'MaxRelativeDiff'});

disp(summary);

set(groot, "defaultAxesFontName", "Arial");
set(groot, "defaultTextFontName", "Arial");
set(groot, "defaultLegendFontName", "Arial");

figure("Color", "w", "Position", [80, 80, 1500, 780]);

subplot(2, 3, 1);
plot(wavelength_nm, ref.ph_sigma_python, "-", "Color", [0.45 0.45 0.45], "LineWidth", 3); hold on;
plot(wavelength_nm, ph_sigma_matlab, "--", "Color", [0 0.45 0.74], "LineWidth", 2);
box on; grid on; xlabel("Wavelength (nm)"); ylabel("Cross section (cm^2)");
title("PH cross section"); legend("Python", "MATLAB", "Location", "best");

subplot(2, 3, 2);
loglog(energy_eV, ref.ei_sigma_python, "-", "Color", [0.45 0.45 0.45], "LineWidth", 3); hold on;
loglog(energy_eV, ei_sigma_matlab, "--", "Color", [0 0.45 0.74], "LineWidth", 2);
box on; grid on; xlabel("Energy (eV)"); ylabel("Cross section (cm^2)");
title("EI cross section"); legend("Python", "MATLAB", "Location", "best");

subplot(2, 3, 3);
loglog(energy_eV, ref.cex_sigma_python, "-", "Color", [0.45 0.45 0.45], "LineWidth", 3); hold on;
loglog(energy_eV, cex_sigma_matlab, "--", "Color", [0 0.45 0.74], "LineWidth", 2);
box on; grid on; xlabel("Energy (eV)"); ylabel("Cross section (cm^2)");
title("CEX cross section"); legend("Python", "MATLAB", "Location", "best");

sample_index = 0:(numel(ph_rate_matlab) - 1);

subplot(2, 3, 4);
plot(sample_index, ref.ph_rate_python(:), "-", "Color", [0.45 0.45 0.45], "LineWidth", 3); hold on;
plot(sample_index, ph_rate_matlab(:), "--", "Color", [0 0.45 0.74], "LineWidth", 2);
box on; grid on; xlabel("Synthetic sample"); ylabel("Rate (s^{-1})");
title("PH production rate"); legend("Python", "MATLAB", "Location", "best");

subplot(2, 3, 5);
plot(sample_index, ref.ei_rate_python(:), "-", "Color", [0.45 0.45 0.45], "LineWidth", 3); hold on;
plot(sample_index, ei_rate_matlab(:), "--", "Color", [0 0.45 0.74], "LineWidth", 2);
box on; grid on; xlabel("Synthetic sample"); ylabel("Rate (s^{-1})");
title("EI production rate"); legend("Python", "MATLAB", "Location", "best");

subplot(2, 3, 6);
plot(sample_index, ref.cex_rate_python(:), "-", "Color", [0.45 0.45 0.45], "LineWidth", 3); hold on;
plot(sample_index, cex_rate_matlab(:), "--", "Color", [0 0.45 0.74], "LineWidth", 2);
box on; grid on; xlabel("Synthetic sample"); ylabel("Rate (s^{-1})");
title("CEX production rate"); legend("Python", "MATLAB", "Location", "best");

out_dir = fullfile(repo_dir, "matlab", "comparison_outputs");
if ~isfolder(out_dir)
    mkdir(out_dir);
end
fig_path = fullfile(out_dir, "matlab_vs_python_py_ionization.png");
pdf_path = fullfile(out_dir, "matlab_vs_python_py_ionization.pdf");
csv_path = fullfile(out_dir, "matlab_vs_python_py_ionization_summary.csv");
exportgraphics(gcf, fig_path, "Resolution", 300);
exportgraphics(gcf, pdf_path, "ContentType", "vector");
writetable(summary, csv_path);

fprintf("Saved figure: %s\n", fig_path);
fprintf("Saved summary: %s\n", csv_path);

function value = max_abs_diff(a, b)
diff_value = abs(double(a(:)) - double(b(:)));
value = max(diff_value(isfinite(diff_value)));
end

function value = max_rel_diff(a, b)
a = double(a(:));
b = double(b(:));
mask = isfinite(a) & isfinite(b) & abs(b) > 0;
value = max(abs(a(mask) - b(mask)) ./ abs(b(mask)));
end

