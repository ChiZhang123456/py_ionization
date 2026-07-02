function sigma_cm2 = cross_section(wavelength_nm, species)
%CROSS_SECTION Interpolate a photoionization cross section.
%
% sigma_cm2 = py_ionization.photoionization.cross_section(wavelength_nm)
% sigma_cm2 = py_ionization.photoionization.cross_section(wavelength_nm, species)
%
% Inputs
%   wavelength_nm  Wavelength in nm.
%   species        Neutral target species. Supported examples include H, O,
%                  CO, CO2, H2O, N2, Na, Mg, and Ca.
%
% Output
%   sigma_cm2      Photoionization cross section in cm^2. Values outside the
%                  tabulated wavelength range are set to zero.

if nargin < 2
    species = "H";
end

filename = species_to_file(species);
path = py_ionization.data_path("photoionization", filename);
data = readmatrix(path, "FileType", "text", "CommentStyle", "#");
data = data(all(isfinite(data(:, 1:2)), 2), 1:2);
[grid_nm, order] = sort(data(:, 1));
sigma_table = data(order, 2);

input_size = size(wavelength_nm);
wavelength_vec = double(wavelength_nm(:));
sigma_vec = interp1(grid_nm, sigma_table, wavelength_vec, "linear", 0);
sigma_cm2 = reshape(sigma_vec, input_size);
end

function filename = species_to_file(species)
key = lower(regexprep(string(species), "\s+", ""));

switch key
    case {"h", "p", "h+"}
        filename = "ph_H_cross_section.txt";
    case {"o", "o+"}
        filename = "ph_O_cross_section.txt";
    case {"co", "co+"}
        filename = "ph_CO_cross_section.txt";
    case {"co2", "co2+"}
        filename = "ph_CO2_cross_section.txt";
    case {"h2o", "h2o+"}
        filename = "ph_H2O_cross_section.txt";
    case {"n2", "n2+"}
        filename = "ph_N2_cross_section.txt";
    case {"na", "na+"}
        filename = "ph_Na_cross_section.txt";
    case {"mg", "mg+"}
        filename = "ph_Mg_cross_section.txt";
    case {"ca", "ca+"}
        filename = "ph_Ca_cross_section.txt";
    otherwise
        error("py_ionization:photoionization:UnsupportedSpecies", ...
            "Unsupported photoionization species: %s", species);
end
end

