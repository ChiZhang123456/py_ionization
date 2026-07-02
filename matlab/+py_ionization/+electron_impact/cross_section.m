function sigma_cm2 = cross_section(energy_eV, species)
%CROSS_SECTION Interpolate an electron impact ionization cross section.
%
% sigma_cm2 = py_ionization.electron_impact.cross_section(energy_eV)
% sigma_cm2 = py_ionization.electron_impact.cross_section(energy_eV, species)
%
% Inputs
%   energy_eV  Electron energy in eV.
%   species    Neutral target species, for example H, O, O2, CO2, N2, Na.
%
% Output
%   sigma_cm2  Cross section in cm^2. Values outside the tabulated energy
%              range are returned as NaN.

if nargin < 2
    species = "H";
end

filename = species_to_file(species);
path = py_ionization.data_path("electron_impact", filename);
data = read_numeric_table(path);
if size(data, 2) == 2
    column_index = 2;
else
    column_index = 3;
end
[grid_eV, sigma_table] = py_ionization.read_x1e16_table(path, column_index);

input_size = size(energy_eV);
energy_vec = double(energy_eV(:));
sigma_vec = interp1(grid_eV, sigma_table, energy_vec, "linear", nan);
sigma_cm2 = reshape(sigma_vec, input_size);
end

function filename = species_to_file(species)
key = lower(regexprep(string(species), "\s+", ""));

switch key
    case {"h", "h+"}
        filename = "electron_impact_H.txt";
    case {"c", "c+"}
        filename = "electron_impact_C.txt";
    case {"o", "o+"}
        filename = "electron_impact_O.txt";
    case {"o2", "o2+"}
        filename = "electron_impact_O2.txt";
    case {"co", "co+"}
        filename = "electron_impact_CO.txt";
    case {"n2", "n2+"}
        filename = "electron_impact_N2.txt";
    case {"co2", "co2+"}
        filename = "electron_impact_CO2.txt";
    case {"ar", "ar+"}
        filename = "electron_impact_Ar.txt";
    case {"na", "na+"}
        filename = "electron_impact_Na.txt";
    case {"mg", "mg+"}
        filename = "electron_impact_Mg.txt";
    case {"al", "al+"}
        filename = "electron_impact_Al.txt";
    case {"si", "si+"}
        filename = "electron_impact_Si.txt";
    case {"k", "k+"}
        filename = "electron_impact_K.txt";
    case {"ca", "ca+"}
        filename = "electron_impact_Ca.txt";
    case {"ti", "ti+"}
        filename = "electron_impact_Ti.txt";
    case {"fe", "fe+"}
        filename = "electron_impact_Fe.txt";
    otherwise
        error("py_ionization:electron_impact:UnsupportedSpecies", ...
            "Unsupported electron impact species: %s", species);
end
end

function data = read_numeric_table(path)
fid = fopen(path, "r");
if fid < 0
    error("py_ionization:electron_impact:OpenFailed", ...
        "Could not open file: %s", path);
end
cleanup = onCleanup(@() fclose(fid));
rows = {};
while true
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    line = strtrim(line);
    if line == "" || startsWith(line, "#")
        continue
    end
    values = sscanf(strrep(line, ",", " "), "%f").';
    if ~isempty(values)
        rows{end + 1, 1} = values; %#ok<AGROW>
    end
end
n_col = max(cellfun(@numel, rows));
data = nan(numel(rows), n_col);
for i = 1:numel(rows)
    values = rows{i};
    data(i, 1:numel(values)) = values;
end
end

