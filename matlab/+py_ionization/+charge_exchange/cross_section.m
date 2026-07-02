function sigma_cm2 = cross_section(energy_eV, varargin)
%CROSS_SECTION Calculate a charge exchange cross section.
%
% sigma_cm2 = py_ionization.charge_exchange.cross_section(energy_eV)
% sigma_cm2 = py_ionization.charge_exchange.cross_section(energy_eV, "Ion", "H+", "Neutral", "H")
%
% Name-value options
%   Ion      Incident ion species, default "H+".
%   Neutral  Neutral target species, default "H".
%
% Output
%   sigma_cm2  Charge exchange cross section in cm^2.

p = inputParser;
p.addParameter("Ion", "H+", @(x) ischar(x) || isstring(x));
p.addParameter("Neutral", "H", @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

ion_key = clean_species(p.Results.Ion);
neutral_key = clean_species(p.Results.Neutral);

if ion_key == "h+" && any(neutral_key == ["na", "mg", "k", "ca", "al", "si", "ti", "other", "others"])
    filename = hp_metal_table(neutral_key);
    sigma_cm2 = table_cross_section(energy_eV, filename);
    return
end

if ion_key == "h+" && neutral_key == "o"
    key = "HpO>OpH";
elseif ion_key == "h+" && neutral_key == "h"
    key = "HpH>HpH";
elseif ion_key == "o+" && neutral_key == "h"
    key = "OpH>HpO";
elseif ion_key == "o+" && neutral_key == "o"
    key = "OpO>OpO";
elseif ion_key == "h+" && neutral_key == "n2"
    key = "HpN2>N2PH";
elseif ion_key == "h+" && neutral_key == "o2"
    key = "HpO2>O2PH";
else
    error("py_ionization:charge_exchange:UnsupportedReaction", ...
        "Unsupported charge exchange reaction for Ion=%s, Neutral=%s.", ...
        p.Results.Ion, p.Results.Neutral);
end

input_size = size(energy_eV);
energy = double(energy_eV(:));
log_ek = log(energy ./ 1e3);

switch key
    case "HpO>OpH"
        a1 = 2.91; a2 = 0.0886; a3 = 50.9;
        a4 = 4.73; a5 = -0.862; a6 = 0.0306;
        term1 = (a1 - a2 .* log_ek).^2 .* (1 - exp(-a3 ./ (energy ./ 1e3))).^2;
        term2 = (a4 - a5 .* log_ek) .* (1 - exp(-a6 ./ (energy ./ 1e3))).^2;
        result = term1 + term2;
    case "HpH>HpH"
        a1 = 4.15; a2 = 0.531; a3 = 67.3;
        result = (a1 - a2 .* log_ek).^2 .* (1 - exp(-a3 ./ (energy ./ 1e3))).^4.5;
    case "OpH>HpO"
        a1 = 3.13; a2 = 0.17; a3 = 87.5;
        result = (a1 - a2 .* log_ek).^2 .* (1 - exp(-a3 ./ (energy ./ 1e3))).^0.8;
    case "OpO>OpO"
        a1 = 4.07; a2 = 0.269; a3 = 415.0;
        result = (a1 - a2 .* log_ek).^2 .* (1 - exp(-a3 ./ (energy ./ 1e3))).^0.8;
    case "HpN2>N2PH"
        a1 = 12.5; a2 = 1.52; a3 = 3.97;
        a4 = 0.36; a5 = -1.2; a6 = 0.208; a7 = 0.741;
        term_a = a1 .* exp(-(log(energy ./ 1e3) - a2).^2 ./ a3) .* ...
            (1 - exp(-(energy ./ 1e3) ./ a4)).^2;
        term_b = (a5 - a6 .* log_ek).^2 .* (1 - exp(-a7 .* 1e3 ./ energy)).^2;
        result = term_a + term_b;
    case "HpO2>O2PH"
        a1 = 1.83; a2 = -0.545; a3 = 15.8;
        a4 = 6.35; a5 = -0.801; a6 = 0.24;
        term1 = (a1 - a2 .* log_ek).^2 .* (1 - exp(-a3 .* 1e3 ./ energy)).^1.5;
        term2 = (a4 - a5 .* log_ek).^2 .* (1 - exp(-a6 .* 1e3 ./ energy));
        result = term1 + term2;
end

sigma_cm2 = reshape(result .* 1e-16, input_size);
end

function key = clean_species(value)
key = lower(regexprep(string(value), "\s+", ""));
end

function filename = hp_metal_table(neutral_key)
switch neutral_key
    case "na"
        filename = "charge_exchange_Hp_Na.txt";
    case "mg"
        filename = "charge_exchange_Hp_Mg.txt";
    case "k"
        filename = "charge_exchange_Hp_K.txt";
    case "ca"
        filename = "charge_exchange_Hp_Ca.txt";
    otherwise
        filename = "charge_exchange_Hp_Others.txt";
end
end

function sigma_cm2 = table_cross_section(energy_eV, filename)
path = py_ionization.data_path("charge_exchange", filename);
[grid_eV, sigma_table] = py_ionization.read_x1e16_table(path, 2);
input_size = size(energy_eV);
energy_vec = double(energy_eV(:));
sigma_vec = interp1(grid_eV, sigma_table, energy_vec, "linear", nan);
sigma_cm2 = reshape(sigma_vec, input_size);
end

