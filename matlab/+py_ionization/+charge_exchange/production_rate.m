function rate_s_1 = production_rate(energy_eV, ion_def, varargin)
%PRODUCTION_RATE Compute a charge exchange production rate.
%
% rate_s_1 = py_ionization.charge_exchange.production_rate(energy_eV, ion_def, "ReactionType", "HpH>HpH")
%
% Name-value options
%   ReactionType  Compact reaction key, for example "HpH>HpH".
%   Ion           Incident ion species.
%   Neutral       Neutral target species.
%   Dphi          Azimuthal width in radians, default 2*pi.
%   ThetaRange    Theta limits in radians, default [-pi/4, pi/4].
%   ThetaType     "auto", "elevation", or "polar", default "auto".
%   MaxEnergyEV   Energies >= this value are ignored, default 4000 eV.

p = inputParser;
p.addParameter("ReactionType", "", @(x) ischar(x) || isstring(x));
p.addParameter("Ion", "", @(x) ischar(x) || isstring(x));
p.addParameter("Neutral", "", @(x) ischar(x) || isstring(x));
p.addParameter("Dphi", 2*pi, @(x) isnumeric(x) && isscalar(x));
p.addParameter("ThetaRange", [-pi/4, pi/4], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter("ThetaType", "auto", @(x) ischar(x) || isstring(x));
p.addParameter("MaxEnergyEV", 4000.0, @(x) isnumeric(x) && (isempty(x) || isscalar(x)));
p.parse(varargin{:});

energy = double(energy_eV(:).');
if any(energy <= 0)
    error("py_ionization:charge_exchange:BadEnergy", ...
        "energy_eV must be positive.");
end

flux = double(ion_def);
scalar_input = isvector(flux);
if scalar_input
    flux = reshape(flux, 1, []);
elseif size(flux, 2) ~= numel(energy) && size(flux, 1) == numel(energy)
    flux = flux.';
end
if size(flux, 2) ~= numel(energy)
    error("py_ionization:charge_exchange:SizeMismatch", ...
        "ion_def must have the same number of energy bins as energy_eV.");
end

[ion, neutral] = infer_reaction_species(p.Results.ReactionType, p.Results.Ion, p.Results.Neutral);
if ion == "" || neutral == ""
    error("py_ionization:charge_exchange:MissingReaction", ...
        "Specify either ReactionType or both Ion and Neutral.");
end

dE = py_ionization.energy_bin_width(energy);
[domega, ~] = py_ionization.solid_angle_from_theta_range( ...
    p.Results.ThetaRange, "Dphi", p.Results.Dphi, "ThetaType", p.Results.ThetaType);

energy_mat = repmat(energy, size(flux, 1), 1);
differential_particle_flux = flux ./ energy_mat;
if ~isempty(p.Results.MaxEnergyEV)
    differential_particle_flux(energy_mat >= p.Results.MaxEnergyEV) = nan;
end

sigma = py_ionization.charge_exchange.cross_section(energy, "Ion", ion, "Neutral", neutral);
integrand = differential_particle_flux .* sigma .* domega .* dE;
integrand(~isfinite(integrand)) = 0;
rate_s_1 = sum(integrand, 2);

if scalar_input
    rate_s_1 = rate_s_1(1);
end
end

function [ion, neutral] = infer_reaction_species(reaction_type, ion, neutral)
ion = string(ion);
neutral = string(neutral);
key = lower(regexprep(string(reaction_type), "\s+", ""));

if key == ""
    return
end

switch key
    case "hpo>oph"
        inferred_ion = "H+"; inferred_neutral = "O";
    case "hph>hph"
        inferred_ion = "H+"; inferred_neutral = "H";
    case "oph>hpo"
        inferred_ion = "O+"; inferred_neutral = "H";
    case "opo>opo"
        inferred_ion = "O+"; inferred_neutral = "O";
    case "hpn2>n2ph"
        inferred_ion = "H+"; inferred_neutral = "N2";
    case "hpo2>o2ph"
        inferred_ion = "H+"; inferred_neutral = "O2";
    otherwise
        error("py_ionization:charge_exchange:UnsupportedReactionType", ...
            "Unsupported ReactionType: %s", reaction_type);
end

if ion == ""
    ion = inferred_ion;
end
if neutral == ""
    neutral = inferred_neutral;
end
end

