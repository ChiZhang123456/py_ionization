function rate_s_1 = production_rate(energy_eV, electron_def, varargin)
%PRODUCTION_RATE Compute an electron impact ionization production rate.
%
% rate_s_1 = py_ionization.electron_impact.production_rate(energy_eV, electron_def)
%
% Name-value options
%   Species                Neutral target species, default "H".
%   SpacecraftPotentialEV  Optional spacecraft potential in eV.
%   Dphi                   Azimuthal width in radians, default 2*pi.
%   ThetaRange             Theta limits in radians, default [-pi/4, pi/4].
%   ThetaType              "auto", "elevation", or "polar", default "auto".
%   MinEnergyEV            Corrected energies <= this value are ignored,
%                          default 10 eV.

p = inputParser;
p.addParameter("Species", "H", @(x) ischar(x) || isstring(x));
p.addParameter("SpacecraftPotentialEV", [], @(x) isnumeric(x) || isempty(x));
p.addParameter("Dphi", 2*pi, @(x) isnumeric(x) && isscalar(x));
p.addParameter("ThetaRange", [-pi/4, pi/4], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter("ThetaType", "auto", @(x) ischar(x) || isstring(x));
p.addParameter("MinEnergyEV", 10.0, @(x) isnumeric(x) && isscalar(x));
p.parse(varargin{:});

energy = double(energy_eV(:).');
if any(energy <= 0)
    error("py_ionization:electron_impact:BadEnergy", ...
        "energy_eV must be positive.");
end

flux = double(electron_def);
scalar_input = isvector(flux);
if scalar_input
    flux = reshape(flux, 1, []);
elseif size(flux, 2) ~= numel(energy) && size(flux, 1) == numel(energy)
    flux = flux.';
end
if size(flux, 2) ~= numel(energy)
    error("py_ionization:electron_impact:SizeMismatch", ...
        "electron_def must have the same number of energy bins as energy_eV.");
end

n_time = size(flux, 1);
energy_mat = repmat(energy, n_time, 1);
scpot = p.Results.SpacecraftPotentialEV;
if isempty(scpot)
    potential = zeros(size(energy_mat));
elseif isscalar(scpot)
    potential = repmat(double(scpot), size(energy_mat));
elseif isvector(scpot)
    potential = repmat(double(scpot(:)), 1, numel(energy));
else
    potential = double(scpot);
end

corrected_energy = energy_mat - potential;
corrected_energy(corrected_energy <= p.Results.MinEnergyEV) = nan;

sigma = py_ionization.electron_impact.cross_section(energy, p.Results.Species);
dE = py_ionization.energy_bin_width(energy);
[domega, ~] = py_ionization.solid_angle_from_theta_range( ...
    p.Results.ThetaRange, "Dphi", p.Results.Dphi, "ThetaType", p.Results.ThetaType);

differential_particle_flux = flux ./ corrected_energy;
differential_particle_flux(~isfinite(corrected_energy)) = nan;
integrand = differential_particle_flux .* sigma .* domega .* dE;
integrand(~isfinite(integrand)) = 0;
rate_s_1 = sum(integrand, 2);

if scalar_input
    rate_s_1 = rate_s_1(1);
end
end

