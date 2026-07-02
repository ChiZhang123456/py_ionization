function rate_s_1 = production_rate(wavelength_nm, irradiance_w_m2_nm, species)
%PRODUCTION_RATE Calculate a photoionization rate from spectral irradiance.
%
% rate_s_1 = py_ionization.photoionization.production_rate(wavelength_nm, irradiance)
% rate_s_1 = py_ionization.photoionization.production_rate(wavelength_nm, irradiance, species)
%
% Inputs
%   wavelength_nm          1D wavelength grid in nm.
%   irradiance_w_m2_nm     Spectral irradiance in W m^-2 nm^-1. Rows are
%                          treated as time samples when a matrix is supplied.
%   species                Neutral target species.
%
% Output
%   rate_s_1               Photoionization rate in s^-1 for one neutral
%                          particle.

if nargin < 3
    species = "H";
end

h_j_s = 6.626e-34;
c_m_s = 3.0e8;

wavelength = double(wavelength_nm(:).');
if numel(wavelength) < 2
    error("py_ionization:photoionization:TooFewSamples", ...
        "wavelength_nm must contain at least two samples.");
end

irradiance = double(irradiance_w_m2_nm);
scalar_input = isvector(irradiance);
if scalar_input
    irradiance = reshape(irradiance, 1, []);
elseif size(irradiance, 2) ~= numel(wavelength) && size(irradiance, 1) == numel(wavelength)
    irradiance = irradiance.';
end

if size(irradiance, 2) ~= numel(wavelength)
    error("py_ionization:photoionization:SizeMismatch", ...
        "The last dimension of irradiance_w_m2_nm must match wavelength_nm.");
end

wavelength_m = wavelength .* 1e-9;
irradiance_w_m3 = irradiance .* 1e9;
photon_energy_j = h_j_s .* c_m_s ./ wavelength_m;
photon_flux_m3_s = irradiance_w_m3 ./ photon_energy_j;
dlambda_m = mean(diff(wavelength_m), "omitnan");
sigma_m2 = py_ionization.photoionization.cross_section(wavelength, species) .* 1e-4;

integrand = photon_flux_m3_s .* sigma_m2 .* dlambda_m;
integrand(~isfinite(integrand)) = 0;
rate_s_1 = sum(integrand, 2);

if scalar_input
    rate_s_1 = rate_s_1(1);
end
end

