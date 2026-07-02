function [domega, theta_weight] = solid_angle_from_theta_range(theta_range, varargin)
%SOLID_ANGLE_FROM_THETA_RANGE Compute solid angle from theta limits.
%
% [domega, theta_weight] = py_ionization.solid_angle_from_theta_range(...)
%
% Name-value options:
%   Dphi      Azimuthal width in radians, default 2*pi.
%   ThetaType "auto", "elevation", or "polar", default "auto".

p = inputParser;
p.addParameter("Dphi", 2*pi, @(x) isnumeric(x) && isscalar(x));
p.addParameter("ThetaType", "auto", @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

theta = double(theta_range(:).');
if numel(theta) ~= 2
    error("py_ionization:solid_angle:BadThetaRange", ...
        "theta_range must contain exactly two limits.");
end

theta1 = theta(1);
theta2 = theta(2);
theta_min = min(theta1, theta2);
theta_max = max(theta1, theta2);
kind = lower(string(p.Results.ThetaType));

if kind == "auto" || kind == ""
    if theta_min < 0
        kind = "elevation";
    else
        kind = "polar";
    end
end

if any(kind == ["elevation", "elev", "latitude", "lat"])
    if theta_min < -pi/2 || theta_max > pi/2
        error("py_ionization:solid_angle:BadElevation", ...
            "Elevation theta_range must be within [-pi/2, pi/2].");
    end
    theta_factor = abs(sin(theta2) - sin(theta1));
    theta_weight = "cos";
elseif any(kind == ["polar", "colatitude", "colat"])
    if theta_min < 0 || theta_max > pi
        error("py_ionization:solid_angle:BadPolar", ...
            "Polar theta_range must be within [0, pi].");
    end
    theta_factor = abs(cos(theta1) - cos(theta2));
    theta_weight = "sin";
else
    error("py_ionization:solid_angle:BadThetaType", ...
        "ThetaType must be auto, elevation, or polar.");
end

domega = double(p.Results.Dphi) .* theta_factor;
end

