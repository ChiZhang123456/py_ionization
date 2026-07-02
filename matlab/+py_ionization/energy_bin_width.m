function dE = energy_bin_width(energy_eV)
%ENERGY_BIN_WIDTH Return positive energy-bin widths in eV.

energy = double(energy_eV(:).');
if numel(energy) < 2
    error("py_ionization:energy_bin_width:TooFewSamples", ...
        "energy_eV must contain at least two bins.");
end

dE = abs(diff(energy));
dE = [dE, dE(end)];
end

