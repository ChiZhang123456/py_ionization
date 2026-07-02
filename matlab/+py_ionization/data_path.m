function path = data_path(process_name, filename)
%DATA_PATH Return the absolute path to a bundled cross-section table.
%
% path = py_ionization.data_path(process_name, filename)
%
% The MATLAB functions read the same data tables that are bundled with the
% Python package in ../py_ionization/<process>/data.

this_dir = fileparts(mfilename("fullpath"));
matlab_dir = fileparts(this_dir);
repo_dir = fileparts(matlab_dir);
path = fullfile(repo_dir, "py_ionization", process_name, "data", filename);

if ~isfile(path)
    error("py_ionization:data_path:FileNotFound", ...
        "Could not find data file: %s", path);
end
end

