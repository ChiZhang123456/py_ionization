function [x, y] = read_x1e16_table(path, column_index)
%READ_X1E16_TABLE Read a table whose cross section is in 1e-16 cm^2.
%
% [x, y] = py_ionization.read_x1e16_table(path, column_index)
%
% column_index is one-based, following MATLAB convention. The returned y is
% converted to cm^2.

if nargin < 2
    column_index = 2;
end

fid = fopen(path, "r");
if fid < 0
    error("py_ionization:read_x1e16_table:OpenFailed", ...
        "Could not open file: %s", path);
end

rows = {};
cleanup = onCleanup(@() fclose(fid));

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

if isempty(rows)
    error("py_ionization:read_x1e16_table:NoData", ...
        "No numeric data lines found in %s", path);
end

n_col = max(cellfun(@numel, rows));
data = nan(numel(rows), n_col);
for i = 1:numel(rows)
    values = rows{i};
    data(i, 1:numel(values)) = values;
end

if size(data, 2) < column_index
    error("py_ionization:read_x1e16_table:BadColumn", ...
        "Column %d is not available in %s", column_index, path);
end

[x, order] = sort(data(:, 1));
y = data(order, column_index) .* 1e-16;
end

