function pass = t63_g3_gate(new_csv, ref_csv, mode, dims)
%T63_G3_GATE  T-63 reproduction gate: cell-by-cell TEXT comparison of two CSV files.
%
%  pass = t63_g3_gate(new_csv, ref_csv, mode, dims)
%
%  mode 'sweep'   (gate G3)
%      new_csv is a sweep file written by run_monte_carlo_paperC
%      (mc_snr_sweep_<stamp>.csv); ref_csv is the published sweep file.
%      The two header lines must be identical. Every cell of every column
%      except the two timing columns (timestamp, runtime_s) is compared as
%      text, row by row.
%  mode 'trials'  (comparison G3-T)
%      new_csv is a per-trial dump that carries the T-63 columns after the
%      earlier ones; ref_csv is a per-trial dump of the same block written
%      before T-63. The header of ref_csv must equal the leading columns of
%      the header of new_csv. Every cell of those leading columns is
%      compared as text, row by row.
%  dims
%      vector of positive integers whose product is the number of data rows
%      both files must hold (for example [13 5]: 13 sweep points, 5 method
%      rows each).
%
%  The function prints ONE verdict line and returns true for PASS. It
%  compares character strings only: no cell is converted to a number, and
%  no cell content is printed (a differing cell is named by its row and
%  column only). It does not stop the run; the caller decides what a FAIL
%  does. A CR before a line feed is ignored, so CRLF and LF files compare
%  alike.
%
%  T-63 (Phase D logging run). Base MATLAB only. 7-bit ASCII.

narginchk(4, 4);
new_csv = char(new_csv);
ref_csv = char(ref_csv);
mode    = char(mode);
assert(any(strcmp(mode, {'sweep', 'trials'})), ...
    't63_g3_gate: mode must be ''sweep'' or ''trials''.');
assert(isnumeric(dims) && ~isempty(dims) && all(dims(:) >= 1) && ...
    all(dims(:) == round(dims(:))), ...
    't63_g3_gate: dims must be a vector of positive integers.');

TIMING = {'timestamp', 'runtime_s'};

[hdr_new, rows_new] = loc_read(new_csv);
[hdr_ref, rows_ref] = loc_read(ref_csv);
n_exp = prod(dims(:));

if strcmp(mode, 'sweep')
    label  = 'G3';
    hdr_ok = isequal(hdr_new, hdr_ref);
    cols   = find(~ismember(hdr_ref, TIMING));
    domain = sprintf(['every cell of the %d columns of the reference ' ...
        'header other than %s, as text, row by row; the two header ' ...
        'lines must be identical'], numel(cols), strjoin(TIMING, ' and '));
else
    label  = 'G3-T';
    hdr_ok = (numel(hdr_new) >= numel(hdr_ref)) && ...
             isequal(hdr_new(1:numel(hdr_ref)), hdr_ref);
    cols   = 1:numel(hdr_ref);
    domain = sprintf(['every cell of the %d columns of the reference ' ...
        'header, as text, row by row; those columns must be the leading ' ...
        'columns of the new header'], numel(cols));
end
cols = cols(:).';

n_new   = numel(rows_new);
n_ref   = numel(rows_ref);
rows_ok = (n_new == n_exp) && (n_ref == n_exp);

n_shape = 0;
for i = 1:n_new
    if numel(rows_new{i}) ~= numel(hdr_new); n_shape = n_shape + 1; end
end
for i = 1:n_ref
    if numel(rows_ref{i}) ~= numel(hdr_ref); n_shape = n_shape + 1; end
end

n_cmp = 0;
n_bad = 0;
where = {};
if hdr_ok && rows_ok && (n_shape == 0)
    for i = 1:n_exp
        a = rows_new{i};
        b = rows_ref{i};
        for j = cols
            n_cmp = n_cmp + 1;
            if ~strcmp(a{j}, b{j})
                n_bad = n_bad + 1;
                if numel(where) < 10
                    where{end+1} = sprintf('row %d column %s', i, hdr_ref{j}); %#ok<AGROW>
                end
            end
        end
    end
end

pass = hdr_ok && rows_ok && (n_shape == 0) && ...
       (n_cmp == n_exp * numel(cols)) && (n_bad == 0);
if pass
    verdict = 'PASS';
else
    verdict = 'FAIL';
end

[~, nm_new, ex_new] = fileparts(new_csv);
[~, nm_ref, ex_ref] = fileparts(ref_csv);
fprintf(['%s [T-63; text-identity domain: %s; new file %s%s against ' ...
    'reference %s%s]: header condition met = %d | data rows: new %d, ' ...
    'reference %d, expected %d | rows whose field count differs from ' ...
    'their header = %d | columns compared %d | cells compared %d | ' ...
    'cells that differ %d: %s\n'], ...
    label, domain, nm_new, ex_new, nm_ref, ex_ref, ...
    double(hdr_ok), n_new, n_ref, n_exp, n_shape, numel(cols), n_cmp, ...
    n_bad, verdict);
for i = 1:numel(where)
    fprintf('  %s differs: %s\n', label, where{i});
end
if n_bad > numel(where)
    fprintf('  %s: %d further differing cells are not listed.\n', ...
        label, n_bad - numel(where));
end

end  % t63_g3_gate


% =========================================================================
%  LOCAL: loc_read -- header cells and data rows of a CSV file, as text
% =========================================================================
function [hdr, rows] = loc_read(f)
assert(isfile(f), 't63_g3_gate: file not found: %s', f);
txt = fileread(f);
txt = strrep(txt, char(13), '');      % drop CR: CRLF and LF compare alike
L   = strsplit(txt, char(10), 'CollapseDelimiters', false);
if ~isempty(L) && isempty(L{end})
    L(end) = [];                      % nothing follows the last line feed
end
assert(numel(L) >= 1, 't63_g3_gate: empty file: %s', f);
hdr  = strsplit(L{1}, ',', 'CollapseDelimiters', false);
rows = cell(numel(L) - 1, 1);
for i = 2:numel(L)
    rows{i-1} = strsplit(L{i}, ',', 'CollapseDelimiters', false);
end
end  % loc_read
