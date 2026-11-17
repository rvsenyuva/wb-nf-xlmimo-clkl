function pass = new5pl_gate(mode, label, new_csv, varargin)
%NEW5PL_GATE  Gates and counts of the pilot of O-3: TEXT comparisons of CSV files.
%
%  pass = new5pl_gate('cols',   label, new_csv, ref_csv, cols, n_exp, skip_method)
%      Both files must hold n_exp data rows and each name of the cell array
%      cols exactly once in their headers. Every cell of those columns is
%      compared as text, row i of one file with row i of the other. If
%      skip_method is not empty, both files must hold a column 'method';
%      rows whose method cell equals skip_method are skipped, and a row
%      whose method cells differ between the files counts as differing.
%  pass = new5pl_gate('cta',    label, new_csv, cta_csv, n_exp, scope)
%      Joins the n_exp rows of new_csv (a dump with the NEW5PL columns) on
%      SNR_dB and mc_seed with ALL rows of cta_csv (the CT-A per-trial file):
%      start_theta_deg against theta_hat_deg, start_r_m against r_hat_m, as
%      text; a row whose start_used is not 4, or a key twice in either file,
%      fails. scope 'full': cta_csv holds exactly n_exp rows; 'subset': at
%      least n_exp; 'self': the self-test of the join (loc_ctaself, below).
%  pass = new5pl_gate('seeds',  label, new_csv, snr_vec, seed_base, N_MC)
%      The dump must hold numel(snr_vec)*N_MC rows in the order the producer
%      writes them: for each sweep value, the trials 1 to N_MC. The cells
%      sweep_value, mc_seed_base, mc_idx and mc_seed are compared with the
%      text of %.6g and %d of the expected values.
%  pass = new5pl_gate('counts', label, new_csv, out_txt)
%      Counts text cells of the dump per sweep value and pooled and writes
%      them to the NEW file out_txt. No count is printed to the console.
%
%  Every mode prints ONE verdict line and returns true for PASS. The
%  function compares character strings only: no cell is converted to a
%  number, and no cell content is printed (a differing cell is named by its
%  row and column only). It does not stop the run; the caller decides what a
%  FAIL does. A CR before a line feed is ignored.
%
%  Model: diagnostics/t63_g3_gate.m (T-63). Pilot of O-3 (paperC-lab
%  sessions/2026-11-13_NEW5PL, SPEC_NEW5PL.md; mode cta: NEW5AM).
%
%  DRAFT. This file was written as text. It has not been parsed or run.
%  Base MATLAB only. 7-bit ASCII.

narginchk(4, 7);
mode    = char(mode);
label   = char(label);
new_csv = char(new_csv);

switch mode
    case 'cols'
        narginchk(7, 7);
        pass = loc_cols(label, new_csv, char(varargin{1}), varargin{2}, ...
                        varargin{3}, char(varargin{4}));
    case 'cta'
        narginchk(6, 6);
        pass = loc_cta(label, new_csv, char(varargin{1}), varargin{2}, char(varargin{3}));
    case 'seeds'
        narginchk(6, 6);
        pass = loc_seeds(label, new_csv, varargin{1}, varargin{2}, varargin{3});
    case 'counts'
        narginchk(4, 4);
        pass = loc_counts(label, new_csv, char(varargin{1}));
    otherwise
        error('new5pl_gate: unknown mode %s.', mode);
end

end  % new5pl_gate


% =========================================================================
%  LOCAL: loc_cols -- named columns of two files, row by row, as text
% =========================================================================
function pass = loc_cols(label, new_csv, ref_csv, cols, n_exp, skip_method)
[hn, rn] = loc_read(new_csv);
[hr, rr] = loc_read(ref_csv);
cols = cellstr(cols);
cols = cols(:).';
jn = loc_find(hn, cols);
jr = loc_find(hr, cols);
hdr_ok = all(jn > 0) && all(jr > 0);
mn = 0;
mr = 0;
if ~isempty(skip_method)
    mn = loc_find(hn, {'method'});
    mr = loc_find(hr, {'method'});
    hdr_ok = hdr_ok && (mn > 0) && (mr > 0);
end
rows_ok = (numel(rn) == n_exp) && (numel(rr) == n_exp);
n_shape = loc_shape(hn, rn) + loc_shape(hr, rr);

n_rows = 0;
n_skip = 0;
n_cmp  = 0;
n_bad  = 0;
where  = {};
if hdr_ok && rows_ok && (n_shape == 0)
    for i = 1:n_exp
        a = rn{i};
        b = rr{i};
        if ~isempty(skip_method)
            if ~strcmp(a{mn}, b{mr})
                n_bad = n_bad + 1;
                if numel(where) < 10
                    where{end+1} = sprintf('row %d column method', i); %#ok<AGROW>
                end
                continue
            end
            if strcmp(b{mr}, skip_method)
                n_skip = n_skip + 1;
                continue
            end
        end
        n_rows = n_rows + 1;
        for c = 1:numel(cols)
            n_cmp = n_cmp + 1;
            if ~strcmp(a{jn(c)}, b{jr(c)})
                n_bad = n_bad + 1;
                if numel(where) < 10
                    where{end+1} = sprintf('row %d column %s', i, cols{c}); %#ok<AGROW>
                end
            end
        end
    end
end

pass = hdr_ok && rows_ok && (n_shape == 0) && (n_cmp > 0) && ...
       (n_rows + n_skip == n_exp) && (n_cmp == n_rows * numel(cols)) && ...
       (n_bad == 0);

[~, nm_new, ex_new] = fileparts(new_csv);
[~, nm_ref, ex_ref] = fileparts(ref_csv);
fprintf(['%s [NEW5PL; text identity of cells, row by row; columns %s; ' ...
    'rows whose method is ''%s'' skipped; new file %s%s against reference ' ...
    '%s%s]: header condition met = %d | data rows: new %d, reference %d, ' ...
    'expected %d | rows whose field count differs from their header = %d ' ...
    '| rows compared %d, skipped %d | cells compared %d | cells that ' ...
    'differ %d: %s\n'], ...
    label, strjoin(cols, ' '), skip_method, nm_new, ex_new, nm_ref, ex_ref, ...
    double(hdr_ok), numel(rn), numel(rr), n_exp, n_shape, n_rows, n_skip, ...
    n_cmp, n_bad, loc_vs(pass));
for i = 1:numel(where)
    fprintf('  %s differs: %s\n', label, where{i});
end
if n_bad > numel(where)
    fprintf('  %s: %d further differing cells are not listed.\n', ...
        label, n_bad - numel(where));
end
end  % loc_cols


% =========================================================================
%  LOCAL: loc_cta -- the start of every row against the CT-A per-trial file
% =========================================================================
function pass = loc_cta(label, new_csv, cta_csv, n_exp, scope)
% NEW5AM: the map holds EVERY row of the reference. 'full': it holds exactly
% n_exp rows (SPEC_NEW5PL.md L802); 'subset': at least n_exp (S-0 part 2); a
% key twice in either file fails. 'self': loc_ctaself; new_csv a CT-A file.
[hc, rc] = loc_read(cta_csv);
[~, nm_ref, ex_ref] = fileparts(cta_csv);
if strcmp(scope, 'self')
    pass = loc_ctaself(label, new_csv, hc, rc, [nm_ref ex_ref], n_exp);
    return
end
[hn, rn] = loc_read(new_csv);
[~, nm_new, ex_new] = fileparts(new_csv);
k = loc_cta_join(hn, rn, hc, rc, n_exp, scope);
loc_cta_line(label, [nm_new ex_new], [nm_ref ex_ref], k);
pass = k.pass;
end  % loc_cta

function k = loc_cta_join(hn, rn, hc, rc, n_exp, scope)
% The counts of the join and its verdict k.pass (text only).
jn = loc_find(hn, {'SNR_dB', 'mc_seed', 'start_used', 'start_theta_deg', 'start_r_m'});
jc = loc_find(hc, {'SNR_dB', 'mc_seed', 'theta_hat_deg', 'r_hat_m'});
k = struct('scope', scope, 'hdr_ok', all(jn > 0) && all(jc > 0), 'n_new', numel(rn), ...
    'n_ref', numel(rc), 'n_exp', n_exp, 'n_shape', loc_shape(hn, rn) + loc_shape(hc, rc), ...
    'rows_ok', false, 'n_keys', 0, 'n_dup', 0, 'n_nokey', 0, 'n_nostart', 0, 'n_cmp', 0, 'n_bad', 0);
if strcmp(scope, 'full')
    k.rows_ok = (k.n_new == n_exp) && (k.n_ref == n_exp);
elseif strcmp(scope, 'subset')
    k.rows_ok = (k.n_new == n_exp) && (k.n_ref >= n_exp);
end
if k.hdr_ok && k.rows_ok && (k.n_shape == 0)
    ref  = containers.Map('KeyType', 'char', 'ValueType', 'any');
    seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    for i = 1:k.n_ref
        f = rc{i};
        ref([f{jc(1)} '|' f{jc(2)}]) = {f{jc(3)}, f{jc(4)}};
    end
    k.n_keys = double(ref.Count);
    for i = 1:n_exp
        a   = rn{i};
        key = [a{jn(1)} '|' a{jn(2)}];
        if isKey(seen, key)
            k.n_dup = k.n_dup + 1;
        elseif ~isKey(ref, key)
            k.n_nokey = k.n_nokey + 1;
        elseif ~strcmp(a{jn(3)}, '4')
            k.n_nostart = k.n_nostart + 1;
        else
            want    = ref(key);
            k.n_cmp = k.n_cmp + 2;
            k.n_bad = k.n_bad + double(~strcmp(a{jn(4)}, want{1})) + double(~strcmp(a{jn(5)}, want{2}));
        end
        seen(key) = true;
    end
end
k.pass = k.hdr_ok && k.rows_ok && (k.n_shape == 0) && (k.n_keys == k.n_ref) && (k.n_dup == 0) ...
    && (k.n_nokey == 0) && (k.n_nostart == 0) && (k.n_cmp == 2 * n_exp) && (k.n_bad == 0);
end  % loc_cta_join


% =========================================================================
%  LOCAL: loc_seeds -- sweep values, seed base, trial index and seed
% =========================================================================
function pass = loc_seeds(label, new_csv, snr_vec, seed_base, N_MC)
[hn, rn] = loc_read(new_csv);
jn = loc_find(hn, {'sweep_value', 'mc_seed_base', 'mc_idx', 'mc_seed'});
hdr_ok  = all(jn > 0);
n_exp   = numel(snr_vec) * N_MC;
rows_ok = (numel(rn) == n_exp);
n_shape = loc_shape(hn, rn);

n_chk = 0;
n_bad = 0;
if hdr_ok && rows_ok && (n_shape == 0)
    i = 0;
    for s = 1:numel(snr_vec)
        for m = 1:N_MC
            i = i + 1;
            a = rn{i};
            ok = strcmp(a{jn(1)}, sprintf('%.6g', snr_vec(s))) ...
              && strcmp(a{jn(2)}, sprintf('%d', seed_base)) ...
              && strcmp(a{jn(3)}, sprintf('%d', m)) ...
              && strcmp(a{jn(4)}, sprintf('%d', seed_base + m));
            n_chk = n_chk + 1;
            n_bad = n_bad + double(~ok);
        end
    end
end

pass = hdr_ok && rows_ok && (n_shape == 0) && (n_chk == n_exp) && (n_bad == 0);

[~, nm_new, ex_new] = fileparts(new_csv);
fprintf(['%s [NEW5PL; text of %%.6g and %%d; each sweep value holds the ' ...
    'trials 1 to N_MC with the seeds base+1 to base+N_MC, in order; file ' ...
    '%s%s]: header condition met = %d | data rows %d, expected %d | rows ' ...
    'whose field count differs from their header = %d | sweep values %d, ' ...
    'trials each %d, seed base %d | rows checked %d | rows that differ ' ...
    '%d: %s\n'], ...
    label, nm_new, ex_new, double(hdr_ok), numel(rn), n_exp, n_shape, ...
    numel(snr_vec), N_MC, seed_base, n_chk, n_bad, loc_vs(pass));
end  % loc_seeds


% =========================================================================
%  LOCAL: loc_counts -- counts of text cells, per sweep value and pooled
% =========================================================================
function pass = loc_counts(label, new_csv, out_txt)
[hn, rn] = loc_read(new_csv);
need = {'sweep_value', 'PhaseD_select', 'L_PhaseD_A', 'L_PhaseD_B', ...
        'L_PhaseD_C', 'converged', 'start_used', 'p_end', ...
        'on_theta_scan', 'on_u_scan', 'on_bpd_rgrid', ...
        'th_eq_start', 'r_eq_start'};
jn = loc_find(hn, need);
hdr_ok  = all(jn > 0);
n_shape = loc_shape(hn, rn);
is_new  = ~isfile(out_txt);
written = false;

NAMES = {'n', 'sel_A', 'sel_B', 'sel_C', 'sel_other', ...
         'start_bpd', 'start_psomp', 'start_midpoint', 'start_comparator', ...
         'th_eq_start', 'r_eq_start', ...
         'on_theta_scan', 'on_u_scan', 'on_bpd_rgrid', ...
         'L_absent', 'L_three_same_text', 'L_two_same_text', ...
         'p_end_zero_text', 'converged_1'};
NC = numel(NAMES);

if hdr_ok && (n_shape == 0) && ~isempty(rn) && is_new
    keys = {};
    C    = zeros(0, NC);
    for i = 1:numel(rn)
        a = rn{i};
        k = find(strcmp(keys, a{jn(1)}));
        if isempty(k)
            keys{end+1} = a{jn(1)}; %#ok<AGROW>
            C = [C; zeros(1, NC)]; %#ok<AGROW>
            k = numel(keys);
        end
        v = zeros(1, NC);
        v(1) = 1;
        sel  = a{jn(2)};
        v(2) = double(strcmp(sel, 'A'));
        v(3) = double(strcmp(sel, 'B'));
        v(4) = double(strcmp(sel, 'C'));
        v(5) = 1 - v(2) - v(3) - v(4);
        su   = a{jn(7)};
        v(6) = double(strcmp(su, '1'));
        v(7) = double(strcmp(su, '2'));
        v(8) = double(strcmp(su, '3'));
        v(9) = double(strcmp(su, '4'));
        v(10) = double(strcmp(a{jn(12)}, '1'));
        v(11) = double(strcmp(a{jn(13)}, '1'));
        v(12) = double(strcmp(a{jn(9)},  '1'));
        v(13) = double(strcmp(a{jn(10)}, '1'));
        v(14) = double(strcmp(a{jn(11)}, '1'));
        LA = a{jn(3)};
        LB = a{jn(4)};
        LC = a{jn(5)};
        absent = strcmp(LA, 'NaN') || strcmp(LB, 'NaN') || strcmp(LC, 'NaN');
        v(15) = double(absent);
        if ~absent
            e = double(strcmp(LA, LB)) + double(strcmp(LA, LC)) + double(strcmp(LB, LC));
            v(16) = double(e == 3);
            v(17) = double(e == 1);
        end
        v(18) = double(strcmp(a{jn(8)}, '0'));
        v(19) = double(strcmp(a{jn(6)}, '1'));
        C(k, :) = C(k, :) + v;
    end

    [~, nm_new, ex_new] = fileparts(new_csv);
    fid = fopen(out_txt, 'w');
    assert(fid > 0, 'new5pl_gate: cannot open %s', out_txt);
    fprintf(fid, ['NEW5PL-COUNTS %s: counts of text cells of %s%s, per ' ...
        'sweep value and pooled; no cell is converted to a number\n'], ...
        label, nm_new, ex_new);
    fprintf(fid, 'sweep_value,%s\n', strjoin(NAMES, ','));
    for k = 1:numel(keys)
        fprintf(fid, '%s', keys{k});
        fprintf(fid, ',%d', C(k, :));
        fprintf(fid, '\n');
    end
    fprintf(fid, 'pooled');
    fprintf(fid, ',%d', sum(C, 1));
    fprintf(fid, '\n');
    fclose(fid);
    written = true;
end

pass = hdr_ok && (n_shape == 0) && ~isempty(rn) && is_new && written;

[~, nm_new, ex_new] = fileparts(new_csv);
[~, nm_out, ex_out] = fileparts(out_txt);
fprintf(['%s [NEW5PL; counts of text cells per sweep value, written to a ' ...
    'file, none printed here; dump %s%s; counts file %s%s]: header ' ...
    'condition met = %d | data rows %d | rows whose field count differs ' ...
    'from their header = %d | counts file did not exist = %d | written = ' ...
    '%d: %s\n'], ...
    label, nm_new, ex_new, nm_out, ex_out, double(hdr_ok), numel(rn), ...
    n_shape, double(is_new), double(written), loc_vs(pass));
end  % loc_counts


% =========================================================================
%  LOCAL helpers
% =========================================================================
function j = loc_find(hdr, names)
% Position of each name in the header; 0 unless it occurs exactly once.
j = zeros(1, numel(names));
for c = 1:numel(names)
    k = find(strcmp(hdr, names{c}));
    if isscalar(k)
        j(c) = k;
    end
end
end  % loc_find

function n = loc_shape(hdr, rows)
% Number of rows whose field count differs from the header's.
n = 0;
for i = 1:numel(rows)
    if numel(rows{i}) ~= numel(hdr)
        n = n + 1;
    end
end
end  % loc_shape

function s = loc_vs(ok)
if ok
    s = 'PASS';
else
    s = 'FAIL';
end
end  % loc_vs

function [hdr, rows] = loc_read(f)
% Header cells and data rows of a CSV file, as text (t63_g3_gate.m L129-143).
assert(isfile(f), 'new5pl_gate: file not found: %s', f);
txt = fileread(f);
txt = strrep(txt, char(13), '');      % drop CR: CRLF and LF compare alike
L   = strsplit(txt, char(10), 'CollapseDelimiters', false);
if ~isempty(L) && isempty(L{end})
    L(end) = [];                      % nothing follows the last line feed
end
assert(numel(L) >= 1, 'new5pl_gate: empty file: %s', f);
hdr  = strsplit(L{1}, ',', 'CollapseDelimiters', false);
rows = cell(numel(L) - 1, 1);
for i = 2:numel(L)
    rows{i-1} = strsplit(L{i}, ',', 'CollapseDelimiters', false);
end
end  % loc_read


% =========================================================================
%  LOCAL: loc_cta_line -- the verdict line of mode cta (NEW5AM)
% =========================================================================
function loc_cta_line(label, f_new, f_ref, k)
% One verdict line (P-8); no cell content is printed.
fprintf(['%s [NEW5PL; text identity at the printed digits; start_theta_deg ' ...
    'against theta_hat_deg and start_r_m against r_hat_m, rows joined on ' ...
    'SNR_dB and mc_seed against every reference row; scope %s; new file %s ' ...
    'against reference %s]: header condition met = %d | data rows: new %d, ' ...
    'reference %d, expected %d | rows whose field count differs from their ' ...
    'header = %d | keys of the reference %d | rows whose key occurred before ' ...
    'in the new file %d | rows without a key in the reference %d | rows ' ...
    'whose start is not the comparator''s %d | cells compared %d | cells ' ...
    'that differ %d: %s\n'], ...
    label, k.scope, f_new, f_ref, double(k.hdr_ok), k.n_new, k.n_ref, k.n_exp, ...
    k.n_shape, k.n_keys, k.n_dup, k.n_nokey, k.n_nostart, k.n_cmp, k.n_bad, ...
    loc_vs(k.pass));
end  % loc_cta_line


% =========================================================================
%  LOCAL: loc_ctaself -- the self-test of mode cta (NEW5AM; G-SELF)
% =========================================================================
function pass = loc_ctaself(label, src_csv, hc, rc, f_ref, n_ref)
% No file is written. Four rows of src_csv (a CT-A per-trial file; the run
% sheet passes the reference itself), two from its start and two from its
% end, are copied as text into the columns of a dump with start_used 4 and
% joined, in scope 'subset', with the reference that loc_cta read. Case
% 'pass': the cells as copied; expected PASS, 8 cells compared, 0 differ.
% Case 'fail': the same rows, one character appended to the angle cell of
% the first; expected FAIL with 8 cells compared, 1 differing and every
% other count 0. A map built from a prefix of four reference rows would
% lack the two rows from the end, so case 'pass' passes only through a map
% formed from the whole reference and joined.
[hs, rs] = loc_read(src_csv);
js = loc_find(hs, {'SNR_dB', 'mc_seed', 'theta_hat_deg', 'r_hat_m'});
hn = {'SNR_dB', 'mc_seed', 'start_used', 'start_theta_deg', 'start_r_m'};
ok_src = all(js > 0) && (numel(rs) == n_ref) && (numel(rc) == n_ref) && (n_ref >= 8) ...
    && (loc_shape(hs, rs) == 0);
ok_p = false;
ok_f = false;
if ok_src
    pick = [1, 2, n_ref - 1, n_ref];
    rn   = cell(numel(pick), 1);
    for i = 1:numel(pick)
        f     = rs{pick(i)};
        rn{i} = {f{js(1)}, f{js(2)}, '4', f{js(3)}, f{js(4)}};
    end
    rf       = rn;
    rf{1}{4} = [rf{1}{4} '0'];
    kp = loc_cta_join(hn, rn, hc, rc, numel(pick), 'subset');
    kf = loc_cta_join(hn, rf, hc, rc, numel(pick), 'subset');
    loc_cta_line([label '-pass'], 'of four rows built in memory', f_ref, kp);
    loc_cta_line([label '-fail'], 'of four rows built in memory, one cell changed', f_ref, kf);
    ok_p = kp.pass && (kp.n_cmp == 8) && (kp.n_bad == 0) && (kp.n_keys == n_ref);
    ok_f = ~kf.pass && kf.hdr_ok && kf.rows_ok && (kf.n_shape == 0) && (kf.n_keys == n_ref) ...
        && (kf.n_dup == 0) && (kf.n_nokey == 0) && (kf.n_nostart == 0) ...
        && (kf.n_cmp == 8) && (kf.n_bad == 1);
end
pass = ok_src && ok_p && ok_f;
[~, nm_src, ex_src] = fileparts(src_csv);
fprintf(['%s [NEW5PL; self-test of mode cta: four rows of %s%s (rows 1, 2 and ' ...
    'the last two) joined in scope subset with every row of %s; case pass ' ...
    'expects PASS with cells compared 8 and cells that differ 0; case fail ' ...
    'expects FAIL with cells compared 8, cells that differ 1 and every other ' ...
    'count 0; no file written]: source and reference hold the expected rows ' ...
    '= %d | case pass as expected = %d | case fail as expected = %d: %s\n'], ...
    label, nm_src, ex_src, f_ref, double(ok_src), double(ok_p), double(ok_f), ...
    loc_vs(pass));
end  % loc_ctaself
