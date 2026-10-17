function new4_plot_m1(results_dir, stamp)
%NEW4_PLOT_M1  M1 (fig:rangeRMSE_SNR, fig8a_rmse_r_snr.pdf): the range accuracy
%  of the four estimators beside two range bound curves of the NEW-4 run.
%
%  new4_plot_m1(results_dir, stamp)
%
%  Written by the NEW-4 code session of paperC-lab (sessions/2026-10-18_NEW4-CODE)
%  from NEW4_SPEC.md as decided in DECISIONS_N4G.md (N4G-D1, D10, D13, D15, D20).
%  Written as text; never parsed or executed by its author (N6).
%
%  Reads only: the committed CSV of M1 and M3 (the estimator columns, for the
%  curves) and the bound table new4_bound_m1_<stamp>.csv of new4_bounds.m in
%  results_dir. Computes no bound, calls no function of the code repository,
%  and forms, draws, prints and writes no quotient or difference of an estimator
%  statistic and a bound (inherited-finding block of Ledger N-195). Writes one
%  PDF into results_dir, named new4_<published PDF base name>_<stamp>.pdf (Y1;
%  N4G-D20: the author's block copies it under the published name later).
%
%  Labels in comments: Yn = requirement of the starter of the code session;
%  N4G-Dn = author decisions; OPEN-n = a choice no record fixes, listed in
%  OPEN_NEW4.md.

% =========================================================================
%  Y1: arguments and named constants; every path is relative to the root of
%  the code repository (D-20; N4G-D17)
% =========================================================================
if nargin ~= 2 || ~ischar(results_dir) || isempty(results_dir) || ~ischar(stamp) || isempty(stamp)
    error('new4_plot_m1:usage', 'new4_plot_m1(results_dir, stamp): two non-empty character vectors are required.');
end
here_dir = fileparts(mfilename('fullpath'));    % Y1: the file stands in plotting/ of the code repository
root_dir = fileparts(here_dir);
CSV_M1M3   = fullfile('results', 'regime_probe_rhi020', 'mc_snr_sweep_20260528_204613.csv');   % CR [R40]: the committed CSV of M1 and M3
TABLE_FILE = ['new4_bound_m1_' stamp '.csv'];                                                    % the bound table of M1 (O1)
PDF_BASE   = 'fig8a_rmse_r_snr';   % COPIED from L:plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L186
M_B1 = 'WB-BPD';             % COPIED from L:plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L53 (N4G-D13: names only)
M_B2 = 'WB-P-SOMP';       % COPIED from L:plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L54
M_B4 = 'WB-CL-KL';         % COPIED from L:plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L55
M_B5 = 'WB-DL-OMP';       % COPIED from L:plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L56
f_csv = fullfile(root_dir, CSV_M1M3);
f_tab = fullfile(results_dir, TABLE_FILE);

% =========================================================================
%  Y3; N4G-D1, N4G-D13: the estimator curves, re-drawn from the committed CSV
% =========================================================================
fprintf('new4_plot_m1: reading %s\n', f_csv);
[hdr_e, body_e] = loc_read_table(f_csv);
[snr_B1, rmse_B1] = loc_series(hdr_e, body_e, M_B1, 'SNR_dB', 'RMSE_r_m');
[snr_B2, rmse_B2] = loc_series(hdr_e, body_e, M_B2, 'SNR_dB', 'RMSE_r_m');
[snr_B4, rmse_B4] = loc_series(hdr_e, body_e, M_B4, 'SNR_dB', 'RMSE_r_m');
[snr_B5, rmse_B5] = loc_series(hdr_e, body_e, M_B5, 'SNR_dB', 'RMSE_r_m');
snr_all = snr_B4;

% Y2; N4G-D10: the published plot window, lines copied from the blob
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L67-71 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
win_mask = snr_B4 <= 17.5 + 1e-9;
snr_B1   = snr_B1(win_mask);   rmse_B1 = rmse_B1(win_mask);
snr_B2   = snr_B2(win_mask);   rmse_B2 = rmse_B2(win_mask);
snr_B4   = snr_B4(win_mask);   rmse_B4 = rmse_B4(win_mask);
snr_B5   = snr_B5(win_mask);   rmse_B5 = rmse_B5(win_mask);
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L74 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
n_pts = numel(snr_B4);

% =========================================================================
%  Y3; N4G-D1, D13; N-32 (a), (b): the two range bound curves, read from O1
% =========================================================================
fprintf('new4_plot_m1: reading %s\n', f_tab);
[hdr_b, body_b] = loc_read_table(f_tab);
ib = loc_col(hdr_b, 'range_bound.');
[snr_fa, y_fa, lab_fa] = loc_family_rows(hdr_b, body_b, 'full_array', ib);
[snr_tm, y_tm, lab_tm] = loc_family_rows(hdr_b, body_b, 'trial_matched', ib);
if numel(snr_fa) ~= numel(snr_all) || numel(snr_tm) ~= numel(snr_all) ...
        || any(abs(snr_fa - snr_all) > 1e-9) || any(abs(snr_tm - snr_all) > 1e-9)
    error('new4_plot_m1:grid', 'the SNR points of the bound table differ from those of the estimator rows.');
end
bound_fa = y_fa(win_mask);                      % the same window as the estimator curves (N4G-D10; OPEN-18)
bound_tm = y_tm(win_mask);
% OPEN-11, decided by N4V-D3 (DECISIONS_N4V.md): the legend text of a bound
% curve is the N-32 (b) observation label read from O1; the compressed curve
% also carries the note of N-32 (b) read from O1. The family, the geometry and
% the N-32 (a) label stay in the label columns of O1 and go to the caption (C2b).
lbl_fa = lab_fa.obs;
lbl_tm = [lab_tm.obs '; ' lab_tm.nb];

% =========================================================================
%  Y2: style, figure and axes, lines copied from the blob
% =========================================================================
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L80-84 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
c_B1  = [0.0000 0.4470 0.7410];   % blue
c_B2  = [0.8500 0.3250 0.0980];   % orange
c_B4  = [0.6350 0.0780 0.1840];   % crimson (proposed)
c_B5  = [0.4660 0.6740 0.1880];   % green
c_CRB = [0.0000 0.0000 0.0000];   % black
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L86-88 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
lw_std = 1.5;  lw_B4 = 2.0;  lw_CRB = 1.0;
ms_std = 8;    ms_B4 = 9;
mi = 1:3:n_pts;   % markers every 3rd point
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L93-94 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
fig_w = 8.8;
fig_h = 8.5;
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L96-101 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
fig = figure('Units', 'centimeters', ...
             'Position', [2 2 fig_w fig_h], ...
             'PaperUnits', 'centimeters', ...
             'PaperSize',  [fig_w fig_h], ...
             'PaperPosition', [0 0 fig_w fig_h], ...
             'Color', 'w');
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L103-105 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ax = axes('Parent', fig, ...
          'Units', 'normalized', ...
          'Position', [0.13 0.37 0.84 0.59]);

% the four estimator curves, as the published script draws them
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L110-130 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
hold(ax, 'on');

h_B1 = semilogy(ax, snr_B1, rmse_B1, '-o', ...
    'Color', c_B1, 'LineWidth', lw_std, ...
    'MarkerSize', ms_std, 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', c_B1, 'MarkerIndices', mi);

h_B2 = semilogy(ax, snr_B2, rmse_B2, '-s', ...
    'Color', c_B2, 'LineWidth', lw_std, ...
    'MarkerSize', ms_std, 'MarkerFaceColor', c_B2, ...
    'MarkerEdgeColor', c_B2, 'MarkerIndices', mi);

h_B4 = semilogy(ax, snr_B4, rmse_B4, '-d', ...
    'Color', c_B4, 'LineWidth', lw_B4, ...
    'MarkerSize', ms_B4, 'MarkerFaceColor', c_B4, ...
    'MarkerEdgeColor', c_B4, 'MarkerIndices', mi);

h_B5 = semilogy(ax, snr_B5, rmse_B5, '-^', ...
    'Color', c_B5, 'LineWidth', lw_std, ...
    'MarkerSize', ms_std, 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', c_B5, 'MarkerIndices', mi);

% the two bound curves: the published bound-curve call (L132-133), once per
% family; the compressed curve takes another line style (OPEN-11)
crb = bound_fa;
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L132-133 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
h_FA = semilogy(ax, snr_B4, crb, '--', ...
    'Color', c_CRB, 'LineWidth', lw_CRB);
crb = bound_tm;
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L132-133 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
h_TM = semilogy(ax, snr_B4, crb, '--', ...
    'Color', c_CRB, 'LineWidth', lw_CRB);
set(h_TM, 'LineStyle', '-.');
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L135 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
hold(ax, 'off');

% Y2; N4G-D10: axes limits, scales and ticks of the published window, labels
% (OPEN-27, decided by N4V-D2 (DECISIONS_N4V.md): the copied limits, then the y-window rule below them)
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L140-143 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ax.XLim    = [-5 17.5];
ax.XTick   = -5:2.5:17.5;
ax.YScale  = 'log';
ax.YLim    = [0.01 3];
% N4V-D2 (OPEN-27): the pre-registered y-window rule. The published lower limit
% is kept unless a plotted value lies below it; then the lower limit is the
% decade at or below the smallest plotted value. The upper limit is kept.
y_all = [rmse_B1(:); rmse_B2(:); rmse_B4(:); rmse_B5(:); bound_fa(:); bound_tm(:)];
y_all = y_all(isfinite(y_all) & y_all > 0);
if ~isempty(y_all) && min(y_all) < ax.YLim(1)
    ax.YLim = [10^floor(log10(min(y_all))), ax.YLim(2)];
end
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L145-153 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ax.FontSize = 8;
ax.FontName = 'Times New Roman';
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'off';
ax.YMinorTick = 'on';
ax.Box        = 'on';
grid(ax, 'on');
ax.GridAlpha     = 0.25;
ax.GridLineStyle = ':';
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L155-156 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
xlabel(ax, 'SNR [dB]',     'Interpreter', 'latex', 'FontSize', 9);
ylabel(ax, 'RMSE$_r$ [m]', 'Interpreter', 'latex', 'FontSize', 9);

% =========================================================================
%  Y3; N4G-D13: the legend (names only for the estimators), then its layout
% =========================================================================
% the option lines of the legend call are copied, except NumColumns, which N4V-D3 sets to 1 (a comment line cannot stand inside the continued statement)
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L168-172 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
leg = legend(ax, ...
    [h_B1, h_B2, h_B4, h_B5, h_FA, h_TM], ...
    {M_B1, M_B2, M_B4, M_B5, lbl_fa, lbl_tm}, ...
    'Interpreter', 'latex', ...
    'FontSize',    6.5, ...
    'NumColumns',  1, ...
    'Box',         'on', ...
    'Units',       'normalized');

% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L174-181 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
drawnow;

ax_pos     = ax.Position;
leg_h      = leg.Position(4);
gap        = 0.10;
leg_bottom = ax_pos(2) - gap - leg_h;
leg_bottom = max(leg_bottom, 0.01);
leg.Position = [ax_pos(1), leg_bottom, ax_pos(3), leg_h];

% =========================================================================
%  Y1, Y8: the PDF into results_dir under the name of N4G-D20, and the lines
%  printed (paths and the completion line; no value of any statistic)
% =========================================================================
out_pdf = fullfile(results_dir, ['new4_' PDF_BASE '_' stamp '.pdf']);
% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L187 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
print(fig, '-dpdf', '-painters', out_pdf);
fprintf('new4_plot_m1: wrote %s\n', out_pdf);
fprintf('NEW4 PLOT M1 COMPLETE: files read: %s ; %s ; file written: %s\n', f_csv, f_tab, out_pdf);
end

function [hdr, body] = loc_read_table(path)
% N4-I7: a CSV (written by new4_bounds.m, by the Monte Carlo driver or by
% wb_crb_d2_e3.m) read as text: one header line, then the rows; the fields of a
% row are split at commas outside double quotes. Reads no estimator output of a
% run: only the files named by the caller.
txt = fileread(path);
L = regexp(txt, '\r?\n', 'split');
L = L(~cellfun(@isempty, L));
hdr = loc_split_csv_line(L{1});
body = cell(numel(L) - 1, numel(hdr));
for i = 2:numel(L)
    f = loc_split_csv_line(L{i});
    if numel(f) ~= numel(hdr)
        error('new4_plot:csv', '%s: row %d has %d fields, the header has %d.', path, i, numel(f), numel(hdr));
    end
    body(i - 1, :) = f;
end
end

function f = loc_split_csv_line(s)
% the fields of one CSV line; a double-quoted field may hold commas and doubled quotes
f = {};
cur = '';
inq = false;
i = 1;
n = numel(s);
while i <= n
    c = s(i);
    if inq
        if c == '"'
            if i < n && s(i + 1) == '"'
                cur(end + 1) = '"';
                i = i + 1;
            else
                inq = false;
            end
        else
            cur(end + 1) = c;
        end
    else
        if c == '"'
            inq = true;
        elseif c == ','
            f{end + 1} = cur;
            cur = '';
        else
            cur(end + 1) = c;
        end
    end
    i = i + 1;
end
f{end + 1} = cur;
end

function j = loc_colname(hdr, name)
% the index of the one column whose name is name
k = find(strcmp(hdr, name));
if numel(k) ~= 1
    error('new4_plot:col', 'column %s: %d found, one required.', name, numel(k));
end
j = k;
end

function j = loc_col(hdr, prefix)
% the index of the one column whose name begins with prefix
k = find(strncmp(hdr, prefix, numel(prefix)));
if numel(k) ~= 1
    error('new4_plot:col', 'column with prefix %s: %d found, one required.', prefix, numel(k));
end
j = k;
end

function [x, y] = loc_series(hdr, body, method, xname, yname)
% the rows of one method of a Monte Carlo CSV, sorted by x
im = loc_colname(hdr, 'method');
sel = strcmp(body(:, im), method);
if ~any(sel)
    error('new4_plot:method', 'no row of method %s.', method);
end
x = str2double(body(sel, loc_colname(hdr, xname)));
y = str2double(body(sel, loc_colname(hdr, yname)));
[x, o] = sort(x);
y = y(o);
end

function [x, y, lab] = loc_family_rows(hdr, body, family, ycol)
% the rows of one family of a bound table of new4_bounds.m, sorted by snr_dB,
% with the label columns of its first row (the labels are read, never typed)
sel = find(strcmp(body(:, loc_colname(hdr, 'family')), family));
if isempty(sel)
    error('new4_plot:family', 'no row of family %s.', family);
end
x = str2double(body(sel, loc_colname(hdr, 'snr_dB')));
y = str2double(body(sel, ycol));
[x, o] = sort(x);
y = y(o);
lab.obs = body{sel(1), loc_colname(hdr, 'label_observation_N32b')};
lab.agg = body{sel(1), loc_colname(hdr, 'label_aggregation_N32a')};
lab.nb  = body{sel(1), loc_colname(hdr, 'label_not_a_bound_N32b')};
lab.geo = body{sel(1), loc_colname(hdr, 'geometry')};
lab.fam = strrep(family, '_', '-');
end
