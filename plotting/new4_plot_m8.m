function new4_plot_m8(results_dir, stamp)
%NEW4_PLOT_M8  M8 (fig:3GPP_UMi_SNR, fig11_riviello_snr_axis.pdf): the range RMSE of
%  WB-CL-KL and the full-array range bound at one stated geometry against SNR,
%  and the histogram of the 3GPP UMi SNR distribution.
%
%  new4_plot_m8(results_dir, stamp)
%
%  Written by the NEW-4 code session of paperC-lab (sessions/2026-10-18_NEW4-CODE)
%  from NEW4_SPEC.md as decided in DECISIONS_N4G.md (N4G-D6, D8 (ii), D10, D13, D15, D20).
%  Written as text; never parsed or executed by its author (N6).
%
%  Reads only: the committed CSV of M8 (the estimator column), the committed
%  distribution file of the lower panel (the data of record, N4G-D8 (ii)) and the
%  bound table new4_bound_m8_<stamp>.csv of new4_bounds.m in results_dir.
%  Computes no bound, calls no function of the code repository. It forms, draws,
%  prints and writes no quotient or difference of an estimator statistic and a
%  bound: the lines of the published script that form and write one (its
%  L40-41 and L89-93) are not copied (inherited-finding block of Ledger N-195).
%  The axis label and the legend entry name the RMSE, the quantity plotted
%  (N4G-D6, N4G-D13). Writes one PDF into results_dir, named
%  new4_<published PDF base name>_<stamp>.pdf (Y1; N4G-D20).
%
%  Labels in comments: Yn = requirement of the starter of the code session;
%  N4G-Dn = author decisions; OPEN-n = a choice no record fixes, listed in
%  OPEN_NEW4.md.

% =========================================================================
%  Y1: arguments and named constants; every path is relative to the root of
%  the code repository (D-20; N4G-D17)
% =========================================================================
if nargin ~= 2 || ~ischar(results_dir) || isempty(results_dir) || ~ischar(stamp) || isempty(stamp)
    error('new4_plot_m8:usage', 'new4_plot_m8(results_dir, stamp): two non-empty character vectors are required.');
end
here_dir = fileparts(mfilename('fullpath'));    % Y1: the file stands in plotting/ of the code repository
root_dir = fileparts(here_dir);
CSV_SWEEP  = fullfile('results', 'riviello_snr', 'mc_snr_sweep_20260610_203626.csv');           % CR Sec. 3.6: the committed CSV of M8
CSV_DIST   = fullfile('results', 'riviello_snr', 'riviello_snr_distribution_2000.csv');         % N4G-D8 (ii): the data of record of the lower panel (header-less)
TABLE_FILE = ['new4_bound_m8_' stamp '.csv'];                                                    % the bound table of M8 (O1)
PDF_BASE   = 'fig11_riviello_snr_axis';   % COPIED from L:plotting/fig11_riviello_snr_axis.m@759a69e4f8df L25
M_B4 = 'WB-CL-KL';          % COPIED from L:plotting/fig11_riviello_snr_axis.m@759a69e4f8df L31
f_csv  = fullfile(root_dir, CSV_SWEEP);
f_tab  = fullfile(results_dir, TABLE_FILE);
DISTRIBUTION_CSV = fullfile(root_dir, CSV_DIST);
OUT_PDF = fullfile(results_dir, ['new4_' PDF_BASE '_' stamp '.pdf']);

% =========================================================================
%  Y7; N4G-D6, N4G-D13: the estimator curve (committed CSV) and the bound (O1)
% =========================================================================
fprintf('new4_plot_m8: reading %s\n', f_csv);
[hdr_e, body_e] = loc_read_table(f_csv);
[snr_pts, rmse_r] = loc_series(hdr_e, body_e, M_B4, 'SNR_dB', 'RMSE_r_m');

fprintf('new4_plot_m8: reading %s\n', f_tab);
[hdr_b, body_b] = loc_read_table(f_tab);
snr_b = str2double(body_b(:, loc_colname(hdr_b, 'snr_dB')));
bound_r = str2double(body_b(:, loc_col(hdr_b, 'range_bound.')));
[snr_b, o] = sort(snr_b);
bound_r = bound_r(o);
if numel(snr_b) ~= numel(snr_pts) || any(abs(snr_b - snr_pts) > 1e-9)
    error('new4_plot_m8:grid', 'the SNR points of the bound table differ from those of the estimator rows.');
end
lab_obs = body_b{1, loc_colname(hdr_b, 'label_observation_N32b')};
lab_sg  = body_b{1, loc_colname(hdr_b, 'label_single_geometry_N32c')};
th_nom  = str2double(body_b{1, loc_colname(hdr_b, 'theta_nom_deg')});
r_nom   = str2double(body_b{1, loc_colname(hdr_b, 'r_nom_m')});
% OPEN-14: the legend entries; the bound entry is the label read from O1, then the
% word of N-32 (c) read from O1 with the geometry of the curve; the estimator
% entry names the RMSE
lbl_bd = [lab_obs ', ' lab_sg ' (' sprintf('%.6g deg, %.6g m', th_nom, r_nom) ')'];
lbl_b4 = [M_B4 ' $\mathrm{RMSE}_r$'];
crb_r = bound_r;                                % the variable name of the published bound-curve call below

fprintf('new4_plot_m8: reading %s\n', DISTRIBUTION_CSV);
% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L36-37 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
snr_dist   = readmatrix(DISTRIBUTION_CSV);
snr_dist   = snr_dist(:);

% Y7, Y2; N4G-D8 (ii): the histogram of the distribution of record, lines copied from the blob
% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L47-50 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
bin_edges = (-20 : 5 : 40);
bin_ctrs  = bin_edges(1:end-1) + 2.5;
counts    = histcounts(snr_dist, bin_edges);
density   = counts / (numel(snr_dist) * 5);   % fraction per dB

% =========================================================================
%  Y2; N4G-D10: layout, window and figure, lines copied from the blob
% =========================================================================
% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L55-56 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
x_lo = -20;   x_hi = 35;
y_lo = 4e-4;  y_hi = 2e0;
% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L58-59 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
fig_w_cm = 8.8;
fig_h_cm = 10.5;
% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L62-63 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ax_main_pos = [0.14  0.38  0.83  0.58];   % top: curves
ax_hist_pos = [0.14  0.11  0.83  0.20];   % bottom: histogram
% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L65-67 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
fig = figure('Units', 'centimeters', ...
    'Position', [2 2 fig_w_cm fig_h_cm], ...
    'Color', 'w');

% =========================================================================
%  Y7: the top panel, the published calls with the curve names changed; the
%  published median line, its label and the annotation of a quotient are not
%  drawn (OPEN-14)
% =========================================================================
% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L72 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ax1 = axes('Parent', fig, 'Position', ax_main_pos);

% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L74-75 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
h_bd = semilogy(ax1, snr_pts, crb_r, 'k-o', ...
    'LineWidth', 1.4, 'MarkerSize', 4, 'MarkerFaceColor', 'k', ...
    'DisplayName', lbl_bd);
% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L78 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
hold(ax1, 'on');

% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L80-82 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
h_b4 = semilogy(ax1, snr_pts, rmse_r, 'b-s', ...
    'LineWidth', 1.4, 'MarkerSize', 4, ...
    'MarkerFaceColor', [0.15 0.35 0.75], ...
    'DisplayName', lbl_b4);

% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L101-109 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ax1.XLim   = [x_lo, x_hi];
ax1.YLim   = [y_lo, y_hi];
ax1.YTick  = [1e-3 1e-2 1e-1 1e0];
ax1.XTick  = x_lo : 5 : x_hi;
ax1.XTickLabel = {};
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 8;
grid(ax1, 'on');
ax1.GridAlpha = 0.25;

% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L112 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ylabel(ax1, '$\mathrm{RMSE}_r$ and range bound (m)', ...
    'Interpreter', 'latex', 'FontSize', 8.5);

% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L116-119 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
legend(ax1, [h_bd, h_b4], ...
    'Location',    'northeast', ...
    'Interpreter', 'latex', ...
    'FontSize',    6.5, ...
    'Box',         'on');

% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L121 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
hold(ax1, 'off');

% =========================================================================
%  Y7: the bottom panel, the histogram of the distribution of record
% =========================================================================
% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L126 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ax2 = axes('Parent', fig, 'Position', ax_hist_pos);

% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L128-132 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
bar(ax2, bin_ctrs, density, 1.0, ...
    'FaceColor', [0.65 0.80 0.92], ...
    'EdgeColor', [0.40 0.60 0.78], ...
    'FaceAlpha', 0.85);
hold(ax2, 'on');

% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L139-143 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
text(ax2, x_hi - 0.5, max(density) * 1.20, ...
    '3GPP UMi SNR density', ...
    'Interpreter', 'latex', 'FontSize', 7, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'Color', [0.20 0.20 0.20]);

% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L145-156 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ax2.XLim   = [x_lo, x_hi];
ax2.YLim   = [0, max(density) * 1.40];
ax2.XTick  = x_lo : 5 : x_hi;
ax2.YTick  = [0 0.01 0.02 0.03];
ax2.TickLabelInterpreter = 'latex';
ax2.FontSize = 8;
grid(ax2, 'on');
ax2.GridAlpha = 0.25;

xlabel(ax2, 'SNR (dB)', 'Interpreter', 'latex', 'FontSize', 8.5);
ylabel(ax2, 'Density (dB$^{-1}$)', 'Interpreter', 'latex', 'FontSize', 8);
hold(ax2, 'off');

% =========================================================================
%  Y1, Y8: the PDF into results_dir under the name of N4G-D20, and the lines
%  printed (paths and the completion line; no value of any statistic)
% =========================================================================
% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L165-167 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
set(fig, 'PaperUnits', 'centimeters', ...
    'PaperSize',     [fig_w_cm fig_h_cm], ...
    'PaperPosition', [0 0 fig_w_cm fig_h_cm]);
% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L168 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
print(fig, OUT_PDF, '-dpdf', '-r300');
fprintf('new4_plot_m8: wrote %s\n', OUT_PDF);
% COPIED from plotting/fig11_riviello_snr_axis.m@759a69e4f8df L170 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
close(fig);
fprintf('NEW4 PLOT M8 COMPLETE: files read: %s ; %s ; %s ; file written: %s\n', f_csv, DISTRIBUTION_CSV, f_tab, OUT_PDF);
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
