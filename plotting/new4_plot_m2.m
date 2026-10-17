function new4_plot_m2(results_dir, stamp)
%NEW4_PLOT_M2  M2 (fig:NMSEr_BW, fig8b_nmse_r_bw.pdf): the range NMSE of the four
%  estimators against bandwidth, re-drawn from the committed CSV, with no bound curve.
%
%  new4_plot_m2(results_dir, stamp)
%
%  Written by the NEW-4 code session of paperC-lab (sessions/2026-10-18_NEW4-CODE)
%  from NEW4_SPEC.md as decided in DECISIONS_N4G.md (N4G-D2, D10, D13, D15, D20).
%  Written as text; never parsed or executed by its author (N6).
%
%  Reads only the committed CSV of M2 (the estimator columns). Draws no bound
%  curve and no mean range (N4G-D2, N4G-D10: the published conversion of a bound
%  to an NMSE level lapses). Computes nothing from the code repository. Writes
%  one PDF into results_dir, named new4_<published PDF base name>_<stamp>.pdf
%  (Y1; N4G-D20).
%
%  Labels in comments: Yn = requirement of the starter of the code session;
%  N4G-Dn = author decisions; OPEN-n = a choice no record fixes, listed in
%  OPEN_NEW4.md.

% =========================================================================
%  Y1: arguments and named constants; every path is relative to the root of
%  the code repository (D-20; N4G-D17)
% =========================================================================
if nargin ~= 2 || ~ischar(results_dir) || isempty(results_dir) || ~ischar(stamp) || isempty(stamp)
    error('new4_plot_m2:usage', 'new4_plot_m2(results_dir, stamp): two non-empty character vectors are required.');
end
here_dir = fileparts(mfilename('fullpath'));    % Y1: the file stands in plotting/ of the code repository
root_dir = fileparts(here_dir);
CSV_M2   = fullfile('results', 'paperC_phase35_fig8b_bw_r020_N600', 'mc_bandwidth_sweep_20260529_210446.csv');   % CR Sec. 3.2: the committed CSV of M2
PDF_BASE = 'fig8b_nmse_r_bw';   % COPIED from L:plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L204
M_B1 = 'WB-BPD';          % COPIED from L:plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L64 (N4G-D13: names only)
M_B2 = 'WB-P-SOMP';    % COPIED from L:plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L65
M_B4 = 'WB-CL-KL';      % COPIED from L:plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L66
M_B5 = 'WB-DL-OMP';    % COPIED from L:plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L67
f_csv = fullfile(root_dir, CSV_M2);

% =========================================================================
%  Y5; N4G-D2, N4G-D13: the estimator curves, from the committed CSV
% =========================================================================
fprintf('new4_plot_m2: reading %s\n', f_csv);
[hdr_e, body_e] = loc_read_table(f_csv);
[bw_B1, nmse_B1] = loc_series(hdr_e, body_e, M_B1, 'B_hz', 'NMSE_r_dB');
[bw_B2, nmse_B2] = loc_series(hdr_e, body_e, M_B2, 'B_hz', 'NMSE_r_dB');
[bw_B4, nmse_B4] = loc_series(hdr_e, body_e, M_B4, 'B_hz', 'NMSE_r_dB');
[bw_B5, nmse_B5] = loc_series(hdr_e, body_e, M_B5, 'B_hz', 'NMSE_r_dB');
bw_B1 = bw_B1 / 1e6;                            % Hz -> MHz, as the published script
bw_B2 = bw_B2 / 1e6;
bw_B4 = bw_B4 / 1e6;
bw_B5 = bw_B5 / 1e6;
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L86 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
n_pts = numel(bw_B4);

% =========================================================================
%  Y2: style, figure and axes, lines copied from the blob
% =========================================================================
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L92-96 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
c_B1  = [0.0000 0.4470 0.7410];
c_B2  = [0.8500 0.3250 0.0980];
c_B4  = [0.6350 0.0780 0.1840];
c_B5  = [0.4660 0.6740 0.1880];
c_CRB = [0.0000 0.0000 0.0000];
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L98-100 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
lw_std = 1.5;  lw_B4 = 2.0;  lw_CRB = 1.0;
ms_std = 8;    ms_B4 = 9;
mi = 1:n_pts;   % all 8 BW points get markers
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L105-106 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
fig_w = 8.8;
fig_h = 8.5;
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L108-113 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
fig = figure('Units', 'centimeters', ...
             'Position', [2 2 fig_w fig_h], ...
             'PaperUnits', 'centimeters', ...
             'PaperSize',  [fig_w fig_h], ...
             'PaperPosition', [0 0 fig_w fig_h], ...
             'Color', 'w');
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L115-117 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ax = axes('Parent', fig, ...
          'Units', 'normalized', ...
          'Position', [0.13 0.37 0.84 0.59]);

% the four estimator curves, as the published script draws them (no bound curve)
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L122-142 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
hold(ax, 'on');

h_B1 = plot(ax, bw_B1, nmse_B1, '-o', ...
    'Color', c_B1, 'LineWidth', lw_std, ...
    'MarkerSize', ms_std, 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', c_B1, 'MarkerIndices', mi);

h_B2 = plot(ax, bw_B2, nmse_B2, '-s', ...
    'Color', c_B2, 'LineWidth', lw_std, ...
    'MarkerSize', ms_std, 'MarkerFaceColor', c_B2, ...
    'MarkerEdgeColor', c_B2, 'MarkerIndices', mi);

h_B4 = plot(ax, bw_B4, nmse_B4, '-d', ...
    'Color', c_B4, 'LineWidth', lw_B4, ...
    'MarkerSize', ms_B4, 'MarkerFaceColor', c_B4, ...
    'MarkerEdgeColor', c_B4, 'MarkerIndices', mi);

h_B5 = plot(ax, bw_B5, nmse_B5, '-^', ...
    'Color', c_B5, 'LineWidth', lw_std, ...
    'MarkerSize', ms_std, 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', c_B5, 'MarkerIndices', mi);
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L147 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
hold(ax, 'off');

% Y2; N4G-D10: axes limits, scales and ticks of the published window; the y limits
% are set by the published rule from the plotted curves, here the four estimator
% curves only, the bound curve being absent (OPEN-26); labels
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L152-154 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ax.XLim    = [50 850];
ax.XTick   = 100:100:800;
ax.YScale  = 'linear';
all_nmse = [nmse_B1; nmse_B2; nmse_B4; nmse_B5];
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L158-161 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
all_nmse = all_nmse(isfinite(all_nmse));
y_lo = floor(min(all_nmse) / 5) * 5 - 5;
y_hi = ceil( max(all_nmse) / 5) * 5 + 5;
ax.YLim = [y_lo y_hi];
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L163-171 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ax.FontSize = 8;
ax.FontName = 'Times New Roman';
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'off';
ax.YMinorTick = 'on';
ax.Box        = 'on';
grid(ax, 'on');
ax.GridAlpha     = 0.25;
ax.GridLineStyle = ':';
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L173-174 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
xlabel(ax, '$B$ [MHz]',      'Interpreter', 'latex', 'FontSize', 9);
ylabel(ax, 'NMSE$_r$ [dB]',  'Interpreter', 'latex', 'FontSize', 9);

% =========================================================================
%  Y5; N4G-D13: the legend (names only), then its layout
% =========================================================================
% the option lines of the legend call are copied (a comment line cannot stand inside the continued statement)
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L186-190 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
leg = legend(ax, ...
    [h_B1, h_B2, h_B4, h_B5], ...
    {M_B1, M_B2, M_B4, M_B5}, ...
    'Interpreter', 'latex', ...
    'FontSize',    6.5, ...
    'NumColumns',  2, ...
    'Box',         'on', ...
    'Units',       'normalized');

% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L192-199 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
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
% COPIED from plotting/fig8b_nmse_r_bw.m@6d85800c4bd3 L205 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
print(fig, '-dpdf', '-painters', out_pdf);
fprintf('new4_plot_m2: wrote %s\n', out_pdf);
fprintf('NEW4 PLOT M2 COMPLETE: file read: %s ; file written: %s\n', f_csv, out_pdf);
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
