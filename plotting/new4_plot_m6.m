function new4_plot_m6(results_dir, stamp)
%NEW4_PLOT_M6  The new M6 (fig:compCRB_SNR, fig6_crb_d2_snr.pdf): the range bound
%  of each of the two paths against SNR, one series per N_RF, from the E3 per-point file.
%
%  new4_plot_m6(results_dir, stamp)
%
%  Written by the NEW-4 code session of paperC-lab (sessions/2026-10-18_NEW4-CODE)
%  from NEW4_SPEC.md Sec. d as decided in DECISIONS_N4G.md (N4G-D13, D14, D15, D20).
%  Written as text; never parsed or executed by its author (N6).
%
%  Reads only the copy of the E3 per-point file F2 in results/paperC_E3_C4 of the
%  code repository (N4G-D14; the run sheet's gate N4-G7 ties its SHA-256 to the
%  committed capture before this function runs). Draws the range aggregate of
%  each path (variance domain, full-Jacobian form, equilibrated exact inverse,
%  as its column names say); draws no angle aggregate, no C_J column and no
%  marking of any point. Computes nothing from the code repository. Writes one
%  PDF into results_dir, named new4_<published PDF base name>_<stamp>.pdf (Y1;
%  N4G-D20).
%
%  Labels in comments: Yn = requirement of the starter of the code session;
%  N4G-Dn = author decisions; OPEN-n = a choice no record fixes, listed in
%  OPEN_NEW4.md.

% =========================================================================
%  Y1: arguments and named constants; every path is relative to the root of
%  the code repository (D-20; N4G-D17)
% =========================================================================
if nargin ~= 2 || ~ischar(results_dir) || isempty(results_dir) || ~ischar(stamp) || isempty(stamp)
    error('new4_plot_m6:usage', 'new4_plot_m6(results_dir, stamp): two non-empty character vectors are required.');
end
here_dir = fileparts(mfilename('fullpath'));    % Y1: the file stands in plotting/ of the code repository
root_dir = fileparts(here_dir);
CSV_F2   = fullfile('results', 'paperC_E3_C4', 'wb_crb_d2_e3_perpoint_20261016_101741.csv');   % N4G-D14, SPEC [S23]: the copy of F2 in the code repository
PDF_BASE = 'fig6_crb_d2_snr';   % COPIED from L:plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L192
LEG_WRAP = 70;                                  % OPEN-13: the label of N-239 is broken into lines of at most this many characters
fig_h_cm = 12;                                  % OPEN-15: the height of the figure in cm (the new M6 has no window of record)
f_f2 = fullfile(root_dir, CSV_F2);

% =========================================================================
%  Y6; N4G-D13, N4G-D14: the columns of F2 (names found by prefix), the series
% =========================================================================
fprintf('new4_plot_m6: reading %s\n', f_f2);
[hdr, body] = loc_read_table(f_f2);
i_nrf = loc_colname(hdr, 'N_RF');
i_snr = loc_colname(hdr, 'snr_dB');
i_M   = loc_colname(hdr, 'M');
i_r1  = loc_col(hdr, 'range_aggregate_path1.');
i_r2  = loc_col(hdr, 'range_aggregate_path2.');
i_lab = loc_col(hdr, 'range_aggregate_label.');
nrf_col = str2double(body(:, i_nrf));
snr_col = str2double(body(:, i_snr));
M_col   = str2double(body(:, i_M));
r1_col  = str2double(body(:, i_r1));
r2_col  = str2double(body(:, i_r2));
lab_txt = body{1, i_lab};                       % the label of N-239, as the run wrote it (read, never typed)
if any(~strcmp(body(:, i_lab), lab_txt))
    error('new4_plot_m6:label', 'the label column of F2 differs between rows.');
end
nrf_vals = unique(nrf_col).';                   % one series per N_RF value of F2, ascending
n_series = numel(nrf_vals);
n_snr    = numel(unique(snr_col));

% =========================================================================
%  Y2: style, figure and axes, lines copied from the blob (the figure height
%  replaced by fig_h_cm; the two panels, one per path, are new: OPEN-15)
% =========================================================================
% COPIED from plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L92-94 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
clr8  = [0.6350 0.0780 0.1840];   % crimson
clr16 = [0.0000 0.3176 0.6275];   % navy
clr32 = [0.1176 0.4706 0.1176];   % forest green
% COPIED from plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L95 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
clr64 = [0.3000 0.3000 0.3000];   % dark gray
% COPIED from plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L97 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
lw_comp = 1.5;    % compressed curves
% COPIED from plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L99 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
msz     = 7;
% COPIED from plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L100 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
mk_idx  = 1 : 3 : n_snr;    % markers at indices 1,4,7,10,13
palette = {clr8, clr16, clr32, clr64};
mk1 = 'o';   % COPIED from L:plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L121 (the marker of the first series)
mk2 = 's';   % COPIED from L:plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L127
mk3 = '^';   % COPIED from L:plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L133
mk4 = 'd';   % COPIED from L:plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L140
markers = {mk1, mk2, mk3, mk4};
if n_series > numel(palette)
    error('new4_plot_m6:series', 'F2 holds %d N_RF values, the style of the published script has %d.', n_series, numel(palette));
end

% COPIED from plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L105-107 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
fig = figure('Units', 'centimeters', 'Position', [2 2 8.8 fig_h_cm]);
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [8.8 fig_h_cm], ...
         'PaperPosition', [0 0 8.8 fig_h_cm]);
ax_top = axes('Parent', fig, 'Position', [0.15 0.70 0.80 0.25]);
ax     = axes('Parent', fig, 'Position', [0.15 0.40 0.80 0.25]);
% COPIED from plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L110-113 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
set(ax_top, 'YScale', 'log', 'Box', 'on', 'FontSize', 7, ...
        'FontName', 'Times New Roman');
hold(ax_top, 'on');
grid(ax_top, 'on');
% COPIED from plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L110-113 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
set(ax, 'YScale', 'log', 'Box', 'on', 'FontSize', 7, ...
        'FontName', 'Times New Roman');
hold(ax, 'on');
grid(ax, 'on');

% =========================================================================
%  Y6; N4G-D13: one curve per N_RF in each panel; the series with N_RF = M is
%  named N_RF = M and by no other name
% =========================================================================
h_top = gobjects(1, n_series);
h_bot = gobjects(1, n_series);
names = cell(1, n_series);
for k = 1:n_series
    sel = find(nrf_col == nrf_vals(k));
    [xs, o] = sort(snr_col(sel));
    y1 = r1_col(sel); y1 = y1(o);
    y2 = r2_col(sel); y2 = y2(o);
    h_top(k) = plot(ax_top, xs, y1, '-', 'Color', palette{k}, 'LineWidth', lw_comp, ...
        'Marker', markers{k}, 'MarkerSize', msz, 'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', palette{k}, 'MarkerIndices', mk_idx);
    h_bot(k) = plot(ax, xs, y2, '-', 'Color', palette{k}, 'LineWidth', lw_comp, ...
        'Marker', markers{k}, 'MarkerSize', msz, 'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', palette{k}, 'MarkerIndices', mk_idx);
    if nrf_vals(k) == M_col(sel(1))
        names{k} = sprintf('$N_{\\rm RF} = %d = M$', nrf_vals(k));
    else
        names{k} = sprintf('$N_{\\rm RF} = %d$', nrf_vals(k));
    end
end

% =========================================================================
%  Y6: axes (logarithmic range axis, set above), limits from the data, labels,
%  panel titles (one per path), then the legend with the label of N-239
% =========================================================================
xlim(ax_top, [min(snr_col) max(snr_col)]);
xlim(ax, [min(snr_col) max(snr_col)]);
% COPIED from plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L172 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
xlabel(ax, 'SNR [dB]', 'FontSize', 7.5, 'FontName', 'Times New Roman');
% COPIED from plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L174 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ylabel(ax_top, 'range bound [m]', ...
    'FontSize', 7.5, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
% COPIED from plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L174 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
ylabel(ax, 'range bound [m]', ...
    'FontSize', 7.5, 'FontName', 'Times New Roman', 'Interpreter', 'latex');
title(ax_top, 'path 1', 'FontSize', 7, 'FontName', 'Times New Roman');
title(ax, 'path 2', 'FontSize', 7, 'FontName', 'Times New Roman');

% COPIED from plotting/fig8a_rmse_r_snr.m@0ebc02da00bd L168-172 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
leg = legend(ax, h_bot, names, ...
    'Interpreter', 'latex', ...
    'FontSize',    6.5, ...
    'NumColumns',  2, ...
    'Box',         'on', ...
    'Units',       'normalized');
leg.Title.String = loc_wrap(lab_txt, LEG_WRAP);   % OPEN-13: the label of N-239 is the title of the legend of the series

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
pdf_path = fullfile(results_dir, ['new4_' PDF_BASE '_' stamp '.pdf']);
% COPIED from plotting/fig6_crb_d2_snr.m@13fae0aa9ec6 L193 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
fprintf('new4_plot_m6: wrote %s\n', pdf_path);
fprintf('NEW4 PLOT M6 COMPLETE: file read: %s ; file written: %s\n', f_f2, pdf_path);
end

function c = loc_wrap(s, w)
% the text s broken at spaces into lines of at most w characters
words = strsplit(s, ' ');
c = {};
cur = '';
for i = 1:numel(words)
    cand = strtrim([cur ' ' words{i}]);
    if numel(cand) > w && ~isempty(cur)
        c{end + 1} = cur;
        cur = words{i};
    else
        cur = cand;
    end
end
c{end + 1} = cur;
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
