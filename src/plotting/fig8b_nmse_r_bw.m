function fig8b_nmse_r_bw()
%FIG8B_NMSE_R_BW  Publication-quality plot of range NMSE vs bandwidth (Fig 8b).
%
%  Produces fig8b_nmse_r_bw.pdf for Paper C (Phase 3.5 / OJ-COMS submission).
%
%  Source CSVs (BOTH required):
%    PR-6 original : mc_bandwidth_sweep_20260520_180317.csv  (25 rows x 24 cols)
%                    BW = {100, 200, 400, 600, 800} MHz
%    PR-9 append   : mc_bandwidth_sweep_20260526_144912.csv  (15 rows x 24 cols)
%                    BW = {300, 500, 700} MHz
%  The script merges both on (B_hz, method) and sorts by B_hz -> 40 rows.
%
%  Output PDF : fig8b_nmse_r_bw.pdf  (single-column, 8.8 cm x 8.5 cm)
%
%  Five curves
%    B1  WB-BPD     (full-array, open blue circle,       1.5 pt)
%    B2  WB-P-SOMP  (compressed, filled orange square,   1.5 pt)
%    B4  WB-CL-KL   (compressed, filled crimson diamond, 2.0 pt) [PROPOSED]
%    B5  WB-DL-OMP  (full-array, open green triangle,    1.5 pt)
%    CRB compressed  (black dashed, no markers,           1.0 pt)
%
%  Markers at every BW point (8 points -- dense enough to show all).
%  x-axis: B in MHz, linear scale (log would compress the rollover region).
%  y-axis: NMSE_r in dB.
%  No title string.
%
%  ASCII compliance: no Unicode, no curly quotes, no Greek in strings.
%  Author  : R. V. Senyuva (Maltepe University)
%  Date    : 2026-05-26
%  Version : v1.0  (Phase 3.5 Subtask 3.5.2)

% =========================================================================
%  1.  LOCATE BOTH CSVs
% =========================================================================
script_dir = fileparts(mfilename('fullpath'));

csv_orig_name   = 'mc_bandwidth_sweep_20260520_180317.csv';
csv_append_name = 'mc_bandwidth_sweep_20260526_144912.csv';

% Original CSV: repo-root results\ (two levels up from src\plotting\)
csv_orig = fullfile(script_dir, '..', '..', 'results', ...
                    'paperC_phase3_fig8b_bw_sweep_v1', csv_orig_name);

% Append CSV: src\results\ (one level up from src\plotting\)
csv_app  = fullfile(script_dir, '..', 'results', ...
                    'fig8b_bw_append_v2', csv_append_name);

% Fallback: same directory as script (manual copy)
if ~isfile(csv_orig), csv_orig = fullfile(script_dir, csv_orig_name); end
if ~isfile(csv_app),  csv_app  = fullfile(script_dir, csv_append_name); end

assert(isfile(csv_orig), ['fig8b: original CSV not found: ' csv_orig]);
assert(isfile(csv_app),  ['fig8b: append CSV not found: '   csv_app]);
fprintf('fig8b: original CSV : %s\n', csv_orig);
fprintf('fig8b: append  CSV  : %s\n', csv_app);

% =========================================================================
%  2.  READ AND MERGE
% =========================================================================
T_orig = readtable(csv_orig, 'TextType', 'string');
T_app  = readtable(csv_app,  'TextType', 'string');
T = [T_orig; T_app];   % vertical concat; sort done per-method below

    function [bv, nv, cv] = extract(T, mname)
        rows = T(strcmp(T.method, mname), :);
        rows = sortrows(rows, 'B_hz');
        bv = rows.B_hz / 1e6;   % Hz -> MHz
        nv = rows.NMSE_r_dB;
        cv = rows.crb_r_m;
    end

[bw_B1, nmse_B1, ~    ] = extract(T, 'WB-BPD');
[bw_B2, nmse_B2, ~    ] = extract(T, 'WB-P-SOMP');
[bw_B4, nmse_B4, crb  ] = extract(T, 'WB-CL-KL');
[bw_B5, nmse_B5, ~    ] = extract(T, 'WB-DL-OMP');

% CRB in dB relative to scene-average r^2 is not directly available;
% plot sqrt-CRB as NMSE: NMSE_CRB = 20*log10(CRB_r / mean(r_true))
% However, the most consistent comparison is to plot the CRB_r column
% converted to NMSE using the same normalisation as the MC results.
% The MC uses NMSE_r_dB = 20*log10(RMSE_r / r_true_mean).
% r_true is drawn uniformly from [r_lo, r_hi] = [1.0631, 21.2625] m
% -> E[r] = (1.0631 + 21.2625)/2 = 11.1628 m.
r_mean = 11.1628;   % m  (scene-average range for NMSE normalisation)
crb_nmse = 20*log10(crb / r_mean);   % compressed CRB as NMSE_r (dB)

n_pts = numel(bw_B4);
fprintf('fig8b: %d BW points loaded after merge.\n', n_pts);

% =========================================================================
%  3.  STYLE  (Phase 3.5 Plan Section 2)
% =========================================================================
c_B1  = [0.0000 0.4470 0.7410];
c_B2  = [0.8500 0.3250 0.0980];
c_B4  = [0.6350 0.0780 0.1840];
c_B5  = [0.4660 0.6740 0.1880];
c_CRB = [0.0000 0.0000 0.0000];

lw_std = 1.5;  lw_B4 = 2.0;  lw_CRB = 1.0;
ms_std = 8;    ms_B4 = 9;

% All 8 BW points get markers (8 pts -- not overcrowded)
mi = 1:n_pts;

% =========================================================================
%  4.  FIGURE SETUP  (identical dimensions to Figs 8a / 9a)
% =========================================================================
fig_w = 8.8;
fig_h = 8.5;

fig = figure('Units', 'centimeters', ...
             'Position', [2 2 fig_w fig_h], ...
             'PaperUnits', 'centimeters', ...
             'PaperSize',  [fig_w fig_h], ...
             'PaperPosition', [0 0 fig_w fig_h], ...
             'Color', 'w');

ax = axes('Parent', fig, ...
          'Units', 'normalized', ...
          'Position', [0.13 0.37 0.84 0.59]);

% =========================================================================
%  5.  PLOT CURVES
% =========================================================================
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

h_CRB = plot(ax, bw_B4, crb_nmse, '--', ...
    'Color', c_CRB, 'LineWidth', lw_CRB);

hold(ax, 'off');

% =========================================================================
%  6.  AXES FORMATTING
% =========================================================================
ax.XLim    = [50 850];
ax.XTick   = [100 200 300 400 500 600 700 800];
ax.YScale  = 'linear';   % NMSE already in dB -- linear axis

% Set y-limits to cover all curves with some margin
all_nmse = [nmse_B1; nmse_B2; nmse_B4; nmse_B5; crb_nmse];
y_lo = floor(min(all_nmse(isfinite(all_nmse))) / 5) * 5 - 5;
y_hi = ceil( max(all_nmse(isfinite(all_nmse))) / 5) * 5 + 5;
ax.YLim = [y_lo y_hi];

ax.FontSize = 8;
ax.FontName = 'Times New Roman';
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'off';
ax.YMinorTick = 'on';
ax.Box        = 'on';
grid(ax, 'on');
ax.GridAlpha     = 0.25;
ax.GridLineStyle = ':';

xlabel(ax, '$B$ [MHz]',      'Interpreter', 'latex', 'FontSize', 9);
ylabel(ax, 'NMSE$_r$ [dB]',  'Interpreter', 'latex', 'FontSize', 9);

% =========================================================================
%  7.  LEGEND  -- 2-column, outside below axes, aligned with x-axis span
% =========================================================================
leg = legend(ax, ...
    [h_B1, h_B2, h_B4, h_B5, h_CRB], ...
    {'WB-BPD (full)', ...
     'WB-P-SOMP (comp.)', ...
     'WB-CL-KL (comp.) [prop.]', ...
     'WB-DL-OMP (full)', ...
     'CRB (comp., fixed $r$)'}, ...
    'Interpreter', 'latex', ...
    'FontSize',    6.5, ...
    'NumColumns',  2, ...
    'Box',         'on', ...
    'Units',       'normalized');

drawnow;

ax_pos   = ax.Position;
ax_left  = ax_pos(1);
ax_width = ax_pos(3);
ax_bot   = ax_pos(2);

leg_h      = leg.Position(4);
gap        = 0.10;
leg_bottom = ax_bot - gap - leg_h;
leg_bottom = max(leg_bottom, 0.01);

leg.Position = [ax_left, leg_bottom, ax_width, leg_h];

% =========================================================================
%  8.  EXPORT
% =========================================================================
out_pdf = 'fig8b_nmse_r_bw.pdf';
print(fig, '-dpdf', '-painters', out_pdf);
fprintf('fig8b: exported to %s\n', out_pdf);

end  % function fig8b_nmse_r_bw
