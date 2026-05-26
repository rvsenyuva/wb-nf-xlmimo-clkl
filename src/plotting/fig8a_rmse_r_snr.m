function fig8a_rmse_r_snr()
%FIG8A_RMSE_R_SNR  Publication-quality plot of range RMSE vs SNR (Fig 8a).
%
%  Produces fig8a_rmse_r_snr.pdf for Paper C (Phase 3.5 / OJ-COMS submission).
%
%  Source CSV : mc_snr_sweep_20260526_090431.csv
%               (65 rows x 24 cols; 13 SNR pts at -5:2.5:25 dB; N_MC=200)
%  Output PDF : fig8a_rmse_r_snr.pdf  (single-column, 8.8 cm x 6.5 cm)
%
%  ASCII compliance: no Unicode, no curly quotes, no Greek in strings.
%  Author  : R. V. Senyuva (Maltepe University)
%  Date    : 2026-05-26
%  Version : v1.3  (legend scaled down, centred below axes, gap increased)

% =========================================================================
%  1.  LOCATE CSV
% =========================================================================
csv_name = 'mc_snr_sweep_20260526_090431.csv';
script_dir = fileparts(mfilename('fullpath'));
csv_path = fullfile(script_dir, '..', 'results', 'fig8a_snr_v2', csv_name);
if ~isfile(csv_path)
    csv_path = fullfile('results', 'fig8a_snr_v2', csv_name);
end
if ~isfile(csv_path)
    csv_path = fullfile(script_dir, csv_name);
end
assert(isfile(csv_path), ...
    ['fig8a: CSV not found. Last tried: ' csv_path]);
fprintf('fig8a: reading %s\n', csv_path);

% =========================================================================
%  2.  READ CSV
% =========================================================================
T = readtable(csv_path, 'TextType', 'string');

    function [sv, rv, cv] = extract(T, mname)
        rows = T(strcmp(T.method, mname), :);
        rows = sortrows(rows, 'SNR_dB');
        sv = rows.SNR_dB;
        rv = rows.RMSE_r_m;
        cv = rows.crb_r_m;
    end

[snr_B1, rmse_B1, ~  ] = extract(T, 'WB-BPD');
[snr_B2, rmse_B2, ~  ] = extract(T, 'WB-P-SOMP');
[snr_B4, rmse_B4, crb] = extract(T, 'WB-CL-KL');
[snr_B5, rmse_B5, ~  ] = extract(T, 'WB-DL-OMP');

n_pts = numel(snr_B4);
fprintf('fig8a: %d SNR points loaded.\n', n_pts);

% =========================================================================
%  3.  STYLE
% =========================================================================
c_B1  = [0.0000 0.4470 0.7410];
c_B2  = [0.8500 0.3250 0.0980];
c_B4  = [0.6350 0.0780 0.1840];
c_B5  = [0.4660 0.6740 0.1880];
c_CRB = [0.0000 0.0000 0.0000];

lw_std = 1.5;  lw_B4 = 2.0;  lw_CRB = 1.0;
ms_std = 8;    ms_B4 = 9;
mi = 1:3:n_pts;

% =========================================================================
%  4.  FIGURE SETUP
%      fig_h = 8.5 cm: axes get ~5.5 cm, gap ~0.6 cm, legend ~1.6 cm,
%      bottom margin ~0.4 cm.
%      Axes top edge at normalized 0.97, bottom at 0.38 -> height 0.59.
%      xlabel sits just below axes bottom (MATLAB default ~0.03 norm units).
%      Legend top placed at 0.28 norm -> gap of ~0.07 norm below xlabel.
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

h_CRB = semilogy(ax, snr_B4, crb, '--', ...
    'Color', c_CRB, 'LineWidth', lw_CRB);

hold(ax, 'off');

% =========================================================================
%  6.  AXES FORMATTING
% =========================================================================
ax.XLim    = [-5 25];
ax.XTick   = -5:5:25;
ax.YScale  = 'log';
ax.YLim    = [0.05 20];

ax.FontSize = 8;
ax.FontName = 'Times New Roman';
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'off';
ax.YMinorTick = 'on';
ax.Box        = 'on';
grid(ax, 'on');
ax.GridAlpha     = 0.25;
ax.GridLineStyle = ':';

xlabel(ax, 'SNR [dB]',       'Interpreter', 'latex', 'FontSize', 9);
ylabel(ax, 'RMSE$_r$ [m]',   'Interpreter', 'latex', 'FontSize', 9);

% =========================================================================
%  7.  LEGEND
%      Alignment strategy:
%        1. Render legend at natural size with drawnow.
%        2. Read the axes Position (normalised) to get the exact left edge
%           and width of the plot area (i.e. the x-axis span).
%        3. Force the legend to the same left edge and width as the axes,
%           so its left and right borders are flush with the x-axis ends.
%        4. Pin the legend top ~0.05 norm units below the axes bottom
%           (provides a clear gap below the x-axis label).
% =========================================================================
leg = legend(ax, ...
    [h_B1, h_B2, h_B4, h_B5, h_CRB], ...
    {'WB-BPD (full)', ...
     'WB-P-SOMP (comp.)', ...
     'WB-CL-KL (comp.) [prop.]', ...
     'WB-DL-OMP (full)', ...
     'CRB (comp., fixed $r$)'}, ...
    'Interpreter', 'latex', ...
    'FontSize',    7, ...
    'NumColumns',  2, ...
    'Box',         'on', ...
    'Units',       'normalized');

% Force a draw so all positions are finalised
drawnow;

% Read axes extent (normalised figure units)
ax_pos   = ax.Position;          % [left bottom width height]
ax_left  = ax_pos(1);
ax_width = ax_pos(3);
ax_bot   = ax_pos(2);

% Natural legend height (keep content-fitted; only override width/position)
leg_h = leg.Position(4);

% Gap: legend top sits 0.05 norm units below the axes bottom edge
gap        = 0.1;
leg_bottom = ax_bot - gap - leg_h;
leg_bottom = max(leg_bottom, 0.01);  % clamp so legend stays inside figure

leg.Position = [ax_left, leg_bottom, ax_width, leg_h];

% =========================================================================
%  8.  EXPORT
% =========================================================================
out_pdf = 'fig8a_rmse_r_snr.pdf';
print(fig, '-dpdf', '-painters', out_pdf);
fprintf('fig8a: exported to %s\n', out_pdf);

end  % function fig8a_rmse_r_snr
