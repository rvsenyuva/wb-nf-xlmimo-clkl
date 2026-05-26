function fig9a_rmse_theta_snr()
%FIG9A_RMSE_THETA_SNR  Publication-quality plot of angle RMSE vs SNR (Fig 9a).
%
%  Produces fig9a_rmse_theta_snr.pdf for Paper C (Phase 3.5 / OJ-COMS).
%
%  Source CSV : mc_snr_sweep_20260526_090431.csv
%               (65 rows x 24 cols; 13 SNR pts at -5:2.5:25 dB; N_MC=200)
%               REUSED from Fig 8a -- no additional MATLAB run needed.
%  Output PDF : fig9a_rmse_theta_snr.pdf  (single-column, 8.8 cm x 8.5 cm)
%
%  Five curves
%    B1  WB-BPD     (full-array, open blue circle,       1.5 pt)
%    B2  WB-P-SOMP  (compressed, filled orange square,   1.5 pt)
%    B4  WB-CL-KL   (compressed, filled crimson diamond, 2.0 pt) [PROPOSED]
%    B5  WB-DL-OMP  (full-array, open green triangle,    1.5 pt)
%    CRB_theta       (black dashed, no markers,           1.0 pt)
%
%  Legend: 2-column, placed outside below x-axis, aligned with axes x-span.
%  No title string.
%
%  ASCII compliance: no Unicode, no curly quotes, no Greek in strings.
%  Author  : R. V. Senyuva (Maltepe University)
%  Date    : 2026-05-26
%  Version : v1.0  (Phase 3.5 Subtask 3.5.4)

% =========================================================================
%  1.  LOCATE CSV  (same file as Fig 8a)
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
    ['fig9a: CSV not found. Last tried: ' csv_path]);
fprintf('fig9a: reading %s\n', csv_path);

% =========================================================================
%  2.  READ CSV
% =========================================================================
T = readtable(csv_path, 'TextType', 'string');

    function [sv, tv, ctv] = extract(T, mname)
        rows = T(strcmp(T.method, mname), :);
        rows = sortrows(rows, 'SNR_dB');
        sv  = rows.SNR_dB;
        tv  = rows.RMSE_theta_deg;
        ctv = rows.crb_theta_deg;
    end

[snr_B1, theta_B1, ~      ] = extract(T, 'WB-BPD');
[snr_B2, theta_B2, ~      ] = extract(T, 'WB-P-SOMP');
[snr_B4, theta_B4, crb_th ] = extract(T, 'WB-CL-KL');
[snr_B5, theta_B5, ~      ] = extract(T, 'WB-DL-OMP');

n_pts = numel(snr_B4);
fprintf('fig9a: %d SNR points loaded.\n', n_pts);
fprintf('fig9a: B4 RMSE_theta floor = %.4f deg (at SNR=%.1f dB)\n', ...
        min(theta_B4), snr_B4(theta_B4 == min(theta_B4)));

% =========================================================================
%  3.  STYLE  (Phase 3.5 Plan Section 2 -- identical to Fig 8a)
% =========================================================================
c_B1  = [0.0000 0.4470 0.7410];
c_B2  = [0.8500 0.3250 0.0980];
c_B4  = [0.6350 0.0780 0.1840];   % crimson
c_B5  = [0.4660 0.6740 0.1880];
c_CRB = [0.0000 0.0000 0.0000];

lw_std = 1.5;  lw_B4 = 2.0;  lw_CRB = 1.0;
ms_std = 8;    ms_B4 = 9;
mi = 1:3:n_pts;   % markers at indices 1, 4, 7, 10, 13

% =========================================================================
%  4.  FIGURE SETUP  (identical dimensions to Fig 8a)
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

h_B1 = semilogy(ax, snr_B1, theta_B1, '-o', ...
    'Color', c_B1, 'LineWidth', lw_std, ...
    'MarkerSize', ms_std, 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', c_B1, 'MarkerIndices', mi);

h_B2 = semilogy(ax, snr_B2, theta_B2, '-s', ...
    'Color', c_B2, 'LineWidth', lw_std, ...
    'MarkerSize', ms_std, 'MarkerFaceColor', c_B2, ...
    'MarkerEdgeColor', c_B2, 'MarkerIndices', mi);

h_B4 = semilogy(ax, snr_B4, theta_B4, '-d', ...
    'Color', c_B4, 'LineWidth', lw_B4, ...
    'MarkerSize', ms_B4, 'MarkerFaceColor', c_B4, ...
    'MarkerEdgeColor', c_B4, 'MarkerIndices', mi);

h_B5 = semilogy(ax, snr_B5, theta_B5, '-^', ...
    'Color', c_B5, 'LineWidth', lw_std, ...
    'MarkerSize', ms_std, 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', c_B5, 'MarkerIndices', mi);

h_CRB = semilogy(ax, snr_B4, crb_th, '--', ...
    'Color', c_CRB, 'LineWidth', lw_CRB);

hold(ax, 'off');

% =========================================================================
%  6.  AXES FORMATTING
% =========================================================================
ax.XLim    = [-5 25];
ax.XTick   = -5:5:25;
ax.YScale  = 'log';
% B4/B1 floor ~0.104 deg; B2 peaks at ~10.5 deg; CRB falls to ~0.0008 deg.
% Set limits to show all curves without wasted space.
ax.YLim    = [5e-4 20];

ax.FontSize = 8;
ax.FontName = 'Times New Roman';
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'off';
ax.YMinorTick = 'on';
ax.Box        = 'on';
grid(ax, 'on');
ax.GridAlpha     = 0.25;
ax.GridLineStyle = ':';

xlabel(ax, 'SNR [dB]',           'Interpreter', 'latex', 'FontSize', 9);
ylabel(ax, 'RMSE$_\theta$ [deg]', 'Interpreter', 'latex', 'FontSize', 9);

% =========================================================================
%  7.  LEGEND  -- 2-column, outside below axes, aligned with x-axis span
% =========================================================================
leg = legend(ax, ...
    [h_B1, h_B2, h_B4, h_B5, h_CRB], ...
    {'WB-BPD (full)', ...
     'WB-P-SOMP (comp.)', ...
     'WB-CL-KL (comp.) [prop.]', ...
     'WB-DL-OMP (full)', ...
     'CRB (comp., fixed $r$, $\theta$)'}, ...
    'Interpreter', 'latex', ...
    'FontSize',    6.5, ...
    'NumColumns',  2, ...
    'Box',         'on', ...
    'Units',       'normalized');

drawnow;

% Align legend left/width with axes x-span; gap = 0.10 norm units below axes
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
out_pdf = 'fig9a_rmse_theta_snr.pdf';
print(fig, '-dpdf', '-painters', out_pdf);
fprintf('fig9a: exported to %s\n', out_pdf);

end  % function fig9a_rmse_theta_snr
