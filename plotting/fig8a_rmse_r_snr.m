function fig8a_rmse_r_snr()
%FIG8A_RMSE_R_SNR  Publication-quality plot of range RMSE vs SNR (Fig 8a).
%
%  Sprint B revision: re-points to PR-11 CSV (r_hi_fac=0.20, N_MC=600).
%  x-axis trimmed to [-5, 17.5] dB; SNR>=20 excluded (grid-edge outliers
%  per R4/L56).
%
%  Source CSV : mc_snr_sweep_20260528_204613.csv
%               (65 rows x 24 cols; 13 SNR pts at -5:2.5:25 dB; N_MC=600)
%               Located in: src/results/regime_probe_rhi020/
%  Output PDF : fig8a_rmse_r_snr.pdf  (single-column, 8.8 cm x 8.5 cm)
%
%  Locked anchors (PR-11, CSV-verified):
%    B4 RMSE_r @ SNR=10 dB = 0.019821 m
%    B4 CRB_r  @ SNR=10 dB = 0.019909 m
%    B4/CRB = 0.996 (CRB-efficient)
%    B4 conv_pct @ SNR=10 dB = 73.0%
%    fail_rate = 0.0% across [-5, 17.5] dB
%
%  ASCII compliance: no Unicode, no curly quotes, no Greek in strings.
%  Author  : R. V. Senyuva (Maltepe University)
%  Date    : June 2026
%  Version : v2.0  (Sprint B -- PR-11 source, x-axis [-5,17.5])

% =========================================================================
%  1.  LOCATE CSV
% =========================================================================
csv_name   = 'mc_snr_sweep_20260528_204613.csv';
script_dir = fileparts(mfilename('fullpath'));

% Primary: src/results/regime_probe_rhi020/
csv_path = fullfile(script_dir, '..', 'results', 'regime_probe_rhi020', csv_name);
if ~isfile(csv_path)
    csv_path = fullfile(script_dir, csv_name);
end
assert(isfile(csv_path), ...
    ['fig8a: CSV not found. Last tried: ' csv_path]);
fprintf('fig8a: reading %s\n', csv_path);

% =========================================================================
%  2.  READ CSV AND ANCHOR VERIFICATION
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

% Anchor verification
idx10 = abs(snr_B4 - 10) < 1e-9;
v_rmse = rmse_B4(idx10);  v_crb = crb(idx10);
fprintf('fig8a: B4 RMSE_r @ SNR=10: %.6f m  (anchor 0.019821) %s\n', ...
    v_rmse, ok_str(abs(v_rmse - 0.019821) < 1e-4));
fprintf('fig8a: B4 CRB_r  @ SNR=10: %.6f m  (anchor 0.019909) %s\n', ...
    v_crb,  ok_str(abs(v_crb  - 0.019909) < 1e-4));

% Apply x-axis window [-5, 17.5] (exclude SNR >= 20 per R4/L56)
win_mask = snr_B4 <= 17.5 + 1e-9;
snr_B1   = snr_B1(win_mask);   rmse_B1 = rmse_B1(win_mask);
snr_B2   = snr_B2(win_mask);   rmse_B2 = rmse_B2(win_mask);
snr_B4   = snr_B4(win_mask);   rmse_B4 = rmse_B4(win_mask);
snr_B5   = snr_B5(win_mask);   rmse_B5 = rmse_B5(win_mask);
crb      = crb(win_mask);

n_pts = numel(snr_B4);
fprintf('fig8a: %d SNR points in window [-5, 17.5] dB.\n', n_pts);

% =========================================================================
%  3.  STYLE
% =========================================================================
c_B1  = [0.0000 0.4470 0.7410];   % blue
c_B2  = [0.8500 0.3250 0.0980];   % orange
c_B4  = [0.6350 0.0780 0.1840];   % crimson (proposed)
c_B5  = [0.4660 0.6740 0.1880];   % green
c_CRB = [0.0000 0.0000 0.0000];   % black

lw_std = 1.5;  lw_B4 = 2.0;  lw_CRB = 1.0;
ms_std = 8;    ms_B4 = 9;
mi = 1:3:n_pts;   % markers every 3rd point

% =========================================================================
%  4.  FIGURE SETUP (8.8 cm x 8.5 cm; legend below axes)
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
ax.XLim    = [-5 17.5];
ax.XTick   = -5:2.5:17.5;
ax.YScale  = 'log';
ax.YLim    = [0.01 3];

ax.FontSize = 8;
ax.FontName = 'Times New Roman';
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'off';
ax.YMinorTick = 'on';
ax.Box        = 'on';
grid(ax, 'on');
ax.GridAlpha     = 0.25;
ax.GridLineStyle = ':';

xlabel(ax, 'SNR [dB]',     'Interpreter', 'latex', 'FontSize', 9);
ylabel(ax, 'RMSE$_r$ [m]', 'Interpreter', 'latex', 'FontSize', 9);

% =========================================================================
%  7.  LEGEND (2-column, below axes, aligned with x-axis span)
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

ax_pos     = ax.Position;
leg_h      = leg.Position(4);
gap        = 0.10;
leg_bottom = ax_pos(2) - gap - leg_h;
leg_bottom = max(leg_bottom, 0.01);
leg.Position = [ax_pos(1), leg_bottom, ax_pos(3), leg_h];

% =========================================================================
%  8.  EXPORT
% =========================================================================
out_pdf = fullfile(script_dir, 'fig8a_rmse_r_snr.pdf');
print(fig, '-dpdf', '-painters', out_pdf);
fprintf('fig8a: exported to %s\n', out_pdf);

end  % function fig8a_rmse_r_snr

function s = ok_str(flag)
if flag, s = 'OK'; else, s = '*** ANCHOR MISMATCH ***'; end
end
