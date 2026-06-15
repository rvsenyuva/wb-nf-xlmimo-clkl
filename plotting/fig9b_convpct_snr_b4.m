function fig9b_convpct_snr_b4()
%FIG9B_CONVPCT_SNR_B4  Publication-quality plot of conv_pct vs SNR (Fig 9b).
%
%  Sprint B revision: re-points to PR-11 CSV (r_hi_fac=0.20, N_MC=600).
%  x-axis [-5, 17.5] dB (same as Figs 8a/9a).
%  REUSES the same CSV as Fig 8a/9a -- no additional MATLAB run needed.
%
%  Source CSV : mc_snr_sweep_20260528_204613.csv
%               Located in: src/results/regime_probe_rhi020/
%
%  Locked anchor (PR-11, CSV-verified):
%    B4 conv_pct @ SNR=10 dB = 73.0%  (N_MC=600; fail_rate=0%)
%
%  FRAMING NOTE (binding per Phase3_5_Plan Sec 3.7):
%    conv_pct is a STRICT-TOLERANCE metric, NOT a success rate.
%    Non-converged trials at SNR=10 still achieve CRB-level RMSE (B4/CRB=0.996).
%    The non-monotonic pattern (100%->100%->84.8%->73%->67.5%->68.2% in window)
%    reflects the optimizer landscape in the strong near-field regime.
%    fail_rate = 0% across the entire [-5, 17.5] dB window.
%
%  ASCII compliance: no Unicode, no curly quotes, no Greek in strings.
%  Author  : R. V. Senyuva (Maltepe University)
%  Date    : June 2026
%  Version : v2.0  (Sprint B -- PR-11 source, x-axis [-5,17.5], strict-tol framing)

% =========================================================================
%  1.  LOCATE CSV  (same file as Figs 8a and 9a)
% =========================================================================
csv_name   = 'mc_snr_sweep_20260528_204613.csv';
script_dir = fileparts(mfilename('fullpath'));

csv_path = fullfile(script_dir, '..', 'results', 'regime_probe_rhi020', csv_name);
if ~isfile(csv_path)
    csv_path = fullfile(script_dir, csv_name);
end
assert(isfile(csv_path), ...
    ['fig9b: CSV not found. Last tried: ' csv_path]);
fprintf('fig9b: reading %s\n', csv_path);

% =========================================================================
%  2.  READ CSV AND ANCHOR VERIFICATION
% =========================================================================
T = readtable(csv_path, 'TextType', 'string');

rows_B4 = T(strcmp(T.method, 'WB-CL-KL'), :);
rows_B4 = sortrows(rows_B4, 'SNR_dB');

snr_all  = rows_B4.SNR_dB;
conv_all = rows_B4.clkl_conv_pct;
fail_all = rows_B4.fail_rate_pct;

% Anchor verification
idx10   = abs(snr_all - 10) < 1e-9;
v_conv  = conv_all(idx10);
v_fail  = fail_all(idx10);
fprintf('fig9b: B4 conv_pct @ SNR=10: %.2f%%  (anchor 73.0) %s\n', ...
    v_conv, ok_str(abs(v_conv - 73.0) < 1.0));
fprintf('fig9b: B4 fail_rate @ SNR=10: %.2f%%  (expected 0.0) %s\n', ...
    v_fail, ok_str(v_fail < 0.5));

% Apply x-axis window [-5, 17.5]
win_mask = snr_all <= 17.5 + 1e-9;
snr_plt  = snr_all(win_mask);
conv_plt = conv_all(win_mask);
fail_plt = fail_all(win_mask);

n_pts = numel(snr_plt);
fprintf('fig9b: %d SNR points in window [-5, 17.5] dB.\n', n_pts);
fprintf('fig9b: conv_pct pattern: ');
fprintf('%.1f%% ', conv_plt);
fprintf('\n');
fprintf('fig9b: fail_rate (all should be 0%%): ');
fprintf('%.1f%% ', fail_plt);
fprintf('\n');

% =========================================================================
%  3.  STYLE  (B4 crimson only; diamond marker as in Fig 8a)
% =========================================================================
c_B4 = [0.6350 0.0780 0.1840];   % crimson
lw_B4 = 2.0;
ms_B4 = 9;
mi = 1:n_pts;   % all points get markers (10 pts -- not overcrowded)

% =========================================================================
%  4.  FIGURE SETUP (8.8 cm x 6.5 cm; self-contained legend inside axes)
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
%  5.  PLOT CURVE
% =========================================================================
hold(ax, 'on');

h_conv = plot(ax, snr_plt, conv_plt, '-d', ...
    'Color', c_B4, 'LineWidth', lw_B4, ...
    'MarkerSize', ms_B4, 'MarkerFaceColor', c_B4, ...
    'MarkerEdgeColor', c_B4, 'MarkerIndices', mi);

hold(ax, 'off');

% =========================================================================
%  6.  AXES FORMATTING
% =========================================================================
ax.XLim    = [-5 17.5];
ax.XTick   = -5:2.5:17.5;
ax.YScale  = 'linear';
ax.YLim    = [0 110];
ax.YTick   = 0:20:100;

ax.FontSize = 8;
ax.FontName = 'Times New Roman';
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.Box        = 'on';
grid(ax, 'on');
ax.GridAlpha     = 0.25;
ax.GridLineStyle = ':';

xlabel(ax, 'SNR [dB]',          'Interpreter', 'latex', 'FontSize', 9);
ylabel(ax, 'Conv. rate [\%]',   'Interpreter', 'latex', 'FontSize', 9);

% =========================================================================
%  7.  LEGEND (single entry; inside axes -- no legend box needed but add for
%      consistency with companion Fig 9a and manuscript style).
%  Note text references the strict-tolerance framing explicitly.
% =========================================================================
lgd = legend(ax, h_conv, ...
    'WB-CL-KL [prop.] (strict tol.)', ...
    'Interpreter', 'latex', ...
    'FontSize', 7, ...
    'Orientation', 'horizontal', ...
    'Box', 'on');
lgd.Units = 'normalized';
%lgd.Position = [0.15, 0.01, 0.70, 0.07];
ax_pos     = ax.Position;
lgd_h      = lgd.Position(4);
gap        = 0.10;
lgd_bottom = ax_pos(2) - gap - lgd_h;
lgd_bottom = max(lgd_bottom, 0.01);
lgd.Position = [ax_pos(1), lgd_bottom, ax_pos(3), lgd_h];

% =========================================================================
%  8.  EXPORT
% =========================================================================
out_pdf = fullfile(script_dir, 'fig9b_convpct_snr_b4.pdf');
print(fig, '-dpdf', '-painters', out_pdf);
fprintf('fig9b: exported to %s\n', out_pdf);

end  % function fig9b_convpct_snr_b4

function s = ok_str(flag)
if flag, s = 'OK'; else, s = '*** ANCHOR MISMATCH ***'; end
end
