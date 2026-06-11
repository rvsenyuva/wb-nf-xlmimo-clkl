%FIG11_RIVIELLO_SNR_AXIS  Plot Fig. 11: two-panel stacked figure.
%  Top panel  : WB-CL-KL RMSE_r and CRB_r vs. SNR (log y-axis).
%               Legend inside top panel (2 entries only).
%  Bottom panel: 3GPP UMi per-UT SNR histogram (linear y-axis).
%               Labelled via text annotation, not legend entry.
%  Both panels share the same x-axis (SNR in dB).
%
%  INPUTS:
%    SNR_SWEEP_CSV    : results/riviello_snr/mc_snr_sweep_20260610_203626.csv
%    DISTRIBUTION_CSV : results/riviello_snr/riviello_snr_distribution_2000.csv
%  OUTPUT:
%    plotting/fig11_riviello_snr_axis.pdf
%
%  Author   : R. V. Senyuva (Maltepe University)
%  Date     : June 2026
%  Ref      : Paper C Phase 4, Session 4.9 (D1 robustness experiment).

% =========================================================================
%  USER PATHS
% =========================================================================
SNR_SWEEP_CSV    = fullfile('results', 'riviello_snr', ...
    'mc_snr_sweep_20260610_203626.csv');
DISTRIBUTION_CSV = fullfile('results', 'riviello_snr', ...
    'riviello_snr_distribution_2000.csv');
OUT_PDF = fullfile('plotting', 'fig11_riviello_snr_axis.pdf');

% =========================================================================
%  LOAD DATA
% =========================================================================
T    = readtable(SNR_SWEEP_CSV);
T_B4 = T(strcmpi(T.method, 'WB-CL-KL'), :);
[snr_pts, sidx] = sort(T_B4.SNR_dB);
rmse_r = T_B4.RMSE_r_m(sidx);
crb_r  = T_B4.crb_r_m(sidx);

snr_dist   = readmatrix(DISTRIBUTION_CSV);
snr_dist   = snr_dist(:);
median_snr = median(snr_dist);

[~, idx_med] = min(abs(snr_pts - median_snr));
b4crb_med    = rmse_r(idx_med) / crb_r(idx_med);

fprintf('[fig11] Median SNR=%.2f dB  nearest sweep=%.0f dB  B4/CRB=%.4f\n', ...
    median_snr, snr_pts(idx_med), b4crb_med);

% Histogram: 5 dB bins over plot range
bin_edges = (-20 : 5 : 40);
bin_ctrs  = bin_edges(1:end-1) + 2.5;
counts    = histcounts(snr_dist, bin_edges);
density   = counts / (numel(snr_dist) * 5);   % fraction per dB

% =========================================================================
%  LAYOUT
% =========================================================================
x_lo = -20;   x_hi = 35;
y_lo = 4e-4;  y_hi = 2e0;

fig_w_cm = 8.8;
fig_h_cm = 10.5;

% Positions: top panel taller, bottom panel shorter, no below-axes space needed
ax_main_pos = [0.14  0.38  0.83  0.58];   % top: curves
ax_hist_pos = [0.14  0.11  0.83  0.20];   % bottom: histogram

fig = figure('Units', 'centimeters', ...
    'Position', [2 2 fig_w_cm fig_h_cm], ...
    'Color', 'w');

% =========================================================================
%  TOP PANEL
% =========================================================================
ax1 = axes('Parent', fig, 'Position', ax_main_pos);

h_crb = semilogy(ax1, snr_pts, crb_r, 'k-o', ...
    'LineWidth', 1.4, 'MarkerSize', 4, 'MarkerFaceColor', 'k', ...
    'DisplayName', ...
    '$\sqrt{\mathrm{CRB}_r}$ (nominal geom., $N_\mathrm{RF}\!=\!8$)');
hold(ax1, 'on');

h_b4 = semilogy(ax1, snr_pts, rmse_r, 'b-s', ...
    'LineWidth', 1.4, 'MarkerSize', 4, ...
    'MarkerFaceColor', [0.15 0.35 0.75], ...
    'DisplayName', 'WB-CL-KL $\sqrt{\mathrm{RMSE}_r}$');

% Median vertical line
semilogy(ax1, [median_snr median_snr], [y_lo*1.5 y_hi*0.85], ...
    'k--', 'LineWidth', 0.9, 'HandleVisibility', 'off');

% B4/CRB annotation: right of the median line, near the B4 curve
text(ax1, median_snr + 0.5, rmse_r(idx_med) * 2.2, ...
    sprintf('B4/CRB$\\!=\\!%.3f$', b4crb_med), ...
    'Interpreter', 'latex', 'FontSize', 7, ...
    'Color', [0.10 0.10 0.65], 'HorizontalAlignment', 'left');

% Median SNR label near below of vertical line
text(ax1, median_snr, y_lo * 0.70, ...
    sprintf('$%.1f$~dB', median_snr), ...
    'Interpreter', 'latex', 'FontSize', 6.5, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

ax1.XLim   = [x_lo, x_hi];
ax1.YLim   = [y_lo, y_hi];
ax1.YTick  = [1e-3 1e-2 1e-1 1e0];
ax1.XTick  = x_lo : 5 : x_hi;
ax1.XTickLabel = {};
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 8;
grid(ax1, 'on');
ax1.GridAlpha = 0.25;

ylabel(ax1, '$\sqrt{\mathrm{RMSE}_r}$, $\sqrt{\mathrm{CRB}_r}$~(m)', ...
    'Interpreter', 'latex', 'FontSize', 8.5);

% Legend: 2 entries only, inside top panel, northeast corner
legend(ax1, [h_crb, h_b4], ...
    'Location',    'northeast', ...
    'Interpreter', 'latex', ...
    'FontSize',    6.5, ...
    'Box',         'on');

hold(ax1, 'off');

% =========================================================================
%  BOTTOM PANEL
% =========================================================================
ax2 = axes('Parent', fig, 'Position', ax_hist_pos);

bar(ax2, bin_ctrs, density, 1.0, ...
    'FaceColor', [0.65 0.80 0.92], ...
    'EdgeColor', [0.40 0.60 0.78], ...
    'FaceAlpha', 0.85);
hold(ax2, 'on');

% Median vertical line
plot(ax2, [median_snr median_snr], [0, max(density) * 1.15], ...
    'k--', 'LineWidth', 0.9);

% Text label in upper-right of histogram panel (replaces legend entry)
text(ax2, x_hi - 0.5, max(density) * 1.20, ...
    '3GPP UMi SNR density', ...
    'Interpreter', 'latex', 'FontSize', 7, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'Color', [0.20 0.20 0.20]);

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
%  SAVE
% =========================================================================
out_dir_fig = fileparts(OUT_PDF);
if ~isempty(out_dir_fig) && ~exist(out_dir_fig, 'dir')
    mkdir(out_dir_fig);
end
set(fig, 'PaperUnits', 'centimeters', ...
    'PaperSize',     [fig_w_cm fig_h_cm], ...
    'PaperPosition', [0 0 fig_w_cm fig_h_cm]);
print(fig, OUT_PDF, '-dpdf', '-r300');
fprintf('[fig11] Saved: %s\n', OUT_PDF);
close(fig);
