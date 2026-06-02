function fig10_kl_convergence()
%FIG10_KL_CONVERGENCE  Publication-quality KL-objective vs iteration (Fig 10).
%
%  Sprint B revision v1.1: plots NORMALIZED objective delta_L(t) = L(t)-L(1)
%  (per-trace offset to zero at iter=1, then median across traces).
%  This matches the Paper B TCOM Fig 9 convention and makes the ~12-23 unit
%  descent visible against the ~800-unit absolute baseline.
%
%  Source CSVs (PR-13, Sprint A -- r=3.0 m, r_hi_fac=0.20, N_MC=50):
%    mc_convergence_Lhist_SNR+0dB_20260531_215757.csv   (9431 rows)
%    mc_convergence_Lhist_SNR+5dB_20260531_215757.csv   (9642 rows)
%    mc_convergence_Lhist_SNR+10dB_20260531_215757.csv  (9621 rows)
%    mc_convergence_Lhist_SNR+15dB_20260531_215757.csv  (9663 rows)
%  Located in: src/results/paperC_sprint_a_fig10_convergence_r3m/
%
%  Fixed convergence scene: theta=35 deg, r=3.0 m (inside locked box
%  [0.2126, 4.2525] m; Fresnel kappa*m_bar^2 = 0.9157 rad > 0.25, PASS).
%  Scenario B (wrong-basin): 47-48/50 traces hit max_iter=200.
%  L is monotone decreasing in ALL trials per Armijo guarantee (Prop 3).
%  delta_L ranges: -11.8 (SNR=0), -23.3 (SNR=5), -17.5 (SNR=10), -21.4 (SNR=15).
%
%  Colors: cool-to-warm sequential (SNR 0->15 dB: blue->cyan->orange->red).
%
%  Output PDF: fig10_kl_convergence.pdf  (single-column, 8.8 cm x 6.5 cm)
%
%  ASCII compliance: no Unicode, no curly quotes, no Greek in strings.
%  Author  : R. V. Senyuva (Maltepe University)
%  Date    : June 2026
%  Version : v1.1  (Sprint B -- normalized delta_L, per-trace offset)

% =========================================================================
%  0.  SNR levels and file names (PR-13 timestamp)
% =========================================================================
snr_levels  = [0, 5, 10, 15];
ts_str      = '20260531_215757';
csv_basename = 'mc_convergence_Lhist_SNR%+.0fdB_%s.csv';

% Colors: blue, cyan-blue, orange, red (cool -> warm = low -> high SNR)
clr_snr = [0.3000 0.3000 0.9500;    % 0 dB:  deep blue
           0.3000 0.7000 0.9500;     % 5 dB:  cyan-blue
           0.9500 0.5500 0.3000;     % 10 dB: orange
           0.8500 0.2000 0.2000];    % 15 dB: red

% =========================================================================
%  1.  LOCATE CSVs
% =========================================================================
script_dir = fileparts(mfilename('fullpath'));
csv_dir_primary = fullfile(script_dir, '..', 'results', ...
    'paperC_sprint_a_fig10_convergence_r3m');
csv_dir_fallback = script_dir;

% =========================================================================
%  2.  LOAD AND COMPUTE MEDIAN NORMALIZED TRAJECTORIES
%      delta_L(t) = L(t) - L(iter=1) computed per trace, then median.
% =========================================================================
max_iter = 200;
med_traj = zeros(max_iter, 4);   % median delta_L per iteration
n_traces = zeros(1, 4);

for si = 1:4
    snr_val  = snr_levels(si);
    csv_name = sprintf(csv_basename, snr_val, ts_str);
    csv_path = fullfile(csv_dir_primary, csv_name);
    if ~isfile(csv_path)
        csv_path = fullfile(csv_dir_fallback, csv_name);
    end
    assert(isfile(csv_path), sprintf('fig10: CSV not found: %s', csv_path));
    fprintf('fig10: loading SNR=%+d dB: %s\n', snr_val, csv_name);

    T = readtable(csv_path, 'TextType', 'string');
    iter_vals = T.iter_idx;
    L_vals    = T.L_val;
    mc_ids    = unique(T.mc_idx);
    n_traces(si) = numel(mc_ids);

    % Build L_by_iter matrix (iter x mc), NaN-filled
    L_by_iter = NaN(max_iter, n_traces(si));
    for mc_ii = 1:numel(mc_ids)
        mc_mask = T.mc_idx == mc_ids(mc_ii);
        iter_mc = iter_vals(mc_mask);
        L_mc    = L_vals(mc_mask);
        for it = 1:numel(iter_mc)
            if iter_mc(it) >= 1 && iter_mc(it) <= max_iter
                L_by_iter(iter_mc(it), mc_ii) = L_mc(it);
            end
        end
    end

    % Per-trace normalization: delta_L(t) = L(t) - L(iter=1)
    % L(iter=1) is the value at the first recorded iteration for each trace
    L_init_per_trace = L_by_iter(1, :);   % 1 x n_traces
    dL_by_iter = L_by_iter - L_init_per_trace;   % broadcast: each col -= scalar

    % For early-stopped traces (iter < max_iter), carry forward the last value
    % (they converged, so delta_L is flat at its final value -- correct behavior)
    for mc_ii = 1:n_traces(si)
        last_valid = find(~isnan(dL_by_iter(:, mc_ii)), 1, 'last');
        if ~isempty(last_valid) && last_valid < max_iter
            dL_by_iter(last_valid+1:end, mc_ii) = dL_by_iter(last_valid, mc_ii);
        end
    end

    % Median across traces (NaN-safe via nanmedian equivalent)
    for it = 1:max_iter
        row = dL_by_iter(it, :);
        row_valid = row(~isnan(row));
        if ~isempty(row_valid)
            med_traj(it, si) = median(row_valid);
        else
            med_traj(it, si) = 0;
        end
    end

    fprintf('fig10:   traces=%d  delta_L: init=0 -> final=%.4f\n', ...
        n_traces(si), med_traj(max_iter, si));
end

iter_axis = (1:max_iter)';

% =========================================================================
%  3.  FIGURE SETUP  (8.8 cm x 6.5 cm; legend inside axes)
% =========================================================================
fig_w = 8.8;
fig_h = 6.5;

fig = figure('Units', 'centimeters', ...
             'Position', [2 2 fig_w fig_h], ...
             'PaperUnits', 'centimeters', ...
             'PaperSize',  [fig_w fig_h], ...
             'PaperPosition', [0 0 fig_w fig_h], ...
             'Color', 'w');

ax = axes('Parent', fig, ...
          'Units', 'normalized', ...
          'Position', [0.14 0.15 0.83 0.79]);

% =========================================================================
%  4.  PLOT CURVES
% =========================================================================
hold(ax, 'on');
h = gobjects(4, 1);
lw_line = 1.5;
for si = 1:4
    h(si) = plot(ax, iter_axis, med_traj(:, si), '-', ...
        'Color', clr_snr(si,:), 'LineWidth', lw_line);
end

% Reference line at delta_L = 0 (grey dashed; shows starting level)
plot(ax, [1 max_iter], [0 0], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);

hold(ax, 'off');

% =========================================================================
%  5.  AXES FORMATTING
% =========================================================================
ax.XLim    = [0 200];
ax.XTick   = 0:40:200;
ax.YScale  = 'linear';

% y-limits: cover full descent with margin
y_min = min(med_traj(:)) * 1.15;   % 15% below deepest point
y_max = max(med_traj(:)) + 2;      % small positive margin above 0
ax.YLim = [y_min y_max];

ax.FontSize = 8;
ax.FontName = 'Times New Roman';
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'off';
ax.YMinorTick = 'on';
ax.Box        = 'on';
grid(ax, 'on');
ax.GridAlpha     = 0.25;
ax.GridLineStyle = ':';

xlabel(ax, 'Iteration',                          'Interpreter', 'latex', 'FontSize', 9);
ylabel(ax, '$\Delta\mathcal{L}$ (median)', 'Interpreter', 'latex', 'FontSize', 9);

% =========================================================================
%  6.  LEGEND (inside axes)
% =========================================================================
leg_str = {'SNR $= 0$ dB', 'SNR $= 5$ dB', ...
           'SNR $= 10$ dB', 'SNR $= 15$ dB'};
legend(ax, h, leg_str, ...
    'Interpreter', 'latex', ...
    'FontSize',    7, ...
    'Location',    'southwest', ...
    'Box',         'on');

% =========================================================================
%  7.  EXPORT
% =========================================================================
out_pdf = fullfile(script_dir, 'fig10_kl_convergence.pdf');
print(fig, '-dpdf', '-painters', out_pdf);
fprintf('fig10: exported to %s\n', out_pdf);

end  % function fig10_kl_convergence
