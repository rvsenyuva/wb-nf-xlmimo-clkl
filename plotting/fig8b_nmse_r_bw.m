function fig8b_nmse_r_bw()
%FIG8B_NMSE_R_BW  Publication-quality plot of range NMSE vs bandwidth (Fig 8b).
%
%  Sprint B revision: single-source PR-12 CSV (r_hi_fac=0.20, N_MC=600, 8 BW pts).
%  The prior 2-CSV merge (PR-6 + PR-9) is SUPERSEDED; use only PR-12.
%
%  Source CSV : mc_bandwidth_sweep_20260529_210446.csv
%               (40 rows x 24 cols; 8 BW pts {100:100:800} MHz; N_MC=600)
%               Located in: src/results/paperC_phase35_fig8b_bw_r020_N600/
%  Output PDF : fig8b_nmse_r_bw.pdf  (single-column, 8.8 cm x 8.5 cm)
%
%  Locked anchors (PR-12, CSV-verified):
%    B4 RMSE_r @ B=400 MHz = 0.019821 m
%    B4 CRB_r  @ B=400 MHz = 0.019909 m;  B4/CRB = 1.00
%    B4 NMSE_r @ B=400 MHz = -43.16 dB
%    CRB data-diversity slope 100->400 MHz = 6.14 dB (Prop1: 6.02)
%    B4 conv_pct: 99.8% @ B=100 -> 47.8% @ B=800; fail_rate=0% throughout
%
%  NMSE normalisation:
%    NMSE_r_dB is read directly from the CSV (already normalised by run_monte_carlo).
%    CRB converted to NMSE using same formula: NMSE_CRB = 20*log10(CRB_r / r_mean)
%    where r_mean = (r_lo + r_hi)/2 = (0.20*0.05 + 0.20)*r_RD/2 is NOT used here;
%    instead NMSE_r_dB column from B4 rows is plotted directly. CRB curve uses
%    the crb_r_m column from B4 rows + the r_mean of the locked scene.
%    r_lo = 0.05*r_hi_fac*r_RD = 0.05*0.20*21.2625 = 0.21263 m (r_hi_fac*r_RD*0.05)
%    r_hi = r_hi_fac*r_RD = 0.20*21.2625 = 4.2525 m
%    r_mean = (r_lo + r_hi)/2 = (0.21263 + 4.2525)/2 = 2.2326 m
%
%  ASCII compliance: no Unicode, no curly quotes, no Greek in strings.
%  Author  : R. V. Senyuva (Maltepe University)
%  Date    : June 2026
%  Version : v2.0  (Sprint B -- single PR-12 CSV, 8-pt, r_hi_fac=0.20)

% =========================================================================
%  1.  LOCATE CSV
% =========================================================================
csv_name   = 'mc_bandwidth_sweep_20260529_210446.csv';
script_dir = fileparts(mfilename('fullpath'));

% Primary: src/results/paperC_phase35_fig8b_bw_r020_N600/
csv_path = fullfile(script_dir, '..', 'results', ...
    'paperC_phase35_fig8b_bw_r020_N600', csv_name);
if ~isfile(csv_path)
    % Fallback: same folder as script
    csv_path = fullfile(script_dir, csv_name);
end
assert(isfile(csv_path), ...
    ['fig8b: CSV not found. Last tried: ' csv_path]);
fprintf('fig8b: reading %s\n', csv_path);

% =========================================================================
%  2.  READ CSV AND ANCHOR VERIFICATION
% =========================================================================
T = readtable(csv_path, 'TextType', 'string');

    function [bv, nv, cv] = extract(T, mname)
        rows = T(strcmp(T.method, mname), :);
        rows = sortrows(rows, 'B_hz');
        bv = rows.B_hz / 1e6;    % Hz -> MHz
        nv = rows.NMSE_r_dB;
        cv = rows.crb_r_m;
    end

[bw_B1, nmse_B1, ~  ] = extract(T, 'WB-BPD');
[bw_B2, nmse_B2, ~  ] = extract(T, 'WB-P-SOMP');
[bw_B4, nmse_B4, crb] = extract(T, 'WB-CL-KL');
[bw_B5, nmse_B5, ~  ] = extract(T, 'WB-DL-OMP');

% Anchor verification at B=400 MHz
idx400 = abs(bw_B4 - 400) < 1e-6;
v_nmse = nmse_B4(idx400);  v_crb = crb(idx400);
fprintf('fig8b: B4 NMSE_r @ B=400: %.4f dB  (anchor -43.16) %s\n', ...
    v_nmse, ok_str(abs(v_nmse - (-43.1561)) < 0.5));
fprintf('fig8b: B4 CRB_r  @ B=400: %.6f m   (anchor 0.019909) %s\n', ...
    v_crb,  ok_str(abs(v_crb - 0.019909) < 1e-4));

% CRB NMSE curve (B4 rows only; normalise by scene-average r)
% r_hi_fac=0.20: r_lo = 0.05*r_hi_fac*r_RD, r_hi = r_hi_fac*r_RD
r_hi_fac = 0.20;  r_RD = 21.2625;  r_hi_frac_lo = 0.05;
r_lo   = r_hi_frac_lo * r_hi_fac * r_RD;   % 0.21263 m
r_hi   = r_hi_fac * r_RD;                   % 4.2525 m
r_mean = (r_lo + r_hi) / 2;                 % 2.23257 m
fprintf('fig8b: r_mean for CRB NMSE normalisation = %.5f m\n', r_mean);
crb_nmse = 20 * log10(crb / r_mean);        % CRB as NMSE_r (dB)

n_pts = numel(bw_B4);
fprintf('fig8b: %d BW points loaded.\n', n_pts);

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
mi = 1:n_pts;   % all 8 BW points get markers

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
ax.XTick   = 100:100:800;
ax.YScale  = 'linear';

% Auto y-limits with margin
all_nmse = [nmse_B1; nmse_B2; nmse_B4; nmse_B5; crb_nmse];
all_nmse = all_nmse(isfinite(all_nmse));
y_lo = floor(min(all_nmse) / 5) * 5 - 5;
y_hi = ceil( max(all_nmse) / 5) * 5 + 5;
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
out_pdf = fullfile(script_dir, 'fig8b_nmse_r_bw.pdf');
print(fig, '-dpdf', '-painters', out_pdf);
fprintf('fig8b: exported to %s\n', out_pdf);

end  % function fig8b_nmse_r_bw

function s = ok_str(flag)
if flag, s = 'OK'; else, s = '*** ANCHOR MISMATCH ***'; end
end
