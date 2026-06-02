function fig7_crb_geodiv_bw()
%FIG7_CRB_GEODIV_BW  Publication-quality plot of geometric diversity CRB (Fig 7).
%
%  sqrt-CRB on range vs OFDM bandwidth B, for M in {32, 64, 128}.
%  d=1, N_RF=8, SNR=10 dB, theta=35 deg, r=5.0 m.
%  Demonstrates: (1) geometric diversity gain (~9 dB over 50->400 MHz),
%  (2) M-doubling gap (~12 dB) dominates the bandwidth gain.
%
%  Source CSV : wb_crb_geodiv_sweep_20260520_184631.csv
%               (18 rows x 7 cols; 3 M x 6 BW pts)
%               Located in: src/results/ (repo root results folder)
%  Output PDF : fig7_crb_geodiv_bw.pdf  (single-column, 8.8 cm x 6.5 cm)
%
%  Locked anchors (CSV-verified):
%    M=32  geodiv gain 50->400 MHz = 9.047 dB
%    M=64  geodiv gain 50->400 MHz = 9.125 dB
%    M=128 geodiv gain 50->400 MHz = 9.174 dB
%    M=32 -> M=64  gap @ B=400 MHz = 11.733 dB
%    M=64 -> M=128 gap @ B=400 MHz = 12.629 dB
%
%  Color scheme (M-specific; distinct from method palette):
%    M=32  : teal        [0.0000 0.5490 0.5490]  solid, circle
%    M=64  : purple      [0.5020 0.0000 0.5020]  solid, square
%    M=128 : dark orange [0.8500 0.4000 0.0000]  solid, triangle-up
%
%  Both axes log scale; x-axis B in MHz.
%  Markers at all 6 BW points (6 pts -- not overcrowded).
%  Legend inside axes (northeast).
%
%  ASCII compliance: no Unicode, no curly quotes, no Greek in strings.
%  Author  : R. V. Senyuva (Maltepe University)
%  Date    : June 2026
%  Version : v1.0  (Sprint B)

% =========================================================================
%  0.  CSV name and path
% =========================================================================
csv_name   = 'wb_crb_geodiv_sweep_20260520_184631.csv';
script_dir = fileparts(mfilename('fullpath'));

% Primary: results/paperC_phase3_fig7_crb_geodiv_v1/
% (script lives at plotting/; CSV lives at results/paperC_.../; one level up)
csv_path = fullfile(script_dir, '..', 'results', ...
                    'paperC_phase3_fig7_crb_geodiv_v1', csv_name);
if ~isfile(csv_path)
    % Fallback: same directory as script (manual copy)
    csv_path = fullfile(script_dir, csv_name);
end
assert(isfile(csv_path), ['fig7: CSV not found. Last tried: ' csv_path]);
fprintf('fig7: reading %s\n', csv_path);

% =========================================================================
%  1.  Load and sort data
% =========================================================================
T = readtable(csv_path);

M_vec = [32, 64, 128];
crb   = struct();

for ii = 1:numel(M_vec)
    mask          = T.M == M_vec(ii);
    sub           = T(mask, :);
    [~, idx]      = sort(sub.B_hz);
    crb(ii).B_mhz = sub.B_hz(idx) / 1e6;
    crb(ii).r     = sub.crb_r_m(idx);
    crb(ii).rRD   = sub.r_RD_m(idx(1));
end

% Anchor verification
tol_dB = 0.05;
for ii = 1:3
    idx50  = abs(crb(ii).B_mhz -  50) < 1e-6;
    idx400 = abs(crb(ii).B_mhz - 400) < 1e-6;
    idx128_lo = abs(crb(ii).B_mhz - 400) < 1e-6;
    gain_50_400 = 20*log10(crb(ii).r(idx50) / crb(ii).r(idx400));
    ref_gain = [9.047, 9.125, 9.174];
    fprintf('fig7: M=%3d  geodiv gain 50->400 MHz = %.3f dB  (ref %.3f)  %s\n', ...
        M_vec(ii), gain_50_400, ref_gain(ii), ...
        ok_str(abs(gain_50_400 - ref_gain(ii)) < tol_dB));
end

% M-doubling gap at B=400
for ii = 1:2
    idx400_lo = abs(crb(ii).B_mhz   - 400) < 1e-6;
    idx400_hi = abs(crb(ii+1).B_mhz - 400) < 1e-6;
    gap = 20*log10(crb(ii).r(idx400_lo) / crb(ii+1).r(idx400_hi));
    ref_gap = [11.733, 12.629];
    fprintf('fig7: M=%3d->%3d gap @ 400 MHz = %.3f dB  (ref %.3f)  %s\n', ...
        M_vec(ii), M_vec(ii+1), gap, ref_gap(ii), ...
        ok_str(abs(gap - ref_gap(ii)) < tol_dB));
end
fprintf('\n');

% =========================================================================
%  2.  Style
% =========================================================================
clr32  = [0.0000 0.5490 0.5490];   % teal
clr64  = [0.5020 0.0000 0.5020];   % purple
clr128 = [0.8500 0.4000 0.0000];   % dark orange

lw = 1.5;
msz = 7;
mi  = 1:6;   % all 6 BW points get markers

% =========================================================================
%  3.  Figure setup (IEEEtran single-column: 8.8 cm x 6.5 cm)
% =========================================================================
fig = figure('Units', 'centimeters', 'Position', [2 2 8.8 6.5]);
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [8.8 6.5], ...
         'PaperPosition', [0 0 8.8 6.5], 'Color', 'w');

ax = axes(fig);
set(ax, 'XScale', 'log', 'YScale', 'log', 'Box', 'on', ...
        'FontSize', 7, 'FontName', 'Times New Roman');
hold(ax, 'on');
grid(ax, 'on');
ax.GridAlpha     = 0.25;
ax.GridLineStyle = ':';

% =========================================================================
%  4.  Plot 3 curves
% =========================================================================
h32 = plot(ax, crb(1).B_mhz, crb(1).r, '-o', ...
    'Color', clr32,  'LineWidth', lw, 'MarkerSize', msz, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr32, ...
    'MarkerIndices', mi);

h64 = plot(ax, crb(2).B_mhz, crb(2).r, '-s', ...
    'Color', clr64,  'LineWidth', lw, 'MarkerSize', msz, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr64, ...
    'MarkerIndices', mi);

h128 = plot(ax, crb(3).B_mhz, crb(3).r, '-^', ...
    'Color', clr128, 'LineWidth', lw, 'MarkerSize', msz, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr128, ...
    'MarkerIndices', mi);

% =========================================================================
%  5.  Gap annotation: M-doubling at B=400 MHz
%      Show one brace between M=64 and M=128 on the right side (B=800 MHz)
%      to avoid crowding the center.
% =========================================================================
% Annotate the M=64->M=128 gap at B=400 MHz with a vertical bracket
ann_B    = 400;
v64_ann  = crb(2).r(abs(crb(2).B_mhz  - ann_B) < 1e-6);
v128_ann = crb(3).r(abs(crb(3).B_mhz  - ann_B) < 1e-6);
gap_dB   = 20*log10(v64_ann / v128_ann);

tick_w = 0.04;   % half-width in log10(B) units
plot(ax, ann_B*[1 1], [v128_ann v64_ann], 'k-', 'LineWidth', 0.8);
plot(ax, ann_B*10.^[-tick_w tick_w], [v64_ann  v64_ann],  'k-', 'LineWidth', 0.8);
plot(ax, ann_B*10.^[-tick_w tick_w], [v128_ann v128_ann], 'k-', 'LineWidth', 0.8);
text(ax, ann_B * 10^(tick_w + 0.06), ...
    exp(0.5*(log(v64_ann)+log(v128_ann))), ...
    sprintf('%.1f dB', gap_dB), ...
    'FontSize', 5.5, 'FontName', 'Times New Roman', ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', 'k');

% =========================================================================
%  6.  Axes formatting
% =========================================================================
xlabel(ax, '$B$ [MHz]', 'FontSize', 7.5, 'FontName', 'Times New Roman', ...
    'Interpreter', 'latex');
ylabel(ax, '$\sqrt{\mathrm{CRB}_r}$ [m]', ...
    'FontSize', 7.5, 'FontName', 'Times New Roman', 'Interpreter', 'latex');

xlim(ax, [40 1000]);
ax.XTick = [50 100 200 400 800];
ax.XTickLabel = {'50','100','200','400','800'};

% y-limits: tight around data with ~15% margin each side in log space
all_r = [crb(1).r; crb(2).r; crb(3).r];
ylim(ax, [min(all_r)*0.82, max(all_r)*1.25]);

ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'off';
ax.YMinorTick = 'on';

% =========================================================================
%  7.  Legend (northeast, inside axes)
% =========================================================================
legend(ax, [h32, h64, h128], ...
    {'$M=32$', '$M=64$', '$M=128$'}, ...
    'Interpreter', 'latex', ...
    'FontSize', 6.5, 'FontName', 'Times New Roman', ...
    'Location', 'northeast', 'Box', 'on');

% =========================================================================
%  8.  Export PDF
% =========================================================================
pdf_path = fullfile(script_dir, 'fig7_crb_geodiv_bw.pdf');
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
fprintf('fig7: saved to %s\n', pdf_path);

% =========================================================================
%  LOCAL HELPER
% =========================================================================
function s = ok_str(flag)
if flag, s = 'OK'; else, s = '*** ANCHOR MISMATCH ***'; end
end

end  % function fig7_crb_geodiv_bw
