% =========================================================================
% fig6_crb_d2_snr.m
% =========================================================================
% Phase 3.5 Subtask 3.5.3 -- Plotting script for Fig 6.
%
% Fig 6: Compressed-domain sqrt-CRB on path-1 range vs SNR
%        for N_RF in {8, 16, 32} plus full-array N_RF=64 reference.
%        d=2 paths, B=400 MHz, M=64.
%
% Key annotations:
%   8.1 dB gap  : N_RF=8  vs N_RF=32  @ SNR=10 dB  (compression cost)
%   11.4 dB gap : N_RF=8  vs N_RF=64  @ SNR=10 dB  (total compression cost)
%
% Source CSV: wb_crb_d2_sweep_20260526_225015_v3.csv  (52 rows x 10 cols)
% Output PDF: fig6_crb_d2_snr.pdf  (same folder as this script)
%
% Color scheme (NRF-specific, no B1/B2/B4/B5 estimator palette):
%   N_RF=8  : crimson    [0.6350 0.0780 0.1840]  solid, circle
%   N_RF=16 : navy       [0.0000 0.3176 0.6275]  solid, square
%   N_RF=32 : forest grn [0.1176 0.4706 0.1176]  solid, triangle-up
%   N_RF=64 : dark gray  [0.3000 0.3000 0.3000]  dashed, diamond (full-array ref)
%
% Markers every 3rd point (indices 1,4,7,10,13 of 13 pts).
% No title. Legend inside axes (southwest).
%
% Script location: <project_root>/src/plotting/fig6_crb_d2_snr.m
% CSV location:    <project_root>/src/results/paperC_phase35_fig6_crb_d2_v2/
%
% Locked anchors (CSV-verified 2026-05-26):
%   crb_r1 @ SNR=10 dB, N_RF=8  = 0.04564 m
%   crb_r1 @ SNR=10 dB, N_RF=32 = 0.01793 m
%   crb_r1 @ SNR=10 dB, N_RF=64 = 0.01228 m
%   Gap N_RF=8 vs N_RF=32 = 8.12 dB
%   Gap N_RF=8 vs N_RF=64 = 11.40 dB
%
% Author : R. V. Senyuva (Maltepe University)
% Date   : May 2026  (Phase 3.5)
% =========================================================================

clear; close all; clc;

% =========================================================================
%  0.  Paths  (script at src/plotting/; CSV and PDF at src/results/...)
% =========================================================================
script_dir = fileparts(mfilename('fullpath'));      % .../src/plotting
src_dir    = fullfile(script_dir, '..');            % .../src
csv_dir    = fullfile(src_dir, 'results', 'paperC_phase35_fig6_crb_d2_v2');
csv_name   = 'wb_crb_d2_sweep_20260526_225015_v3.csv';
csv_path   = fullfile(csv_dir, csv_name);

if ~isfile(csv_path)
    listing = dir(fullfile(csv_dir, 'wb_crb_d2_sweep_*_v3.csv'));
    if isempty(listing)
        error('fig6_crb_d2_snr: CSV not found in %s', csv_dir);
    end
    csv_path = fullfile(csv_dir, listing(end).name);
    fprintf('  Using CSV: %s\n', csv_path);
end

% =========================================================================
%  1.  Load and sort data
% =========================================================================
T = readtable(csv_path);

nrf_vec = [8, 16, 32, 64];
n_snr   = numel(unique(T.snr_db));   % 13

crb = struct();
for ii = 1:numel(nrf_vec)
    mask = T.N_RF == nrf_vec(ii);
    sub  = T(mask, :);
    [~, idx]    = sort(sub.snr_db);
    crb(ii).snr = sub.snr_db(idx);
    crb(ii).r1  = sub.crb_r1_m(idx);
end

% Anchor verification
v8  = crb(1).r1(abs(crb(1).snr - 10) < 1e-9);
v32 = crb(3).r1(abs(crb(3).snr - 10) < 1e-9);
v64 = crb(4).r1(abs(crb(4).snr - 10) < 1e-9);
gap_8_32 = 20*log10(v8/v32);
gap_8_64 = 20*log10(v8/v64);
fprintf('  Anchor N_RF=8  @ SNR=10: %.5f m  %s\n', v8,  ok_str(abs(v8 -0.04564)<=1e-4));
fprintf('  Anchor N_RF=32 @ SNR=10: %.5f m  %s\n', v32, ok_str(abs(v32-0.01793)<=1e-4));
fprintf('  Anchor N_RF=64 @ SNR=10: %.5f m  %s\n', v64, ok_str(abs(v64-0.01228)<=1e-4));
fprintf('  Gap 8 vs 32: %.4f dB  %s\n', gap_8_32, ok_str(abs(gap_8_32-8.12)<=0.05));
fprintf('  Gap 8 vs 64: %.4f dB  %s\n\n', gap_8_64, ok_str(abs(gap_8_64-11.40)<=0.05));

% =========================================================================
%  2.  Style
% =========================================================================
clr8  = [0.6350 0.0780 0.1840];   % crimson
clr16 = [0.0000 0.3176 0.6275];   % navy
clr32 = [0.1176 0.4706 0.1176];   % forest green
clr64 = [0.3000 0.3000 0.3000];   % dark gray (full-array reference)

lw_comp = 1.5;    % compressed curves
lw_full = 1.2;    % full-array reference (slightly thinner, dashed)
msz     = 7;
mk_idx  = 1 : 3 : n_snr;    % markers at indices 1,4,7,10,13

% =========================================================================
%  3.  Figure setup (IEEEtran single-column: 8.8 cm x 6.5 cm)
% =========================================================================
fig = figure('Units', 'centimeters', 'Position', [2 2 8.8 6.5]);
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [8.8 6.5], ...
         'PaperPosition', [0 0 8.8 6.5]);

ax = axes(fig);
set(ax, 'YScale', 'log', 'Box', 'on', 'FontSize', 7, ...
        'FontName', 'Times New Roman');
hold(ax, 'on');
grid(ax, 'on');

% =========================================================================
%  4.  Plot 4 curves
% =========================================================================
% Compressed curves (solid)
h8 = plot(ax, crb(1).snr, crb(1).r1, '-', ...
    'Color', clr8,  'LineWidth', lw_comp, ...
    'Marker', 'o',  'MarkerSize', msz, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr8, ...
    'MarkerIndices', mk_idx);

h16 = plot(ax, crb(2).snr, crb(2).r1, '-', ...
    'Color', clr16, 'LineWidth', lw_comp, ...
    'Marker', 's',  'MarkerSize', msz, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr16, ...
    'MarkerIndices', mk_idx);

h32 = plot(ax, crb(3).snr, crb(3).r1, '-', ...
    'Color', clr32, 'LineWidth', lw_comp, ...
    'Marker', '^',  'MarkerSize', msz, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr32, ...
    'MarkerIndices', mk_idx);

% Full-array reference (dashed, dark gray, diamond markers)
h64 = plot(ax, crb(4).snr, crb(4).r1, '--', ...
    'Color', clr64, 'LineWidth', lw_full, ...
    'Marker', 'd',  'MarkerSize', msz, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr64, ...
    'MarkerIndices', mk_idx);

% =========================================================================
%  5.  Gap annotations at SNR=10 dB
% =========================================================================
snr_ann  = 10;
tick_len = 0.25;   % half-width of end ticks [dB]

% --- Gap 1: N_RF=8 vs N_RF=32 (8.1 dB) ---
plot(ax, [snr_ann snr_ann],              [v32 v8],   'k-', 'LineWidth', 0.8);
plot(ax, snr_ann + [-tick_len tick_len], [v8  v8],   'k-', 'LineWidth', 0.8);
plot(ax, snr_ann + [-tick_len tick_len], [v32 v32],  'k-', 'LineWidth', 0.8);
text(ax, snr_ann + 0.6, exp(0.5*(log(v8)+log(v32))), '8.1 dB', ...
    'FontSize', 6, 'FontName', 'Times New Roman', ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', 'k');

% --- Gap 2: N_RF=32 vs N_RF=64 (3.3 dB) annotated to the right ---
snr_ann2 = 17.5;
v32_s2 = crb(3).r1(abs(crb(3).snr - snr_ann2) < 1e-9);
v64_s2 = crb(4).r1(abs(crb(4).snr - snr_ann2) < 1e-9);
plot(ax, [snr_ann2 snr_ann2],              [v64_s2 v32_s2], 'k-', 'LineWidth', 0.8);
plot(ax, snr_ann2 + [-tick_len tick_len],  [v32_s2 v32_s2], 'k-', 'LineWidth', 0.8);
plot(ax, snr_ann2 + [-tick_len tick_len],  [v64_s2 v64_s2], 'k-', 'LineWidth', 0.8);
text(ax, snr_ann2 + 0.6, exp(0.5*(log(v32_s2)+log(v64_s2))), '3.3 dB', ...
    'FontSize', 6, 'FontName', 'Times New Roman', ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', 'k');

% =========================================================================
%  6.  Axes formatting
% =========================================================================
xlabel(ax, 'SNR [dB]', 'FontSize', 7.5, 'FontName', 'Times New Roman');
ylabel(ax, '$\sqrt{\mathrm{CRB}_r}$ [m]', ...
    'FontSize', 7.5, 'FontName', 'Times New Roman', 'Interpreter', 'latex');

xlim(ax, [-5 25]);
ylim(ax, [min(crb(4).r1)*0.88,  max(crb(1).r1)*1.15]);

% =========================================================================
%  7.  Legend
% =========================================================================
legend(ax, [h8, h16, h32, h64], ...
    {'$N_{\rm RF}=8$', '$N_{\rm RF}=16$', '$N_{\rm RF}=32$', ...
     '$N_{\rm RF}=64$ (full array)'}, ...
    'Interpreter', 'latex', ...
    'FontSize', 6.5, 'FontName', 'Times New Roman', ...
    'Location', 'northeast', 'Box', 'on');

% =========================================================================
%  8.  Export PDF  (same folder as this script)
% =========================================================================
pdf_path = fullfile(script_dir, 'fig6_crb_d2_snr.pdf');
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
fprintf('  Saved: %s\n', pdf_path);

% =========================================================================
%  LOCAL HELPER
% =========================================================================
function s = ok_str(flag)
if flag, s = 'OK'; else, s = '*** ANCHOR MISMATCH ***'; end
end
