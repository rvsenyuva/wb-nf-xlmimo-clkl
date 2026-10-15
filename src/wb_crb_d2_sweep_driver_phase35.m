% =========================================================================
% wb_crb_d2_sweep_driver_phase35.m
% =========================================================================
% Phase 3.5 Subtask 3.5.3 -- Fig 6 CRB rerun driver (v2: adds N_RF=64).
%
% Produces wb_crb_d2_sweep_<date>_v3.csv  (52 rows x 10 cols)
%   = 13 SNR pts x 4 N_RF values {8, 16, 32, 64}
%
% TC-4 canonical scene (confirmed from driver_script_CRB_sweep.m):
%   theta_true = [30; 50] deg
%   r_true     = [3.0; 7.0] m
%   p_true     = [1.0; 1.0]  (equal power)
%   N0         = p_true(1) / 10^(SNR/10)  [single-path SNR reference]
%
% Changes from v2 driver:
%   - nrf_vec = [8, 16, 32, 64]  (added N_RF=64 full-array reference)
%   - output filename suffix: _v3.csv
%   - expected dims: 52 rows x 10 cols
%
% Anchors to preserve (within +/- 1e-4 m):
%   crb_r1_m @ SNR=10 dB, N_RF=8  = 0.04564 m
%   crb_r1_m @ SNR=10 dB, N_RF=32 = 0.01793 m
%   Compression gap (8 vs 32) = 8.12 +/- 0.05 dB
%
% Author : R. V. Senyuva (Maltepe University)
% Date   : May 2026  (Phase 3.5 Subtask 3.5.3, v2 driver)
% =========================================================================

clear; close all; clc;

fprintf('============================================================\n');
fprintf('  WB_CRB_D2_SWEEP -- Phase 3.5 Subtask 3.5.3 (v3 CSV)\n');
fprintf('============================================================\n\n');

% =========================================================================
%  1.  Base parameter struct
% =========================================================================
P = setup_production_P_v4('snr', 'full');

P.d     = 2;
P.d_max = 2;

fprintf('  d      = %d paths\n', P.d);
fprintf('  N_seed = %d  (random-W draws for CRB envelope)\n', P.N_seed);

% =========================================================================
%  2.  Sweep parameters
% =========================================================================
snr_vec  = -5 : 2.5 : 25;        % 13 pts  (unchanged)
nrf_vec  = [8, 16, 32, 64];       % added N_RF=64 full-array reference
n_snr    = numel(snr_vec);
n_nrf    = numel(nrf_vec);

fprintf('  snr_vec  = [%s] dB  (%d pts)\n', num2str(snr_vec,'%.1f '), n_snr);
fprintf('  nrf_vec  = [%d %d %d %d]\n', nrf_vec(1), nrf_vec(2), nrf_vec(3), nrf_vec(4));
fprintf('  Total sweep points: %d x %d = %d\n\n', n_snr, n_nrf, n_snr*n_nrf);

% =========================================================================
%  3.  Scene (TC-4 canonical -- unchanged)
% =========================================================================
theta_true = [30; 50] * pi / 180;
r_true     = [3.0; 7.0];
p_true     = [1.0; 1.0];

fprintf('  Scene: theta=[%.0f, %.0f] deg, r=[%.1f, %.1f] m\n', ...
    theta_true(1)*180/pi, theta_true(2)*180/pi, r_true(1), r_true(2));
fprintf('  p_true=[%.1f, %.1f],  N0=p(1)/10^(SNR/10)\n\n', p_true(1), p_true(2));

% =========================================================================
%  4.  Allocate output arrays
% =========================================================================
n_rows = n_snr * n_nrf;
out_snr_db     = zeros(n_rows, 1);
out_N_RF       = zeros(n_rows, 1);
out_crb_r1     = zeros(n_rows, 1);
out_crb_r2     = zeros(n_rows, 1);
out_crb_theta1 = zeros(n_rows, 1);
out_crb_theta2 = zeros(n_rows, 1);
out_crb_r1_std = zeros(n_rows, 1);
out_crb_r2_std = zeros(n_rows, 1);
out_cond_J     = zeros(n_rows, 1);
out_n_sing     = zeros(n_rows, 1);

% =========================================================================
%  5.  Main sweep loop
% =========================================================================
row_idx = 0;

for i_nrf = 1 : n_nrf

    P_local      = P;
    P_local.N_RF = nrf_vec(i_nrf);
    fprintf('------ N_RF = %d ------\n', P_local.N_RF);

    for i_snr = 1 : n_snr

        snr_dB = snr_vec(i_snr);
        N0     = p_true(1) / (10^(snr_dB/10));

        [crb_r, crb_theta, info] = wb_crb_compressed( ...
            theta_true, r_true, p_true, N0, P_local);

        row_idx = row_idx + 1;
        out_snr_db(row_idx)     = snr_dB;
        out_N_RF(row_idx)       = P_local.N_RF;
        out_crb_r1(row_idx)     = crb_r(1);
        out_crb_r2(row_idx)     = crb_r(2);
        out_crb_theta1(row_idx) = crb_theta(1);
        out_crb_theta2(row_idx) = crb_theta(2);
        out_crb_r1_std(row_idx) = info.crb_r_std(1);
        out_crb_r2_std(row_idx) = info.crb_r_std(2);
        out_cond_J(row_idx)     = mean(info.cond_J_all);
        out_n_sing(row_idx)     = mean(info.n_singular_all);

        fprintf('  SNR=%5.1f dB | crb_r1=%.5f m | crb_r2=%.5f m | cond=%.2e\n', ...
            snr_dB, crb_r(1), crb_r(2), out_cond_J(row_idx));
    end
    fprintf('\n');
end

% =========================================================================
%  6.  Anchor check
% =========================================================================
fprintf('============================================================\n');
fprintf('  ANCHOR VERIFICATION\n');
fprintf('============================================================\n');

tol_m   = 1e-4;
tol_gap = 0.05;

snr10_nrf8  = (abs(out_snr_db - 10) < 1e-9) & (out_N_RF == 8);
snr10_nrf32 = (abs(out_snr_db - 10) < 1e-9) & (out_N_RF == 32);
snr10_nrf64 = (abs(out_snr_db - 10) < 1e-9) & (out_N_RF == 64);

v8  = out_crb_r1(snr10_nrf8);
v32 = out_crb_r1(snr10_nrf32);
v64 = out_crb_r1(snr10_nrf64);

gap_8_32  = 20*log10(v8  / v32);
gap_8_64  = 20*log10(v8  / v64);
gap_32_64 = 20*log10(v32 / v64);

fprintf('  crb_r1 @ SNR=10, N_RF=8  : %.5f m  (expected 0.04564)  %s\n', v8,  pass_str(abs(v8 -0.04564)<=tol_m));
fprintf('  crb_r1 @ SNR=10, N_RF=32 : %.5f m  (expected 0.01793)  %s\n', v32, pass_str(abs(v32-0.01793)<=tol_m));
fprintf('  crb_r1 @ SNR=10, N_RF=64 : %.5f m  (new, informational)\n',   v64);
fprintf('  Gap N_RF=8  vs N_RF=32   : %.4f dB  (expected 8.12)    %s\n', gap_8_32,  pass_str(abs(gap_8_32 -8.12)<=tol_gap));
fprintf('  Gap N_RF=8  vs N_RF=64   : %.4f dB  (informational)\n',       gap_8_64);
fprintf('  Gap N_RF=32 vs N_RF=64   : %.4f dB  (informational)\n',       gap_32_64);
fprintf('  Dimensions               : %d rows x 10 cols  %s\n', row_idx, pass_str(row_idx == 52));
fprintf('============================================================\n\n');

% =========================================================================
%  7.  CSV export
% =========================================================================
ts      = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = fullfile('src', 'results', 'paperC_phase35_fig6_crb_d2_v2');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
csv_name = fullfile(out_dir, ['wb_crb_d2_sweep_' ts '_v3.csv']);

T_out = table( ...
    out_snr_db(1:row_idx), ...
    out_N_RF(1:row_idx), ...
    out_crb_r1(1:row_idx), ...
    out_crb_r2(1:row_idx), ...
    out_crb_theta1(1:row_idx), ...
    out_crb_theta2(1:row_idx), ...
    out_crb_r1_std(1:row_idx), ...
    out_crb_r2_std(1:row_idx), ...
    out_cond_J(1:row_idx), ...
    out_n_sing(1:row_idx), ...
    'VariableNames', { ...
        'snr_db', 'N_RF', ...
        'crb_r1_m', 'crb_r2_m', ...
        'crb_theta1_deg', 'crb_theta2_deg', ...
        'crb_r1_std_m', 'crb_r2_std_m', ...
        'cond_J_mean', 'n_singular_mean'});

writetable(T_out, csv_name);
fprintf('  Exported: %s\n', csv_name);
fprintf('  Rows: %d  Cols: %d\n\n', height(T_out), width(T_out));

% =========================================================================
%  LOCAL HELPER
% =========================================================================
function s = pass_str(flag)
if flag, s = 'PASS'; else, s = 'FAIL'; end
end
