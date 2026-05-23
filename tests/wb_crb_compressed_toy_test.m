% =========================================================================
% wb_crb_compressed_toy_test.m
% =========================================================================
% Task 11.5 toy-test and acceptance verification for wb_crb_compressed.m.
%
% Acceptance tests:
%   T1 -- FIM positive definiteness (min eigenvalue > 0)
%   T2 -- CRB decreases with SNR
%   T3 -- Data diversity scaling vs. narrowband (>= 10 dB at K_s=16)
%   T4 -- Random-W envelope sanity (nonzero std, std < 50% of mean)
%
% Output: wb_crb_compressed_toy_diag.csv
%
% Scene (per starter prompt Sec.8.1):
%   M=64, N_RF=8, N=64, fc=28 GHz, d=2, B=400 MHz, K_s=16, Delta_f=25 MHz
%   theta_true = [35; 45] deg, r_true = [5; 8] m, SNR=10 dB
%
% Author:  R. V. Senyuva (Maltepe University)
% Date:    May 2026
% =========================================================================

clear; close all; clc;
rng(42, 'twister');   % single seed call before data generation (Sec.8.1)

fprintf('=============================================================\n');
fprintf('  WB_CRB_COMPRESSED TOY TEST  (Task 11.5)\n');
fprintf('=============================================================\n\n');

%% ========================================================================
%  SECTION 1: PARAMETER SETUP (Sec.8.1 of starter prompt)
% =========================================================================

P.c       = 3e8;
P.c0      = P.c;
P.fc      = 28e9;
P.lambda  = P.c / P.fc;
P.lambda_c = P.lambda;
P.M       = 64;
P.d_ant   = P.lambda_c / 2;
P.N_RF    = 8;
P.N       = 64;          % snapshots
P.d       = 2;           % d=2 paths for Task 11.5 toy test

% Range bounds (Phase 2 convention)
D_ap    = (P.M - 1) * P.d_ant;
P.r_RD  = 2 * D_ap^2 / P.lambda_c;
P.u_min = 1 / (P.r_RD * 2.0);
P.u_max = 2.0 / (0.05 * P.r_RD);

% Wideband: B=400 MHz, K_s=16, Delta_f=25 MHz (Lesson L29)
P.K_s      = 16;
P.B        = 400e6;
P.Delta_f  = P.B / P.K_s;                        % 25 MHz (NOT 120 kHz)
k_idx      = (-(P.K_s/2) : (P.K_s/2 - 1)).';
P.alpha_k_vec = 1 + k_idx * P.Delta_f / P.fc;    % K_s x 1

% N_seed for random-W envelope (C1 deliverable)
P.N_seed = 50;

% True path parameters (d=2)
SNR_dB    = 10;
theta_true = [35; 45] * pi/180;   % [rad]
r_true     = [5; 8];              % [m]  -- well inside near-field
p_true     = [1.0; 0.8];
N0_true    = sum(p_true) / (10^(SNR_dB/10));

fprintf('  Array:   M=%d, N_RF=%d, fc=%.0f GHz\n', P.M, P.N_RF, P.fc/1e9);
fprintf('  Scene:   d=%d paths, theta=[%.0f, %.0f] deg, r=[%.0f, %.0f] m\n', ...
    P.d, theta_true(1)*180/pi, theta_true(2)*180/pi, r_true(1), r_true(2));
fprintf('  Wideband: B=%.0f MHz, K_s=%d, Delta_f=%.0f MHz\n', ...
    P.B/1e6, P.K_s, P.Delta_f/1e6);
fprintf('  SNR=%d dB, N0=%.4f, N=%d snapshots, N_seed=%d\n', ...
    SNR_dB, N0_true, P.N, P.N_seed);
fprintf('  r_RD = %.2f m\n', P.r_RD);
fprintf('  r_true(1)/r_RD = %.3f,  r_true(2)/r_RD = %.3f\n', ...
    r_true(1)/P.r_RD, r_true(2)/P.r_RD);
fprintf('-------------------------------------------------------------\n\n');

%% ========================================================================
%  T1: FIM POSITIVE DEFINITENESS
% =========================================================================
fprintf('=============================================================\n');
fprintf('  T1: FIM POSITIVE DEFINITENESS\n');
fprintf('=============================================================\n\n');

opts_J.return_J = true;
opts_J.verbose  = false;

[crb_r_nom, crb_theta_nom, info_nom] = wb_crb_compressed( ...
    theta_true, r_true, p_true, N0_true, P, opts_J);

J_test    = info_nom.J_WB;
eig_J     = real(eig(J_test));
min_eig   = min(eig_J);
T1_pass   = (min_eig > 0);
cond_J    = info_nom.cond_J_all(end);
n_sing    = max(info_nom.n_singular_all);

fprintf('  min(eig(J_WB)) = %.4e\n', min_eig);
fprintf('  cond(J_WB)     = %.4e\n', cond_J);
fprintf('  n_singular     = %d\n',   n_sing);
fprintf('  T1 PASS = %d  (min eigenvalue > 0: %d > 0)\n\n', T1_pass, T1_pass);

fprintf('  sqrt(CRB_r):     path 1 = %.4f m,  path 2 = %.4f m\n', ...
    crb_r_nom(1), crb_r_nom(2));
fprintf('  sqrt(CRB_theta): path 1 = %.4f deg, path 2 = %.4f deg\n\n', ...
    crb_theta_nom(1), crb_theta_nom(2));

%% ========================================================================
%  T2: CRB DECREASES WITH SNR
% =========================================================================
fprintf('=============================================================\n');
fprintf('  T2: CRB DECREASES WITH SNR\n');
fprintf('=============================================================\n\n');

SNR_vals = [0, 10, 20];   % dB
crb_r_snr     = zeros(P.d, numel(SNR_vals));
crb_theta_snr = zeros(P.d, numel(SNR_vals));

opts_t2.verbose = false;
for is = 1:numel(SNR_vals)
    N0_s = sum(p_true) / (10^(SNR_vals(is)/10));
    [crb_r_snr(:,is), crb_theta_snr(:,is)] = wb_crb_compressed( ...
        theta_true, r_true, p_true, N0_s, P, opts_t2);
    fprintf('  SNR=%3d dB | sqrt_crb_r = [%.4f, %.4f] m | sqrt_crb_theta = [%.4f, %.4f] deg\n', ...
        SNR_vals(is), crb_r_snr(1,is), crb_r_snr(2,is), ...
        crb_theta_snr(1,is), crb_theta_snr(2,is));
end

% Verify monotone decrease for both paths
T2_pass = all(crb_r_snr(:,1) > crb_r_snr(:,2)) && ...
          all(crb_r_snr(:,2) > crb_r_snr(:,3));

fprintf('\n  T2 PASS = %d', T2_pass);
if T2_pass
    fprintf('  (CRB_r monotone decreasing with SNR for both paths)\n\n');
else
    fprintf('  *** FAIL: CRB_r not monotone in SNR ***\n\n');
end

%% ========================================================================
%  T3: DATA DIVERSITY SCALING (K_s=16 vs K_s=1)
% =========================================================================
fprintf('=============================================================\n');
fprintf('  T3: DATA DIVERSITY SCALING (K_s=16 vs K_s=1)\n');
fprintf('=============================================================\n\n');

% Narrowband: single subcarrier at carrier frequency (K_s=1, alpha_k=1)
P_nb           = P;
P_nb.K_s       = 1;
P_nb.alpha_k_vec = 1.0;   % alpha_k = 1 at carrier frequency

opts_t3.verbose = false;
[crb_r_nb, ~] = wb_crb_compressed(theta_true, r_true, p_true, N0_true, P_nb, opts_t3);
[crb_r_wb, ~] = wb_crb_compressed(theta_true, r_true, p_true, N0_true, P,    opts_t3);

% Data diversity gain in dB (K_s=16 => theoretical ~10*log10(16) = 12.04 dB)
gain_dB_path1 = 20*log10(crb_r_nb(1) / crb_r_wb(1));
gain_dB_path2 = 20*log10(crb_r_nb(2) / crb_r_wb(2));
threshold_dB  = 10*log10(P.K_s) - 2;   % >= 10 dB for K_s=16 (12.04 - 2 = 10.04 dB)

fprintf('  K_s=1  sqrt_crb_r = [%.4f, %.4f] m\n', crb_r_nb(1), crb_r_nb(2));
fprintf('  K_s=16 sqrt_crb_r = [%.4f, %.4f] m\n', crb_r_wb(1), crb_r_wb(2));
fprintf('  Diversity gain:  path 1 = %.2f dB,  path 2 = %.2f dB\n', ...
    gain_dB_path1, gain_dB_path2);
fprintf('  Threshold: >= %.2f dB (= 10*log10(%d) - 2)\n', threshold_dB, P.K_s);

T3_pass = (gain_dB_path1 >= threshold_dB) && (gain_dB_path2 >= threshold_dB);

fprintf('  T3 PASS = %d', T3_pass);
if T3_pass
    fprintf('  (both paths show >= %.2f dB diversity gain)\n\n', threshold_dB);
else
    fprintf('  *** FAIL: diversity gain below threshold ***\n\n');
end

%% ========================================================================
%  T4: RANDOM-W ENVELOPE SANITY (N_seed = 50)
% =========================================================================
fprintf('=============================================================\n');
fprintf('  T4: RANDOM-W ENVELOPE SANITY\n');
fprintf('=============================================================\n\n');

% Use the nominal result (already computed with N_seed=50)
crb_r_std1 = info_nom.crb_r_std(1);
crb_r_std2 = info_nom.crb_r_std(2);
crb_r_mean1 = info_nom.crb_r_mean(1);
crb_r_mean2 = info_nom.crb_r_mean(2);
cv1 = crb_r_std1 / crb_r_mean1;   % coefficient of variation, path 1
cv2 = crb_r_std2 / crb_r_mean2;   % coefficient of variation, path 2

fprintf('  Path 1: sqrt_crb_r mean=%.4f m, std=%.4f m, CV=%.3f\n', ...
    crb_r_mean1, crb_r_std1, cv1);
fprintf('  Path 2: sqrt_crb_r mean=%.4f m, std=%.4f m, CV=%.3f\n', ...
    crb_r_mean2, crb_r_std2, cv2);
fprintf('  10th percentile: [%.4f, %.4f] m\n', ...
    info_nom.crb_r_p10(1), info_nom.crb_r_p10(2));
fprintf('  90th percentile: [%.4f, %.4f] m\n', ...
    info_nom.crb_r_p90(1), info_nom.crb_r_p90(2));

% T4 checks: nonzero variation AND CV < 0.5
T4_pass = (crb_r_std1 > 0) && (cv1 < 0.5) && ...
          (crb_r_std2 > 0) && (cv2 < 0.5);

fprintf('  T4 PASS = %d', T4_pass);
if T4_pass
    fprintf('  (nonzero std, CV < 0.5 for both paths)\n\n');
else
    fprintf('  *** FAIL: std=0 or CV >= 0.5 ***\n\n');
    if crb_r_std1 == 0 || crb_r_std2 == 0
        fprintf('  --> std is zero: check RNG seeds and N_seed\n');
    end
    if cv1 >= 0.5 || cv2 >= 0.5
        fprintf('  --> CV too large: combiner variation dominates the bound\n');
    end
end

%% ========================================================================
%  SUMMARY TABLE
% =========================================================================
fprintf('=============================================================\n');
fprintf('  ACCEPTANCE SUMMARY\n');
fprintf('=============================================================\n');
fprintf('  T1 (FIM positive definite):          %s\n', pass_str(T1_pass));
fprintf('  T2 (CRB decreases with SNR):         %s\n', pass_str(T2_pass));
fprintf('  T3 (data diversity >= %.1f dB):      %s\n', threshold_dB, pass_str(T3_pass));
fprintf('  T4 (random-W envelope sanity):       %s\n', pass_str(T4_pass));
fprintf('=============================================================\n\n');

all_pass = T1_pass && T2_pass && T3_pass && T4_pass;
if all_pass
    fprintf('  ALL TESTS PASSED -- Task 11.5 acceptance criteria met.\n\n');
else
    fprintf('  *** ONE OR MORE TESTS FAILED -- inspect diagnostics above. ***\n\n');
end

%% ========================================================================
%  CSV OUTPUT: wb_crb_compressed_toy_diag.csv
% =========================================================================

T_out = table( ...
    theta_true(1)*180/pi, ...
    r_true(1), ...
    theta_true(2)*180/pi, ...
    r_true(2), ...
    SNR_dB, ...
    P.K_s, ...
    P.d, ...
    P.N_RF, ...
    P.M, ...
    crb_theta_nom(1), ...
    crb_r_nom(1), ...
    crb_theta_nom(2), ...
    crb_r_nom(2), ...
    crb_r_std1, ...
    crb_r_std2, ...
    T1_pass, ...
    T2_pass, ...
    T3_pass, ...
    T4_pass, ...
    cond_J, ...
    n_sing, ...
    'VariableNames', { ...
        'theta_true_1_deg', 'r_true_1_m', ...
        'theta_true_2_deg', 'r_true_2_m', ...
        'SNR_dB', 'K_s', 'd', 'N_RF', 'M', ...
        'sqrt_crb_theta_1_deg', 'sqrt_crb_r_1_m', ...
        'sqrt_crb_theta_2_deg', 'sqrt_crb_r_2_m', ...
        'crb_r_std_1_m', 'crb_r_std_2_m', ...
        'T1_pass', 'T2_pass', 'T3_pass', 'T4_pass', ...
        'cond_J', 'n_singular'});

writetable(T_out, 'wb_crb_compressed_toy_diag.csv');
fprintf('  Exported: wb_crb_compressed_toy_diag.csv\n\n');

%% ========================================================================
%  LOCAL HELPER
% =========================================================================
function s = pass_str(flag)
if flag
    s = 'PASS';
else
    s = 'FAIL';
end
end
