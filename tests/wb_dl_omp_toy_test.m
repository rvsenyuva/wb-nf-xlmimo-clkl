%WB_DL_OMP_TOY_TEST  Standalone acceptance test for wb_dl_omp.m.
%
%  Verifies four acceptance criteria (T1-T4) for the wideband DL-OMP
%  estimator (Paper C Baseline B5, Task 11.7).
%
%  Scene: M=64, N=64, K_s=16, B=400 MHz, Delta_f=25 MHz, d=1 path,
%         SNR=10 dB, r_hi_fac=0.30 (strong near-field, Lesson L33).
%
%  Acceptance tests:
%    T1 -- angle error < P.dtheta_tol  (5 deg, tight; full-array access)
%    T2 -- relative range error < P.dr_fac_tol  (50%)
%    T3 -- law-of-sines range: info.r_los > 0  AND  < 2*r_hi_fac*r_RD
%    T4 -- info struct: K_iter_used>=1, d_hat==P.d, delta_m>0
%
%  Writes: wb_dl_omp_toy_diag.csv (17 fields)
%
%  Author  : R. V. Senyuva (Maltepe University)
%  Date    : May 2026
%  Task    : Paper C Phase 2, Task 11.7

clear; clc;
fprintf('=== WB_DL_OMP TOY TEST ===\n\n');

% ---- Random seed (single call, Lesson L29 discipline) -----------------
rng(42, 'twister');

% ---- Scene parameters -------------------------------------------------
P.c        = 3e8;
P.c0       = P.c;
P.fc       = 28e9;
P.lambda   = P.c / P.fc;
P.lambda_c = P.lambda;
P.M        = 64;
P.d_ant    = P.lambda_c / 2;
P.N_RF     = 8;              % not used by DL-OMP; included for struct compat.
P.N        = 64;
P.d        = 1;

% Range and angle bounds (Phase 2 convention)
D_ap       = (P.M - 1) * P.d_ant;
P.r_RD     = 2 * D_ap^2 / P.lambda_c;
P.r_lo_fac = 0.05;
% TOY-TEST ONLY: r_hi_fac = 0.30 (strong near-field, Lesson L33).
% Law-of-sines most reliable when kappa*m_bar_max^2 > 0.1 rad.
% Production Monte Carlo uses r_hi_fac = 1.0.
P.r_hi_fac = 0.30;
P.u_margin = 2.0;
P.u_min    = 1 / (P.r_hi_fac * P.r_RD * P.u_margin);
P.u_max    = 1 / (P.r_lo_fac * P.r_RD / P.u_margin);
P.theta_lo = 20 * pi / 180;
P.theta_hi = 60 * pi / 180;

% Dictionary sizing (reduced from 256 for speed in toy test)
P.Q_theta  = 64;

% Wideband: B=400 MHz, K_s=16, Delta_f=25 MHz (Lesson L29)
P.K_s      = 16;
P.B        = 400e6;
P.Delta_f  = P.B / P.K_s;                              % 25 MHz
P.K        = P.K_s;
k_idx      = (-(P.K_s/2):(P.K_s/2-1)).';
P.alpha_k_vec = 1 + k_idx * P.Delta_f / P.fc;         % K_s x 1
P.k_indices   = k_idx;

SNR_dB     = 10;

% Tolerance parameters for acceptance tests
P.dtheta_tol  = 5 * pi / 180;   % T1: 5 deg angle tolerance (full-array)
P.dr_fac_tol  = 0.50;           % T2: 50% relative range tolerance

% ---- Data generation --------------------------------------------------
% Suppress channel-gen CSV output in toy test
P.wb_gen_write_csv = false;

fprintf('Generating wideband channel data ...\n');
[X_full, ~, ~, theta_gen, r_gen, ~, ~, ~] = wb_channel_gen_ofdm_nf(P, SNR_dB);

% Option (b): adopt drawn scene as ground truth (seed-independent pass/fail)
theta_true = theta_gen(1);   % [rad]
r_true     = r_gen(1);       % [m]

fprintf('Ground truth: theta = %.2f deg,  r = %.4f m\n', ...
    theta_true * 180/pi, r_true);
fprintf('r_RD = %.4f m,  r_true/r_RD = %.4f\n', P.r_RD, r_true/P.r_RD);
fprintf('X_full size: %d x %d x %d\n\n', size(X_full,1), size(X_full,2), size(X_full,3));

% ---- Run estimator ----------------------------------------------------
fprintf('Running wb_dl_omp ...\n');
t_start = tic;
[theta_hat, r_hat, info] = wb_dl_omp(X_full, P);
t_elapsed = toc(t_start);

theta_hat_deg = theta_hat(1) * 180/pi;
r_hat_m       = r_hat(1);
theta_err_deg = abs(theta_hat(1) - theta_true) * 180/pi;
r_err_m       = abs(r_hat_m - r_true);
r_err_frac    = r_err_m / r_true;

fprintf('Estimated: theta = %.2f deg,  r = %.4f m\n', theta_hat_deg, r_hat_m);
fprintf('Errors:    dtheta = %.3f deg,  dr = %.4f m (%.1f%%)\n', ...
    theta_err_deg, r_err_m, r_err_frac * 100);
fprintf('Elapsed time: %.2f s\n\n', t_elapsed);

% ---- Acceptance tests T1-T4 -------------------------------------------
fprintf('--- Acceptance Tests ---\n');

% T1: angle accuracy
T1_pass = (theta_err_deg < (P.dtheta_tol * 180/pi));
fprintf('T1 (angle < %.1f deg):  err = %.3f deg  -->  %s\n', ...
    P.dtheta_tol * 180/pi, theta_err_deg, pass_str(T1_pass));

% T2: relative range accuracy
T2_pass = (r_err_frac < P.dr_fac_tol);
fprintf('T2 (range < %.0f%%):       err = %.1f%%       -->  %s\n', ...
    P.dr_fac_tol * 100, r_err_frac * 100, pass_str(T2_pass));

% T3: law-of-sines range sanity
r_los_val = info.r_los;
r_los_lo  = 0;
r_los_hi  = 2 * P.r_hi_fac * P.r_RD;
T3_pass   = (r_los_val > r_los_lo) && (r_los_val < r_los_hi);
fprintf('T3 (r_los in (0, %.4f)):  r_los = %.4f  -->  %s\n', ...
    r_los_hi, r_los_val, pass_str(T3_pass));

% T4: info struct populated correctly
T4_pass = (info.K_iter_used >= 1) && ...
          (info.d_hat == P.d) && ...
          (info.delta_m > 0);
fprintf('T4 (info struct):  K_iter_used=%d, d_hat=%d, delta_m=%.4f  -->  %s\n', ...
    info.K_iter_used, info.d_hat, info.delta_m, pass_str(T4_pass));

all_pass = T1_pass && T2_pass && T3_pass && T4_pass;
fprintf('\n=== OVERALL: %s ===\n\n', pass_str(all_pass));

% ---- Write diagnostic CSV ---------------------------------------------
csv_path = 'wb_dl_omp_toy_diag.csv';

fid = fopen(csv_path, 'w');
if fid < 0
    error('wb_dl_omp_toy_test: could not open %s for writing', csv_path);
end

% Header
fprintf(fid, ['theta_true_deg,r_true_m,theta_hat_deg,r_hat_m,' ...
    'theta_err_deg,r_err_m,r_err_frac,' ...
    'SNR_dB,K_s,d,N_RF,M,Q_theta,K_iter,' ...
    'T1_pass,T2_pass,T3_pass,T4_pass\n']);

% Data row
fprintf(fid, ['%.6f,%.6f,%.6f,%.6f,' ...
    '%.6f,%.6f,%.6f,' ...
    '%d,%d,%d,%d,%d,%d,%d,' ...
    '%d,%d,%d,%d\n'], ...
    theta_true * 180/pi, r_true, theta_hat_deg, r_hat_m, ...
    theta_err_deg, r_err_m, r_err_frac, ...
    SNR_dB, P.K_s, P.d, P.N_RF, P.M, P.Q_theta, info.K_iter_used, ...
    T1_pass, T2_pass, T3_pass, T4_pass);

fclose(fid);
fprintf('Diagnostic CSV written to: %s\n', csv_path);

% ---- Additional diagnostics -------------------------------------------
fprintf('\n--- Additional Info ---\n');
fprintf('  alpha_k range:  [%.6f, %.6f]\n', ...
    min(P.alpha_k_vec), max(P.alpha_k_vec));
fprintf('  Subarray sep. delta = %.4f m  (%.1f lambda_c)\n', ...
    info.delta_m, info.delta_m / P.lambda_c);
fprintf('  r_los (raw law-of-sines) = %.4f m\n', info.r_los);
fprintf('  r_hat (clamped estimate) = %.4f m\n', r_hat_m);
fprintf('  r_min = %.4f m,  r_max = %.4f m\n', ...
    P.r_lo_fac * P.r_RD, P.r_hi_fac * P.r_RD);


% ---- Helper -----------------------------------------------------------
function s = pass_str(flag)
    if flag
        s = 'PASS';
    else
        s = 'FAIL';
    end
end
