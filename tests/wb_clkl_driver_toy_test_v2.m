%WB_CLKL_DRIVER_TOY_TEST_V2  Task 11.3 (rev3) C2 verification for
%                            wb_clkl_driver.m (v2, multi-start).
%
%  Paper C Phase 3, Task 11.3 (rev3) Component 2.
%  Verifies the v2 4-phase WB-CL-KL driver against six acceptance checks:
%    T1 -- Convergence check          (info.converged == true)
%    T2 -- Monotone decrease          (Armijo guarantee)
%    T3 -- Improvement over BPD warm-start
%    T4 -- N0 estimate sanity         (within factor-of-3 of N0_true)
%    T5 -- Multi-start activation     (info.best_start in {1,...,n_ms_starts})
%    T6 -- Multi-start non-regression (L_ms_final <= L_ss_final + tol)
%
%  T1-T4 run with P.use_multi_start = false (single-start regression);
%  the existing v1 behaviour must be preserved bit-for-bit.
%  T5-T6 run with P.use_multi_start = true and use the TRUE (theta, r) as
%  warm-start, per starter prompt Sec.6.2 Changes 2-3.
%
%  SCENE PARAMETERS (from Task 11.4 starter prompt Sec. 8.1, unchanged):
%    M=64, N=64, N_RF=8, K_s=16, d=1, SNR=10 dB
%    rng(42, 'twister')
%    use_exact = true  (USW truth model, Branch B)
%
%  USAGE
%    Run from Paper C repo root:  >> wb_clkl_driver_toy_test_v2
%    All output is printed to the console.
%    Writes wb_clkl_driver_toy_test_v2_diag.csv to current directory.
%
%  DEPENDENCIES
%    wb_channel_gen_ofdm_nf.m  (Task 11.2)
%    bpd_baseline.m            (Task 11.1 / Phase 2)
%    wb_clkl_driver.m          (Task 11.3 rev3 C2 -- under test, from
%                               wb_clkl_driver_v2.m placed in repo root)
%    wb_clkl_estimator.m       (Task 11.3 rev2)
%    wb_nf_fresnel_steer.m     (Task 11.1)
%    nf_usw_steer.m            (Paper B)
%
%  Author : R. V. Senyuva (Maltepe University)
%  Date   : May 2026 (v2, Task 11.3 rev3 Component 2)

clear; clc;

fprintf('=============================================================\n');
fprintf(' wb_clkl_driver_toy_test_v2 -- Task 11.3 (rev3) C2 verification\n');
fprintf('=============================================================\n');

% =========================================================================
%  SECTION 1: PARAMETER STRUCT
%  Build a minimal Paper C parameter struct consistent with Task 11.2/11.3
%  toy tests.  Lessons applied: set BOTH P.c0=P.c AND P.lambda_c=P.lambda.
%
%  v2 NOTE: P.use_multi_start defaults to true in the v2 driver.  T1-T4
%  test the single-start regression and therefore explicitly set
%  P.use_multi_start = false here.  T5-T6 use a P_ms copy with
%  use_multi_start = true.
% =========================================================================

P.c       = 3e8;
P.c0      = P.c;           % Lesson from Task 11.2: both aliases required
P.fc      = 28e9;
P.lambda  = P.c / P.fc;
P.lambda_c = P.lambda;     % Lesson from Task 11.2: both aliases required

P.M       = 64;
P.d_ant   = P.lambda / 2;
P.N_RF    = 8;
P.N       = 64;
P.d       = 1;

P.theta_lo = 20 * pi/180;
P.theta_hi = 60 * pi/180;

% Range bounds (from nf_params v8 + nf_update_derived logic)
D_ap    = (P.M - 1) * P.d_ant;
P.r_RD  = 2 * D_ap^2 / P.lambda_c;
P.r_lo_fac  = 0.05;
P.r_hi_fac  = 1.0;
P.u_margin  = 2.0;
r_min   = P.r_lo_fac * P.r_RD;
r_max   = P.r_hi_fac * P.r_RD;
P.u_min = 1 / (r_max * P.u_margin);
P.u_max = P.u_margin / r_min;

% Wideband parameters: K_s=16, B=400 MHz.
%   Delta_f = B / K_s = 400 MHz / 16 = 25 MHz.
%   This gives alpha_k spread of ~0.013, sufficient for range discrimination
%   via cross-subcarrier quadratic phase (kappa) diversity.
P.K        = 128;
P.K_s      = 16;
P.B        = 400e6;                          % target bandwidth [Hz]
P.Delta_f  = P.B / P.K_s;                   % 25 MHz per subcarrier
k_idx_half = (-(P.K_s/2) : (P.K_s/2 - 1)).'; % K_s x 1, centred on dc
P.k_indices    = k_idx_half + P.K/2 + 1;    % 1-based subcarrier indices
P.alpha_k_vec  = 1 + k_idx_half * P.Delta_f / P.fc;  % f_k/fc, K_s x 1

% Fresnel coherence parameter
P.beta_delta = 1.2;

% CL-KL hyper-parameters (from nf_params)
P.lambda_reg = 1e-3;
P.max_iter   = 150;
P.tol_clkl   = 5e-4;
P.alpha_p    = 1.0;
P.ls_beta    = 0.5;
P.ls_sigma   = 1e-4;

% BPD parameters
P.G_theta = 128;
P.G_r     = 64;
P.d_max   = P.d;     % d=1 for toy test, so d_max=d (no truncation needed)
P.r_EBRD_fac = 0.058;

% CSV output control
P.wb_gen_write_csv   = false;   % suppress generator CSV in toy test
P.bpd_write_csv      = false;   % suppress BPD CSV

% v2 multi-start defaults (explicit; T1-T4 use single-start)
P.use_multi_start = false;
P.n_ms_starts     = 4;
P.short_iter_ms   = 20;

SNR_dB   = 10;
use_exact = true;

fprintf('Scene: M=%d, N=%d, N_RF=%d, K_s=%d, d=%d, SNR=%d dB, seed 42\n', ...
    P.M, P.N, P.N_RF, P.K_s, P.d, SNR_dB);

% =========================================================================
%  SECTION 2: GENERATE DATA
%  Single rng(42) call; one wb_channel_gen_ofdm_nf call supplies both
%  X_full (for BPD) and Y_full (for driver).
%  IMPORTANT: do NOT re-call the generator with a fresh seed (spec Sec. 8.1).
% =========================================================================

rng(42, 'twister');
[X_full, Y_full, ~, theta_true, r_true, p_true, N0_true, W_comb] = ...
    wb_channel_gen_ofdm_nf(P, SNR_dB, use_exact);

fprintf('True path: theta=%.4f deg, r=%.4f m\n', ...
    theta_true(1)*180/pi, r_true(1));

% -- Form compressed sample covariances (raw, NOT whitened)
R_hat_cell = cell(P.K_s, 1);
for k = 1:P.K_s
    Yk            = Y_full(:, :, k);
    R_hat_cell{k} = (1/P.N) * (Yk * Yk');
    R_hat_cell{k} = (R_hat_cell{k} + R_hat_cell{k}') / 2;
end

% =========================================================================
%  SECTION 3: BPD WARM-START
% =========================================================================

[theta_bpd_full, r_bpd_full, ~, ~] = bpd_baseline(X_full, P);

% Truncate to d entries (d_max = d = 1 here, so no truncation needed;
% but the pattern is explicit for the general case)
theta_bpd = theta_bpd_full(1:P.d);
r_bpd     = r_bpd_full(1:P.d);

theta_err_bpd = abs(theta_bpd(1) - theta_true(1)) * 180/pi;
r_err_bpd     = abs(r_bpd(1)     - r_true(1));

fprintf('\n--- BPD warm-start ---\n');
fprintf('  theta_bpd = %.4f deg  (err = %.4f deg)\n', theta_bpd(1)*180/pi, theta_err_bpd);
fprintf('  r_bpd     = %.4f m    (err = %.4f m)\n',   r_bpd(1),            r_err_bpd);

% =========================================================================
%  SECTION 4: WB-CL-KL DRIVER -- SINGLE-START (T1-T4 baseline)
%  P.use_multi_start = false ensures v1-equivalent behaviour for the
%  regression block.
% =========================================================================

assert(P.use_multi_start == false, ...
    'Section 4 requires P.use_multi_start = false for T1-T4 regression.');

[theta_hat, r_hat, p_hat, N0_hat, info] = ...
    wb_clkl_driver(R_hat_cell, W_comb, theta_bpd, r_bpd, [], P);

theta_err_driver = abs(theta_hat(1) - theta_true(1)) * 180/pi;
r_err_driver     = abs(r_hat(1)     - r_true(1));

fprintf('\n--- WB-CL-KL driver output (single-start) ---\n');
fprintf('  theta_hat = %.4f deg  (err = %.4f deg)\n', theta_hat(1)*180/pi, theta_err_driver);
fprintf('  r_hat     = %.4f m    (err = %.4f m)\n',   r_hat(1),            r_err_driver);
fprintf('  p_hat     = %.6f\n',   p_hat(1));
fprintf('  N0_hat    = %.4e  (true N0 = %.4e)\n', N0_hat, N0_true);
fprintf('  Converged : %s  (n_iter = %d)\n', ...
    mat2str(info.converged), info.n_iter);
if numel(info.L_hist) >= 1
    fprintf('  L_hist[1] = %.6f   L_hist[end] = %.6f\n', ...
        info.L_hist(1), info.L_hist(end));
end
fprintf('  best_start (single-start mode) = %d (expect 1)\n', info.best_start);

% =========================================================================
%  SECTION 5: ACCEPTANCE TESTS T1-T4
% =========================================================================

% -- T1: Convergence
T1_pass = info.converged && (info.n_iter <= P.max_iter);

fprintf('\n--- T1: Convergence check ---\n');
fprintf('  info.converged = %s, n_iter = %d (max %d)\n', ...
    mat2str(info.converged), info.n_iter, P.max_iter);
if T1_pass; fprintf('  PASS\n'); else; fprintf('  FAIL\n'); end

% -- T2: Monotone decrease (Armijo guarantee)
if numel(info.L_hist) >= 2
    dL = diff(info.L_hist(1:info.n_iter));
    max_increase = max(dL);
else
    dL = [];
    max_increase = 0;
end
T2_pass = isempty(dL) || (max_increase <= 1e-9);

fprintf('\n--- T2: Monotone decrease ---\n');
fprintf('  Max dL increase: %.3e\n', max_increase);
if T2_pass; fprintf('  PASS\n'); else; fprintf('  FAIL\n'); end

% -- T3: Improvement over BPD warm-start
T3_theta_pass = (theta_err_driver < theta_err_bpd) || (theta_err_driver < 0.1);
T3_r_pass     = (r_err_driver     < r_err_bpd)     || (r_err_driver     < 0.05);

fprintf('\n--- T3: Improvement over warm-start ---\n');
fprintf('  theta: BPD %.4f deg -> driver %.4f deg  %s\n', ...
    theta_err_bpd, theta_err_driver, pass_str(T3_theta_pass));
fprintf('  r:     BPD %.4f m   -> driver %.4f m   %s\n', ...
    r_err_bpd, r_err_driver, pass_str(T3_r_pass));

% -- T4: N0 sanity
T4_pass = (N0_hat > N0_true/3) && (N0_hat < N0_true*3);

fprintf('\n--- T4: N0 sanity ---\n');
fprintf('  N0_hat = %.4e  in (N0_true/3=%.4e, 3*N0_true=%.4e)?  %s\n', ...
    N0_hat, N0_true/3, N0_true*3, pass_str(T4_pass));

% =========================================================================
%  SECTION 6: WB-CL-KL DRIVER -- MULTI-START (T5-T6 verification)
%
%  Per starter prompt Sec.6.2 Changes 2-3, T5/T6 are evaluated with
%  (theta_true, r_true) as warm-start.  This is a controlled-environment
%  check: with the truth as input, Candidate 1 (BPD slot, here = truth)
%  starts at the optimum, and the alternatives test that multi-start
%  does not regress relative to single-start.
% =========================================================================

P_ms = P;
P_ms.use_multi_start = true;
P_ms.n_ms_starts     = 4;
P_ms.short_iter_ms   = 20;

fprintf('\n--- Multi-start driver call (T5/T6, theta_true/r_true input) ---\n');

% Multi-start call
[~, ~, ~, ~, info_ms] = wb_clkl_driver( ...
    R_hat_cell, W_comb, theta_true, r_true, [], P_ms);

% --- v2 driver detection (catches v1-on-path silently) ----------------
% v2's info struct contains L_hist_ms (empty when use_multi_start=false,
% short_iter_ms x n_ms_starts matrix when true).  v1 has neither field.
% If this assertion fires, MATLAB resolved wb_clkl_driver to a file named
% wb_clkl_driver.m on the path that is NOT the v2 file.  Fix: rename
% wb_clkl_driver_v2.m -> wb_clkl_driver.m (archive the v1 file first).
assert(isfield(info_ms, 'L_hist_ms'), ...
    ['wb_clkl_driver_toy_test_v2:WrongDriverVersion -- ' ...
     'Driver returned info without L_hist_ms.  v1 appears to be on the ' ...
     'path.  Rename wb_clkl_driver_v2.m to wb_clkl_driver.m (archive ' ...
     'v1 first), then run "clear functions" before re-running.']);

% Single-start call (P.use_multi_start = false from Section 1)
[~, ~, ~, ~, info_ss] = wb_clkl_driver( ...
    R_hat_cell, W_comb, theta_true, r_true, [], P);

L_ss_final = info_ss.L_hist(end);
L_ms_final = info_ms.L_hist(end);

fprintf('  single-start: converged=%s, n_iter=%d, L_final=%.6f\n', ...
    mat2str(info_ss.converged), info_ss.n_iter, L_ss_final);
fprintf('  multi-start : converged=%s, n_iter=%d, L_final=%.6f, best_start=%d\n', ...
    mat2str(info_ms.converged), info_ms.n_iter, L_ms_final, info_ms.best_start);

% Brief peek at the L_hist_ms matrix (initial vs short-trial-final per cand.)
if ~isempty(info_ms.L_hist_ms)
    L_ms_init  = nan(1, P_ms.n_ms_starts);
    L_ms_end   = nan(1, P_ms.n_ms_starts);
    for si = 1:P_ms.n_ms_starts
        col = info_ms.L_hist_ms(:, si);
        col = col(~isnan(col));
        if ~isempty(col)
            L_ms_init(si) = col(1);
            L_ms_end(si)  = col(end);
        end
    end
    fprintf('  L_hist_ms initial (per candidate) :');
    fprintf(' %.4f', L_ms_init); fprintf('\n');
    fprintf('  L_hist_ms final   (per candidate) :');
    fprintf(' %.4f', L_ms_end);  fprintf('\n');
end

% =========================================================================
%  SECTION 7: ACCEPTANCE TESTS T5-T6
% =========================================================================

% -- T5: best_start_idx is in valid range {1,...,n_ms_starts}
T5_pass = (info_ms.best_start >= 1) && (info_ms.best_start <= P_ms.n_ms_starts);

fprintf('\n--- T5: Multi-start activation ---\n');
fprintf('  best_start_idx = %d (expect in {1,...,%d})  %s\n', ...
    info_ms.best_start, P_ms.n_ms_starts, pass_str(T5_pass));

% -- T6: multi-start non-regression vs single-start
T6_pass = (L_ms_final <= L_ss_final + 1e-6);

fprintf('\n--- T6: Multi-start non-regression ---\n');
fprintf('  L_ss_final = %.6f, L_ms_final = %.6f, gap = %+.3e  %s\n', ...
    L_ss_final, L_ms_final, L_ms_final - L_ss_final, pass_str(T6_pass));

% =========================================================================
%  SECTION 8: SUMMARY
% =========================================================================

all_pass = T1_pass && T2_pass && T3_theta_pass && T3_r_pass && ...
           T4_pass && T5_pass && T6_pass;

fprintf('\n=============================================================\n');
fprintf(' SUMMARY:\n');
fprintf('   T1 %s | T2 %s | T3-theta %s | T3-r %s | T4 %s\n', ...
    pass_str(T1_pass), pass_str(T2_pass), ...
    pass_str(T3_theta_pass), pass_str(T3_r_pass), pass_str(T4_pass));
fprintf('   T5 %s | T6 %s\n', ...
    pass_str(T5_pass), pass_str(T6_pass));
if all_pass
    fprintf(' wb_clkl_driver (v2, Component 2) is verified.\n');
else
    fprintf(' wb_clkl_driver (v2) has FAILURES -- check output above.\n');
end
fprintf('=============================================================\n');

% =========================================================================
%  SECTION 9: CSV DIAGNOSTIC OUTPUT
% =========================================================================

csv_path = 'wb_clkl_driver_toy_test_v2_diag.csv';
fid = fopen(csv_path, 'w');
fprintf(fid, 'field,value\n');
fprintf(fid, 'theta_true_deg,%.8f\n',  theta_true(1)*180/pi);
fprintf(fid, 'r_true_m,%.8f\n',        r_true(1));
fprintf(fid, 'N0_true,%.8e\n',         N0_true);
fprintf(fid, 'theta_bpd_deg,%.8f\n',   theta_bpd(1)*180/pi);
fprintf(fid, 'r_bpd_m,%.8f\n',         r_bpd(1));
fprintf(fid, 'theta_err_bpd_deg,%.8f\n', theta_err_bpd);
fprintf(fid, 'r_err_bpd_m,%.8f\n',     r_err_bpd);
fprintf(fid, 'theta_hat_deg,%.8f\n',   theta_hat(1)*180/pi);
fprintf(fid, 'r_hat_m,%.8f\n',         r_hat(1));
fprintf(fid, 'theta_err_driver_deg,%.8f\n', theta_err_driver);
fprintf(fid, 'r_err_driver_m,%.8f\n',  r_err_driver);
fprintf(fid, 'p_hat,%.8f\n',           p_hat(1));
fprintf(fid, 'N0_hat,%.8e\n',          N0_hat);
fprintf(fid, 'converged,%d\n',         info.converged);
fprintf(fid, 'n_iter,%d\n',            info.n_iter);
if numel(info.L_hist) >= 1
    fprintf(fid, 'L_hist_first,%.8f\n', info.L_hist(1));
    fprintf(fid, 'L_hist_last,%.8f\n',  info.L_hist(end));
end
% v2 additions
fprintf(fid, 'best_start_idx,%d\n', info_ms.best_start);
fprintf(fid, 'n_ms_starts,%d\n',    P_ms.n_ms_starts);
fprintf(fid, 'short_iter_ms,%d\n',  P_ms.short_iter_ms);
fprintf(fid, 'L_ss_final,%.8f\n',   L_ss_final);
fprintf(fid, 'L_ms_final,%.8f\n',   L_ms_final);
fprintf(fid, 'L_gap_ms_minus_ss,%.8e\n', L_ms_final - L_ss_final);
fprintf(fid, 'T1_pass,%d\n',       T1_pass);
fprintf(fid, 'T2_pass,%d\n',       T2_pass);
fprintf(fid, 'T3_theta_pass,%d\n', T3_theta_pass);
fprintf(fid, 'T3_r_pass,%d\n',     T3_r_pass);
fprintf(fid, 'T4_pass,%d\n',       T4_pass);
fprintf(fid, 'T5_pass,%d\n',       T5_pass);
fprintf(fid, 'T6_pass,%d\n',       T6_pass);
fclose(fid);
fprintf('\nDiagnostic CSV written to: %s\n', csv_path);

% =========================================================================
%  LOCAL FUNCTION: pass_str
% =========================================================================
function s = pass_str(flag)
if flag; s = 'PASS'; else; s = 'FAIL'; end
end
