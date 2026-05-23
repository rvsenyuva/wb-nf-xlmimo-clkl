%WB_CLKL_TOY_TEST  Standalone verification suite for wb_clkl_estimator.m (v2).
%
%  Paper C Phase 2, Task 11.3 (rev2) toy test.
%
%  Run from the MATLAB command window:   wb_clkl_toy_test
%
%  Scene: M=64, N=16, N_RF=8, K_s=16, d=1, SNR=10 dB, seed 42
%         use_exact = true  (USW truth model, Branch B)
%
%  Five test blocks:
%    T1 -- Output shape check: L is real finite scalar; grad_L is [(3d+1) x 1] real
%          NOTE: L = log det(R) + trace(R^{-1} R_hat) can be negative (no sign
%          constraint on this two-term form); correct check is isfinite, not L>0.
%    T2 -- Finite-difference gradient check (eps=1e-6); pass: max_rel_err < 1e-4
%    T3 -- Symmetry: L(eta_true) < L(eta_true + 0.01*randn(...))
%    T4 -- N0 gradient sign: gradient > 0 when N0 is underestimated
%    T5 -- Conditioning gate: min rcond(Ry_reg) > 1e-10 at all SNR in
%          [-5, 0, 5, 10, 15, 20] dB; no RCOND singularity warnings.
%          This test is the PRIMARY VERIFICATION for Task 11.3 (rev2).
%
%  Expected output: all five tests PASS.
%  Key numerical results: T2 max relative gradient error, T5 min RCOND.
%
%  Author : R. V. Senyuva (Maltepe University)
%  Date   : May 2026 (v2: T5 added for Task 11.3 rev2 conditioning gate)

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' wb_clkl_toy_test -- Task 11.3 (rev2) verification\n');
fprintf('=============================================================\n');

% =========================================================================
%  BUILD PARAMETER STRUCT
% =========================================================================
P = toy_params_clkl();

fprintf('Scene: M=%d, N=%d, N_RF=%d, K_s=%d, d=%d\n', ...
    P.M, P.N, P.N_RF, P.K_s, P.d);

% =========================================================================
%  GENERATE CHANNEL DATA  (SNR=10 dB for T1-T4)
% =========================================================================
rng(42, 'twister');

[~, Y_full, ~, theta_true, r_true, p_true, N0_true, W_comb] = ...
    wb_channel_gen_ofdm_nf(P, 10, true);

% Form raw compressed sample covariances  (no whitening, Design Decision D3)
R_hat_cell = cell(P.K_s, 1);
for k = 1:P.K_s
    Yk              = Y_full(:, :, k);
    R_hat_cell{k}   = (1/P.N) * (Yk * Yk');
    R_hat_cell{k}   = (R_hat_cell{k} + R_hat_cell{k}') / 2;   % symmetrise
end

% =========================================================================
%  BUILD eta_true FROM GROUND TRUTH
% =========================================================================
omega_true = (2*pi*P.d_ant/P.lambda_c) * cos(theta_true);     % d x 1
kappa_true = (pi*P.d_ant^2/P.lambda_c) * sin(theta_true).^2 ./ r_true;  % d x 1
eta_true   = [omega_true; kappa_true; p_true; N0_true];        % (3d+1) x 1

fprintf('eta_true:\n');
fprintf('  omega = %.8f\n', omega_true);
fprintf('  kappa = %.8f\n', kappa_true);
fprintf('  p     = %.8f\n', p_true);
fprintf('  N0    = %.4e\n', N0_true);
fprintf('\n');

% =========================================================================
%  EVALUATE OBJECTIVE AND GRADIENT AT eta_true
% =========================================================================
[L_true, grad_L_true, ~] = wb_clkl_estimator(eta_true, R_hat_cell, W_comb, P);

% =========================================================================
%  T1 -- Output shape check
% =========================================================================
fprintf('--- T1: Output shape check ---\n');

t1_L_pass    = isscalar(L_true) && isreal(L_true) && isfinite(L_true) && ~isinf(L_true);
t1_grad_pass = isequal(size(grad_L_true), [3*P.d+1, 1]) && isreal(grad_L_true);

if t1_L_pass
    fprintf('  PASS L: scalar, real, finite  (L = %.6f)\n', L_true);
else
    fprintf('  FAIL L: scalar=%d, real=%d, finite=%d  (L = %g)\n', ...
        isscalar(L_true), isreal(L_true), isfinite(L_true), L_true);
end

if t1_grad_pass
    fprintf('  PASS grad_L: [%d x 1] real\n', numel(grad_L_true));
else
    fprintf('  FAIL grad_L: size=[%s], isreal=%d\n', ...
        num2str(size(grad_L_true)), isreal(grad_L_true));
end
t1_pass = t1_L_pass && t1_grad_pass;

% =========================================================================
%  T2 -- Finite-difference gradient check
% =========================================================================
fprintf('\n--- T2: Finite-difference gradient check (eps=1e-6) ---\n');

eps_fd   = 1e-6;
n_eta    = numel(eta_true);
grad_fd  = zeros(n_eta, 1);

for i = 1:n_eta
    ep    = eta_true;  ep(i) = ep(i) + eps_fd;
    em    = eta_true;  em(i) = em(i) - eps_fd;
    Lp    = wb_clkl_estimator(ep, R_hat_cell, W_comb, P);
    Lm    = wb_clkl_estimator(em, R_hat_cell, W_comb, P);
    grad_fd(i) = (Lp - Lm) / (2 * eps_fd);
end

rel_err     = abs(grad_L_true - grad_fd) ./ (abs(grad_fd) + 1e-10);
max_rel_err = max(rel_err);

% Build component labels
comp_labels = cell(n_eta, 1);
for l = 1:P.d
    comp_labels{l}           = sprintf('omega_%d', l);
    comp_labels{P.d + l}     = sprintf('kappa_%d', l);
    comp_labels{2*P.d + l}   = sprintf('p_%d',     l);
end
comp_labels{3*P.d + 1} = 'N0';

fprintf('  %-12s  %-14s  %-14s  %-12s\n', ...
    'Component', 'analytic', 'FD', 'rel_err');
fprintf('  %-12s  %-14s  %-14s  %-12s\n', ...
    repmat('-',1,12), repmat('-',1,14), repmat('-',1,14), repmat('-',1,12));
for i = 1:n_eta
    fprintf('  %-12s  %+14.6e  %+14.6e  %.4e\n', ...
        comp_labels{i}, grad_L_true(i), grad_fd(i), rel_err(i));
end
fprintf('\n  Max relative error: %.4e\n', max_rel_err);

t2_pass = (max_rel_err < 1e-4);
if t2_pass
    fprintf('  PASS (< 1e-4)\n');
else
    fprintf('  FAIL (>= 1e-4)  -- investigate gradient formula\n');
end

% =========================================================================
%  T3 -- Symmetry: L(eta_true) should be less than L at a perturbed point
% =========================================================================
fprintf('\n--- T3: Symmetry ---\n');

rng(123, 'twister');
eta_perturb = eta_true + 0.01 * randn(n_eta, 1);
% Ensure positive p and N0 after perturbation
eta_perturb(2*P.d+1 : 3*P.d) = max(eta_perturb(2*P.d+1 : 3*P.d), 0);
eta_perturb(3*P.d+1)          = max(eta_perturb(3*P.d+1),          1e-12);

L_perturb = wb_clkl_estimator(eta_perturb, R_hat_cell, W_comb, P);

fprintf('  L(eta_true)         = %.6f\n', L_true);
fprintf('  L(eta_true+0.01*n)  = %.6f  (should be larger)\n', L_perturb);

t3_pass = (L_true < L_perturb);
if t3_pass
    fprintf('  PASS\n');
else
    fprintf('  FAIL  -- L increased at true parameters; check signal model\n');
end

% =========================================================================
%  T4 -- N0 gradient sign: gradient positive when N0 underestimated
% =========================================================================
fprintf('\n--- T4: N0 gradient sign ---\n');

eta_N0_low       = eta_true;
eta_N0_low(end)  = eta_true(end) * 0.5;   % halve N0
L_N0_low = wb_clkl_estimator(eta_N0_low, R_hat_cell, W_comb, P);

fprintf('  grad_L(N0) at eta_true   = %+.6e\n', grad_L_true(end));
fprintf('  L(eta_true)              = %.6f\n', L_true);
fprintf('  L(N0 * 0.5)             = %.6f  (should be larger)\n', L_N0_low);

t4_pass = (grad_L_true(end) > 0) && (L_N0_low > L_true);
if t4_pass
    fprintf('  PASS  (gradient positive AND L increases when N0 halved)\n');
else
    if grad_L_true(end) <= 0
        fprintf('  FAIL  (gradient non-positive: %.4e)\n', grad_L_true(end));
    else
        fprintf('  FAIL  (L did not increase when N0 halved; check noise term)\n');
    end
end

% =========================================================================
%  T5 -- Conditioning gate: min rcond(Ry_reg) > 1e-10 across SNR range
%        PRIMARY VERIFICATION for Task 11.3 (rev2)
% =========================================================================
fprintf('\n--- T5: Conditioning gate (SNR-adaptive regularisation) ---\n');
fprintf('  Criterion: min rcond(Ry_reg_k) > 1e-10 at ALL subcarriers k\n');
fprintf('             for ALL SNR in [-5, 0, 5, 10, 15, 20] dB\n');
fprintf('  P.eps_reg = %.4g (default 1e-3)\n', P.eps_reg);

snr_t5_vec    = [-5, 0, 5, 10, 15, 20];   % dB
t5_pass_all   = true;
min_rcond_all = Inf;
min_rcond_snr = NaN;   % SNR at which minimum RCOND occurs

WtW_t5  = W_comb' * W_comb;    % N_RF x N_RF (reuse from above scope)
I_t5    = eye(P.N_RF);

for snr_t5 = snr_t5_vec
    % Generate fresh channel at this SNR
    rng(42, 'twister');   % fixed seed for reproducibility
    [~, Y_t5, ~, theta_t5, r_t5, p_t5, N0_t5, ~] = ...
        wb_channel_gen_ofdm_nf(P, snr_t5, true);

    % Form compressed sample covariances
    Rhat_t5 = cell(P.K_s, 1);
    for k = 1:P.K_s
        Yk_t5        = Y_t5(:, :, k);
        Rhat_t5{k}   = (1/P.N) * (Yk_t5 * Yk_t5');
        Rhat_t5{k}   = (Rhat_t5{k} + Rhat_t5{k}') / 2;
    end

    % Build eta at true parameters
    om_t5  = (2*pi*P.d_ant/P.lambda_c) * cos(theta_t5);
    kp_t5  = (pi*P.d_ant^2/P.lambda_c) * sin(theta_t5).^2 ./ r_t5;
    eta_t5 = [om_t5; kp_t5; p_t5; N0_t5];

    % Rebuild model covariance at eta_true and compute Ry_reg per subcarrier
    %   Replicate Phase A + inner loop from wb_clkl_estimator to get Ry_reg
    c_lin_t5  = 2 * pi * P.d_ant / P.lambda_c;
    c_quad_t5 =     pi * P.d_ant^2 / P.lambda_c;
    arg_t5    = om_t5 / c_lin_t5;
    arg_t5    = min(1-1e-6, max(1e-6, arg_t5));
    th_t5     = acos(arg_t5);
    sin2_t5   = sin(th_t5).^2;
    u_t5      = kp_t5 ./ (c_quad_t5 * sin2_t5);
    u_t5      = min(P.u_max, max(P.u_min, u_t5));
    m_bar_t5  = ((0:P.M-1).' - (P.M-1)/2);

    min_rcond_snr_k = Inf;
    for k = 1:P.K_s
        ak_t5 = zeros(P.M, P.d);
        for l = 1:P.d
            ak_t5(:,l) = wb_nf_fresnel_steer(th_t5(l), u_t5(l), ...
                P.alpha_k_vec(k), P) * sqrt(P.M);
        end
        Dk_t5  = W_comb' * ak_t5;
        Ry_t5  = Dk_t5 * diag(max(p_t5, 0)) * Dk_t5' + max(N0_t5, 1e-12) * WtW_t5;
        Ry_t5  = (Ry_t5 + Ry_t5') / 2;
        eps_t5 = P.eps_reg * trace(Rhat_t5{k}) / P.N_RF;
        Ry_reg_t5 = Ry_t5 + eps_t5 * I_t5;
        rc_k = rcond(Ry_reg_t5);
        if rc_k < min_rcond_snr_k
            min_rcond_snr_k = rc_k;
        end
    end

    if min_rcond_snr_k < min_rcond_all
        min_rcond_all = min_rcond_snr_k;
        min_rcond_snr = snr_t5;
    end

    if min_rcond_snr_k < 1e-10
        fprintf('  SNR=%+3d dB: min_rcond=%.2e -- FAIL\n', snr_t5, min_rcond_snr_k);
        t5_pass_all = false;
    else
        fprintf('  SNR=%+3d dB: min_rcond=%.2e -- PASS\n', snr_t5, min_rcond_snr_k);
    end
end

t5_pass = t5_pass_all;
if t5_pass
    fprintf('  T5 PASS: min RCOND = %.2e at SNR=%+d dB (threshold 1e-10)\n', ...
        min_rcond_all, min_rcond_snr);
else
    fprintf('  T5 FAIL: min RCOND = %.2e -- below 1e-10 threshold.\n', min_rcond_all);
    fprintf('  Action: increase P.eps_reg (try 1e-2) and re-run T5.\n');
end

% =========================================================================
%  WRITE DIAGNOSTIC CSV
% =========================================================================
diag_file = 'wb_clkl_toy_test_v2_diag.csv';
fid = fopen(diag_file, 'w');
fprintf(fid, 'field,value\n');
fprintf(fid, 'max_rel_grad_err_T2,%.6e\n', max_rel_err);
fprintf(fid, 'L_true,%.6f\n',              L_true);
fprintf(fid, 'L_perturb,%.6f\n',           L_perturb);
fprintf(fid, 'grad_N0_true,%.6e\n',        grad_L_true(end));
fprintf(fid, 'L_N0_low,%.6f\n',            L_N0_low);
fprintf(fid, 'min_rcond_all,%.6e\n',       min_rcond_all);
fprintf(fid, 'min_rcond_snr_dB,%d\n',      min_rcond_snr);
fprintf(fid, 'eps_reg,%.6e\n',             P.eps_reg);
fprintf(fid, 'T1_pass,%d\n', double(t1_pass));
fprintf(fid, 'T2_pass,%d\n', double(t2_pass));
fprintf(fid, 'T3_pass,%d\n', double(t3_pass));
fprintf(fid, 'T4_pass,%d\n', double(t4_pass));
fprintf(fid, 'T5_pass,%d\n', double(t5_pass));
fclose(fid);
fprintf('\nDiagnostic CSV written: %s\n', diag_file);

% =========================================================================
%  SUMMARY
% =========================================================================
pass_all = t1_pass && t2_pass && t3_pass && t4_pass && t5_pass;
pass_str = @(b) subsref({'FAIL','PASS'}, struct('type','{}','subs',{{1+b}}));

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' SUMMARY T1-T5: T1 %s | T2 %s | T3 %s | T4 %s | T5 %s\n', ...
    pass_str(t1_pass), pass_str(t2_pass), pass_str(t3_pass), ...
    pass_str(t4_pass), pass_str(t5_pass));
if pass_all
    fprintf(' wb_clkl_estimator v2 verified.\n');
else
    fprintf(' ATTENTION: one or more tests FAILED.  See details above.\n');
end
fprintf(' Max relative gradient error (T2): %.4e\n', max_rel_err);
fprintf(' Min RCOND across SNR sweep  (T5): %.4e  (at SNR=%+d dB)\n', ...
    min_rcond_all, min_rcond_snr);
fprintf('=============================================================\n\n');


% =========================================================================
%  LOCAL HELPER: toy_params_clkl
% =========================================================================

function P = toy_params_clkl()
%TOY_PARAMS_CLKL  Minimal parameter struct for wb_clkl_toy_test (v2).
%
%  Scene: M=64, N=16, N_RF=8, K_s=16, d=1
%  System: fc=28 GHz, d_ant=lambda_c/2, B=16*Delta_f (Delta_f=120 kHz)
%
%  Follows Task 11.2 lesson: set BOTH P.c0 and P.lambda_c (and their
%  aliases P.c and P.lambda) to avoid field-not-found errors downstream.
%  v2 addition: P.eps_reg = 1e-3 (Task 11.3 rev2 SNR-adaptive loading).

P = struct();

% --- Physical constants -------------------------------------------------
P.c0        = 3e8;
P.c         = P.c0;                   % Paper B alias (keep both)
P.fc        = 28e9;

% --- Wavelengths (BOTH aliases required, Task 11.2 lesson) -------------
P.lambda_c  = P.c0 / P.fc;            % Paper C canonical name
P.lambda    = P.lambda_c;             % Paper B alias: must equal lambda_c

% --- Array --------------------------------------------------------------
P.M         = 64;
P.d_ant     = P.lambda_c / 2;         % half-wavelength spacing

% --- Hybrid architecture -----------------------------------------------
P.N_RF      = 8;

% --- Snapshots ----------------------------------------------------------
P.N         = 16;

% --- Channel ------------------------------------------------------------
P.d         = 1;
P.theta_lo  = 20 * pi/180;
P.theta_hi  = 60 * pi/180;

% --- Range bounds (required by wb_channel_gen_ofdm_nf) -----------------
P.r_lo_fac  = 0.05;
P.r_hi_fac  = 0.40;
P.D_ap      = (P.M - 1) * P.d_ant;
P.r_RD      = 2 * P.D_ap^2 / P.lambda_c;
P.u_margin  = 2.0;
r_min       = P.r_lo_fac * P.r_RD;
r_max       = P.r_hi_fac * P.r_RD;
P.u_min     = 1 / (r_max * P.u_margin);
P.u_max     = P.u_margin / r_min;

% --- Wideband subcarrier parameters ------------------------------------
P.K         = 2048;                   % total OFDM subcarriers
P.Delta_f   = 120e3;                  % subcarrier spacing [Hz] (5G NR FR2)
P.K_s       = 16;                     % used subcarriers (toy: 16)
% Centre K_s subcarriers symmetrically around fc
k0          = round(P.K / 2) - floor(P.K_s / 2) + 1;
P.k_indices = (k0 : k0 + P.K_s - 1).';   % K_s x 1
f_k_vec     = (P.k_indices - 1) * P.Delta_f;
P.B         = P.K_s * P.Delta_f;          % occupied bandwidth [Hz]
P.alpha_k_vec = (P.fc + f_k_vec) / P.fc; % K_s x 1  (f_k/fc, near 1.0)

% --- Optimiser hyper-parameters ----------------------------------------
P.lambda_reg = 1e-3;
P.max_iter   = 150;
P.tol_clkl   = 5e-4;
P.alpha_p    = 1.0;
P.ls_beta    = 0.5;
P.ls_sigma   = 1e-4;

% --- SNR-adaptive regularisation (Task 11.3 rev2) ----------------------
P.eps_reg    = 1e-3;   % default; range [1e-4, 1e-3]

% --- CSV diagnostics (disable for toy test) ----------------------------
P.wb_gen_write_csv = false;

end  % toy_params_clkl
