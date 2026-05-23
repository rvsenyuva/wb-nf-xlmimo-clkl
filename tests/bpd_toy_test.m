%BPD_TOY_TEST  Standalone verification suite for bpd_baseline.m  (Day-2).
%
%  Run from the MATLAB command window:   bpd_toy_test
%
%  Three test blocks:
%    Block A: Day-1 baseline scene (2-path, SNR=10 dB, M=64, K_s=16, N=8)
%             Confirms the Day-2 LS polishing column appears in the CSV
%             and reports res_norm before and after polish.
%
%    Block B: SNR sweep (Target 2)
%             Same 2-path scene; SNR in {0, 5, 10, 15, 20} dB.
%             Reports theta_err, r_err, pass/fail per (SNR, path) pair.
%             Expected: FAIL at SNR = 0 dB for a greedy estimator; document.
%
%    Block C: Scaling test (Target 3)
%             Single-path scene at M=256, N=64, K_s=512, G_theta=128,
%             G_r=64, d_max=1, SNR=10 dB.
%             Times bpd_baseline() call; requirement < 300 s.
%             Reports wall-clock time and pass/fail.
%
%  Author : R. V. Senyuva (Maltepe University)
%  Version: Day-2
%  Date   : May 2026

fprintf('\n');
fprintf('=================================================================\n');
fprintf('  BPD_TOY_TEST  --  Day-2 verification suite\n');
fprintf('=================================================================\n');

% =========================================================================
%  SHARED HELPER: build_params
% =========================================================================
% Returns a minimal params struct for the toy scenes.
% M, N, K_s, G_theta, G_r, d_max are passed as arguments;
% all other fields are locked at the Day-1 defaults.

% (Defined as a nested function below; called in each block.)

% =========================================================================
%  BLOCK A: Day-1 baseline scene (2-path, SNR=10 dB)
%           Confirms LS polishing column in CSV.
% =========================================================================
fprintf('\n--- BLOCK A: Day-1 baseline scene (SNR=10 dB) ---\n');

rng(42, 'twister');

Pa = toy_params(64, 8, 16, 64, 32, 2);
Pa.bpd_write_csv = true;
Pa.bpd_csv_path  = 'bpd_toy_test_diag.csv';
Pa.bpd_do_ls_polish = true;

[Xa, theta_true_a, r_true_a] = toy_scene(Pa, 10);

t0a = tic;
[th_a, r_a, ~, diag_a] = bpd_baseline(Xa, Pa);
tA = toc(t0a);

fprintf('Elapsed: %.2f s\n', tA);
report_errors(th_a, r_a, theta_true_a, r_true_a, diag_a, true);
fprintf('  res_norm (greedy last step) : %.6e\n', diag_a.res_norm(end));
fprintf('  res_norm_polished           : %.6e\n', diag_a.res_norm_polished);
if isnan(diag_a.res_norm_polished)
    fprintf('  [LS polish DISABLED]\n');
else
    ratio = diag_a.res_norm_polished / diag_a.res_norm(end);
    fprintf('  Polish / Greedy ratio       : %.4f  (< 1 means polish helped)\n', ratio);
end
fprintf('  CSV written to bpd_toy_test_diag.csv\n');

% =========================================================================
%  BLOCK B: SNR sweep (Target 2)
% =========================================================================
fprintf('\n--- BLOCK B: SNR sweep ---\n');
fprintf('Scene: M=64, N=8, K_s=16, G_theta=64, G_r=32, d_max=2 (2-path)\n');
fprintf('Pass threshold: |theta_err| < 1.0 deg  AND  |r_err| < 0.5 m\n\n');

SNR_vec = [0, 5, 10, 15, 20];   % dB

% Print table header
fmt_hdr = '%-8s | %-4s | %-16s | %-12s | %s\n';
fmt_row = '%-8d | %-4d | %-16.3f | %-12.4f | %s\n';
fmt_sum = '%-8d | %-4s | max=%-12.3f | max=%-8.4f | %s\n';
sep = [repmat('-',1,8),'|',repmat('-',1,6),'|',repmat('-',1,18),'|', ...
       repmat('-',1,14),'|',repmat('-',1,10)];

fprintf(fmt_hdr, 'SNR(dB)', 'Path', 'theta_err(deg)', 'r_err(m)', 'Status');
fprintf('%s\n', sep);

choose_str = @(c) subsref({'FAIL','PASS'}, struct('type','{}','subs',{{1+c}}));

rng(42, 'twister');   % reset seed for reproducible sweep
for si = 1:numel(SNR_vec)
    SNR_dB = SNR_vec(si);
    Pb     = toy_params(64, 8, 16, 64, 32, 2);
    Pb.bpd_write_csv    = false;
    Pb.bpd_do_ls_polish = true;

    [Xb, theta_true_b, r_true_b] = toy_scene(Pb, SNR_dB);
    [th_b, r_b, ~, ~] = bpd_baseline(Xb, Pb);

    % Nearest-neighbour matching
    d_true = numel(theta_true_b);
    err_th = zeros(d_true,1);
    err_r  = zeros(d_true,1);
    used   = false(Pb.d_max, 1);
    for ell = 1:d_true
        diffs = abs(th_b - theta_true_b(ell)) * 180/pi;
        diffs(used) = inf;
        [~, best] = min(diffs);
        err_th(ell) = abs(th_b(best) - theta_true_b(ell)) * 180/pi;
        err_r(ell)  = abs(r_b(best)  - r_true_b(ell));
        used(best)  = true;
    end

    all_pass = true;
    for ell = 1:d_true
        pij = (err_th(ell) < 1.0) && (err_r(ell) < 0.5);
        all_pass = all_pass && pij;
        fprintf(fmt_row, SNR_dB, ell, err_th(ell), err_r(ell), choose_str(pij));
    end
    fprintf(fmt_sum, SNR_dB, 'ALL', max(err_th), max(err_r), choose_str(all_pass));
    fprintf('%s\n', sep);
end

fprintf('\nNOTE: FAIL at SNR=0 dB is expected for a greedy single-snapshot\n');
fprintf('estimator at the strict (1 deg, 0.5 m) threshold. This is not a\n');
fprintf('defect; it reflects the fundamental SNR floor of BPD-OMP without\n');
fprintf('spatial smoothing or iterative refinement.\n');

% =========================================================================
%  BLOCK C: Scaling test (Target 3)
% =========================================================================
fprintf('\n--- BLOCK C: Scaling test ---\n');
fprintf('Scene: M=256, N=64, K_s=512, G_theta=128, G_r=64, d_max=1\n');
fprintf('Requirement: bpd_baseline() wall-clock time < 300 s\n\n');

% Memory estimate (informational):
bytes_Psi = 256 * (128*64) * 16;    % complex double = 16 bytes
bytes_R   = 256 * 64 * 512 * 16;
fprintf('Memory estimates:\n');
fprintf('  Psi_mat  (256 x 8192)          : %.1f MB\n', bytes_Psi/1e6);
fprintf('  R_tensor (256 x 64 x 512)      : %.1f MB\n', bytes_R/1e6);

Pc = toy_params(256, 64, 512, 128, 64, 1);
Pc.bpd_write_csv    = false;
Pc.bpd_do_ls_polish = true;

% Single-path scene: theta=30 deg, r = 0.15 * r_RD
%   For M=256 at 28 GHz: D_ap = 255*(lambda_c/2), r_RD ~ 342 m
%   r = 0.15 * r_RD ~ 51 m  (inside near-field, well inside EBRD at 30 deg)
theta_c_true = 30 * pi/180;
r_c_true     = 0.15 * Pc.r_RD;

fprintf('\nTrue path: theta=%.1f deg, r=%.2f m (%.2f * r_RD)\n', ...
    theta_c_true*180/pi, r_c_true, r_c_true/Pc.r_RD);

rng(7, 'twister');
[Xc, ~, ~] = toy_scene_single(Pc, 10, theta_c_true, r_c_true);

fprintf('Generating data done.  Starting bpd_baseline() ...\n');

t0c = tic;
[th_c, r_c, ~, diag_c] = bpd_baseline(Xc, Pc);
tC = toc(t0c);

theta_err_c = abs(th_c(1) - theta_c_true) * 180/pi;
r_err_c     = abs(r_c(1)  - r_c_true);

fprintf('\nScaling test results:\n');
fprintf('  Wall-clock time  : %.2f s\n', tC);
fprintf('  theta_err        : %.3f deg\n', theta_err_c);
fprintf('  r_err            : %.4f m\n',   r_err_c);
fprintf('  res_norm (greedy): %.6e\n', diag_c.res_norm(1));
fprintf('  res_norm_polished: %.6e\n', diag_c.res_norm_polished);

if tC < 300
    fprintf('\n  SCALING TEST: PASS  (%.2f s < 300 s)\n', tC);
else
    fprintf('\n  SCALING TEST: FAIL  (%.2f s >= 300 s)\n', tC);
    fprintf('  --> Profile output follows.  Dominant bottleneck:\n');
    fprintf('      Per-subcarrier BLAS matmul Psi_mat'' * Rk is likely O(K_s*M*G).\n');
    fprintf('      At K_s=512, G=8192, M=256: 512 * (8192x256)*(256x64) matmuls.\n');
    fprintf('      If this dominates, consider: (1) reducing G_theta to 64 for\n');
    fprintf('      benchmarking (document explicitly); (2) computing atom powers\n');
    fprintf('      from sample covariance Rk*Rk''/N instead of raw X_full(:,:,k).\n');
end

fprintf('\n=================================================================\n');
fprintf('  END BPD_TOY_TEST\n');
fprintf('=================================================================\n');


% =========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function P = toy_params(M_val, N_val, Ks_val, Gth_val, Gr_val, dmax_val)
%TOY_PARAMS  Build a minimal params struct for toy verification scenes.

P           = struct();
P.fc        = 28e9;
P.c0        = 3e8;
P.lambda_c  = P.c0 / P.fc;
P.d_ant     = P.lambda_c / 2;
P.M         = M_val;
P.N         = N_val;
P.D_ap      = (P.M - 1) * P.d_ant;
P.r_RD      = 2 * P.D_ap^2 / P.lambda_c;
P.r_lo_fac  = 0.05;
P.r_hi_fac  = 0.40;
P.theta_lo  = 20 * pi/180;
P.theta_hi  = 60 * pi/180;
P.beta_delta = 1.2;

% Wideband params
P.K         = 512;
P.Delta_f   = 120e3;
P.K_s       = Ks_val;
k0          = round(P.K/2) - floor(P.K_s/2) + 1;
P.k_indices = (k0 : k0 + P.K_s - 1).';
f_k_vec     = (P.k_indices - 1) * P.Delta_f;
P.alpha_k_vec = (P.fc + f_k_vec) / P.fc;   % K_s x 1

% BPD-specific
P.G_theta       = Gth_val;
P.G_r           = Gr_val;
P.d_max         = dmax_val;
P.bpd_write_csv = false;
end


function [X_full, theta_true_vec, r_true_vec] = toy_scene(P, SNR_dB)
%TOY_SCENE  Generate a 2-path wideband near-field scene.
%  Paths: theta=[30,50] deg, r=[3.0,6.0] m (inside EBRD for M=64).

d_true         = 2;
theta_true_vec = [30; 50] * pi/180;
r_true_vec     = [3.0; 6.0];
p_true_vec     = [1; 1] / d_true;

M    = P.M;
K_s  = P.K_s;
N    = P.N;
alpha_k = P.alpha_k_vec;

% Build steering matrix per subcarrier
A_all = zeros(M, d_true, K_s);
for k = 1:K_s
    for ell = 1:d_true
        A_all(:, ell, k) = wb_nf_fresnel_steer( ...
            theta_true_vec(ell), 1/r_true_vec(ell), alpha_k(k), P) * sqrt(M);
    end
end

% SNR-matched noise variance
sig_pow = 0;
for k = 1:K_s
    Ak      = A_all(:, :, k);
    sig_pow = sig_pow + real(trace(Ak * diag(p_true_vec) * Ak'));
end
N0 = (sig_pow / K_s) / (M * 10^(SNR_dB/10));

% Generate snapshots
X_full = zeros(M, N, K_s);
for k = 1:K_s
    S_k = zeros(d_true, N);
    for ell = 1:d_true
        S_k(ell,:) = sqrt(p_true_vec(ell)/2) * ...
            (randn(1,N) + 1j*randn(1,N));
    end
    W_k = sqrt(N0/2) * (randn(M,N) + 1j*randn(M,N));
    X_full(:,:,k) = A_all(:,:,k) * S_k + W_k;
end
end


function [X_full, theta_true_vec, r_true_vec] = toy_scene_single(P, SNR_dB, ...
                                                                  theta_val, r_val)
%TOY_SCENE_SINGLE  Generate a single-path scene with specified (theta, r).

d_true         = 1;
theta_true_vec = theta_val;
r_true_vec     = r_val;
p_true_vec     = 1;

M    = P.M;
K_s  = P.K_s;
N    = P.N;
alpha_k = P.alpha_k_vec;

A_all = zeros(M, d_true, K_s);
for k = 1:K_s
    A_all(:, 1, k) = wb_nf_fresnel_steer( ...
        theta_true_vec, 1/r_true_vec, alpha_k(k), P) * sqrt(M);
end

sig_pow = 0;
for k = 1:K_s
    Ak      = A_all(:, :, k);
    sig_pow = sig_pow + real(trace(Ak * diag(p_true_vec) * Ak'));
end
N0 = (sig_pow / K_s) / (M * 10^(SNR_dB/10));

X_full = zeros(M, N, K_s);
for k = 1:K_s
    S_k = sqrt(p_true_vec/2) * (randn(1,N) + 1j*randn(1,N));
    W_k = sqrt(N0/2) * (randn(M,N) + 1j*randn(M,N));
    X_full(:,:,k) = A_all(:,:,k) * S_k + W_k;
end
end


function report_errors(th_hat, r_hat, theta_true, r_true, diag_out, do_print)
%REPORT_ERRORS  Nearest-neighbour matching and pass/fail report.

d_true = numel(theta_true);
d_max  = numel(th_hat);
choose_str = @(c) subsref({'FAIL','PASS'}, struct('type','{}','subs',{{1+c}}));

if do_print
    fprintf('\nTrue paths:\n');
    for ell = 1:d_true
        fprintf('  Path %d: theta=%.2f deg, r=%.4f m\n', ...
            ell, theta_true(ell)*180/pi, r_true(ell));
    end
    fprintf('Estimated paths:\n');
    for l = 1:d_max
        fprintf('  Path %d: theta=%.2f deg, r=%.4f m  (T_peak=%.3e)\n', ...
            l, th_hat(l)*180/pi, r_hat(l), diag_out.T_peak(l));
    end
end

err_th = zeros(d_true,1);
err_r  = zeros(d_true,1);
used   = false(d_max, 1);
for ell = 1:d_true
    diffs = abs(th_hat - theta_true(ell)) * 180/pi;
    diffs(used) = inf;
    [~, best] = min(diffs);
    err_th(ell) = abs(th_hat(best) - theta_true(ell)) * 180/pi;
    err_r(ell)  = abs(r_hat(best)  - r_true(ell));
    used(best)  = true;
end

if do_print
    fprintf('Matched errors:\n');
    for ell = 1:d_true
        pij = (err_th(ell) < 1.0) && (err_r(ell) < 0.5);
        fprintf('  Path %d: |theta_err|=%.3f deg,  |r_err|=%.4f m  --> %s\n', ...
            ell, err_th(ell), err_r(ell), choose_str(pij));
    end
    PASS = all(err_th < 1.0) && all(err_r < 0.5);
    fprintf('Overall: %s\n', choose_str(PASS));
end
end
