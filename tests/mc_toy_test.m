%MC_TOY_TEST  Smoke-test for run_monte_carlo_paperC.m (Task 11.8).
%
%  Runs a single SNR sweep point (SNR=10 dB) with N_MC=5 and verifies:
%    T1 -- All five methods (B1, B2, B4, B5) execute without error; B6 (CRB)
%           also executes without error.
%    T2 -- Output dimensions are correct for all methods.
%    T3 -- mc_toy_diag.csv is written with correct header and N_MC rows;
%           no NaN values in the B1, B2, B4, B5 estimate columns.
%    T4 -- B1 (full-array BPD) at SNR=10 dB, r_hi_fac=0.30 achieves:
%           median |theta_err| * 180/pi < 5.0 deg
%           median |r_err| / r_true     < 0.50
%
%  Acceptance: all four tests must report PASS.
%
%  Author  : R. V. Senyuva (Maltepe University)
%  Date    : May 2026
%  Ref     : Paper C Phase 2, Task 11.8 (Sec. 9 of starter prompt)

clear; close all;
fprintf('=== mc_toy_test: Task 11.8 smoke-test ===\n\n');

% =========================================================================
%  SCENE PARAMETERS (Sec. 9.1 of starter prompt -- exact spec)
% =========================================================================
rng(42, 'twister');

P.c        = 3e8;
P.c0       = P.c;
P.fc       = 28e9;
P.lambda   = P.c / P.fc;
P.lambda_c = P.lambda;
P.M        = 64;
P.d_ant    = P.lambda_c / 2;
P.N_RF     = 8;
P.N        = 64;
P.d        = 1;

D_ap       = (P.M - 1) * P.d_ant;
P.r_RD     = 2 * D_ap^2 / P.lambda_c;
P.r_lo_fac = 0.05;
P.r_hi_fac = 0.30;       % strong near-field (Lesson L33); production uses 1.0
P.u_margin = 2.0;
P.u_min    = 1 / (P.r_hi_fac * P.r_RD * P.u_margin);
P.u_max    = 1 / (P.r_lo_fac * P.r_RD / P.u_margin);
P.theta_lo = 20 * pi/180;
P.theta_hi = 60 * pi/180;
P.Q_theta  = 64;          % reduced for speed; production: 256

P.K_s      = 16;
P.B        = 400e6;
P.Delta_f  = P.B / P.K_s;
P.K        = P.K_s;
k_idx      = (-(P.K_s/2) : (P.K_s/2 - 1)).';
P.alpha_k_vec = 1 + k_idx * P.Delta_f / P.fc;
P.k_indices   = k_idx;

% BPD-specific parameters
P.G_theta       = 64;
P.G_r           = 32;
P.d_max         = P.d;
P.beta_delta    = 1.2;

% WB-CL-KL parameters (defaults matching wb_clkl_driver)
P.lambda_reg    = 1e-4;
P.max_iter      = 80;
P.tol_clkl      = 1e-5;
P.alpha_p       = 0.5;
P.ls_beta       = 0.5;
P.ls_sigma      = 1e-4;

% Tolerances for nf_metrics (used in compute_matched_errors)
P.dtheta_tol    = 5.0 * pi/180;   % 5 deg
P.dr_fac_tol    = 0.50;

SNR_dB     = 10;
N_MC_toy   = 5;
mc_seed_base = 1000;
out_dir    = '.';   % write CSV to current directory

fprintf('Scene: M=%d, N_RF=%d, N=%d, d=%d, K_s=%d, B=%d MHz\n', ...
    P.M, P.N_RF, P.N, P.d, P.K_s, round(P.B/1e6));
fprintf('r_RD = %.2f m   r_hi_fac = %.2f   SNR = %d dB   N_MC = %d\n\n', ...
    P.r_RD, P.r_hi_fac, SNR_dB, N_MC_toy);

% =========================================================================
%  COLLECT PER-REALISATION RESULTS
% =========================================================================
th_true_all = nan(P.d, N_MC_toy);
r_true_all  = nan(P.d, N_MC_toy);
th_B1_all   = nan(P.d, N_MC_toy);
r_B1_all    = nan(P.d, N_MC_toy);
th_B2_all   = nan(P.d, N_MC_toy);
r_B2_all    = nan(P.d, N_MC_toy);
th_B4_all   = nan(P.d, N_MC_toy);
r_B4_all    = nan(P.d, N_MC_toy);
th_B5_all   = nan(P.d, N_MC_toy);
r_B5_all    = nan(P.d, N_MC_toy);

err_B1_flag = false;
err_B2_flag = false;
err_B4_flag = false;
err_B5_flag = false;
err_B6_flag = false;

dim_fail_B1 = false;
dim_fail_B2 = false;
dim_fail_B4 = false;
dim_fail_B5 = false;
dim_fail_crb = false;

% --- CRB: called ONCE outside MC loop (B6) ---
theta_nom = (P.theta_lo + P.theta_hi) / 2;
r_nom     = P.r_RD * (P.r_lo_fac + P.r_hi_fac) / 2;
p_nom     = ones(P.d, 1) / P.d;
N0_nom    = 10^(-SNR_dB/10);

try
    [crb_r, crb_theta, ~] = wb_crb_compressed( ...
        repmat(theta_nom, P.d, 1), repmat(r_nom, P.d, 1), ...
        p_nom, N0_nom, P);
    crb_r_val     = mean(crb_r(:));       % scalar [m]
    crb_theta_val = mean(crb_theta(:));   % scalar [deg]
    fprintf('[B6 CRB] sqrt(CRB_r)=%.4e m   sqrt(CRB_theta)=%.4f deg\n', ...
        crb_r_val, crb_theta_val);
    % Dimension check: d x 1
    if ~isequal(size(crb_r), [P.d, 1]) || ~isequal(size(crb_theta), [P.d, 1])
        dim_fail_crb = true;
        fprintf('  [WARNING] CRB output dimension mismatch.\n');
    end
catch ME
    err_B6_flag = true;
    crb_r_val     = NaN;
    crb_theta_val = NaN;
    fprintf('[B6 CRB] ERROR: %s\n', ME.message);
end

% --- MC loop (serial for toy test) ---
fprintf('\nRunning %d MC realisations...\n', N_MC_toy);

for mc = 1 : N_MC_toy

    rng(mc_seed_base + mc, 'twister');

    % Data generation
    [X_full, Y_full, ~, theta_gen, r_gen, ~, ~, W_comb] = ...
        wb_channel_gen_ofdm_nf(P, SNR_dB);

    th_true_all(:, mc) = theta_gen;
    r_true_all(:, mc)  = r_gen;

    % R_hat_cell formation (locked data-flow, Sec. 6)
    R_hat_cell = cell(P.K_s, 1);
    for k = 1 : P.K_s
        Yk            = Y_full(:, :, k);
        R_hat_cell{k} = (Yk * Yk') / P.N;
        R_hat_cell{k} = (R_hat_cell{k} + R_hat_cell{k}') / 2;
    end

    % B1: WB-BPD
    try
        [th_B1, r_B1, ~, ~] = bpd_baseline(X_full, P);
        th_B1 = th_B1(:);  r_B1 = r_B1(:);
        th_B1_all(:, mc) = th_B1;
        r_B1_all(:, mc)  = r_B1;
        % Dimension check
        if ~isequal(size(th_B1), [P.d, 1]) || ~isequal(size(r_B1), [P.d, 1])
            dim_fail_B1 = true;
        end
    catch ME
        err_B1_flag = true;
        fprintf('  mc=%d B1 error: %s\n', mc, ME.message);
    end

    % B5: WB-DL-OMP
    try
        [th_B5, r_B5, ~] = wb_dl_omp(X_full, P);
        th_B5 = th_B5(:);  r_B5 = r_B5(:);
        th_B5_all(:, mc) = th_B5;
        r_B5_all(:, mc)  = r_B5;
        if ~isequal(size(th_B5), [P.d, 1]) || ~isequal(size(r_B5), [P.d, 1])
            dim_fail_B5 = true;
        end
    catch ME
        err_B5_flag = true;
        fprintf('  mc=%d B5 error: %s\n', mc, ME.message);
    end

    % B2: WB-P-SOMP
    try
        [th_B2, r_B2, ~] = wb_psomp(R_hat_cell, W_comb, P);
        th_B2 = th_B2(:);  r_B2 = r_B2(:);
        th_B2_all(:, mc) = th_B2;
        r_B2_all(:, mc)  = r_B2;
        if ~isequal(size(th_B2), [P.d, 1]) || ~isequal(size(r_B2), [P.d, 1])
            dim_fail_B2 = true;
        end
    catch ME
        err_B2_flag = true;
        fprintf('  mc=%d B2 error: %s\n', mc, ME.message);
    end

    % B4: WB-CL-KL (warm-started by B1)
    if all(isfinite(th_B1_all(:, mc))) && all(isfinite(r_B1_all(:, mc)))
        theta_init = th_B1_all(:, mc);
        r_init     = r_B1_all(:, mc);
    else
        theta_init = repmat((P.theta_lo + P.theta_hi)/2, P.d, 1);
        r_init     = repmat(r_nom, P.d, 1);
    end
    p_init = ones(P.d, 1) / P.d;

    try
        [th_B4, r_B4, ~, ~, ~] = wb_clkl_driver( ...
            R_hat_cell, W_comb, theta_init, r_init, p_init, P);
        th_B4 = th_B4(:);  r_B4 = r_B4(:);
        th_B4_all(:, mc) = th_B4;
        r_B4_all(:, mc)  = r_B4;
        if ~isequal(size(th_B4), [P.d, 1]) || ~isequal(size(r_B4), [P.d, 1])
            dim_fail_B4 = true;
        end
    catch ME
        err_B4_flag = true;
        fprintf('  mc=%d B4 error: %s\n', mc, ME.message);
    end

    fprintf('  mc=%d  theta_true=%.2f deg  r_true=%.3f m\n', ...
        mc, theta_gen(1)*180/pi, r_gen(1));

end  % for mc

fprintf('\n');

% =========================================================================
%  WRITE mc_toy_diag.csv  (30 columns, one row per MC realisation)
% =========================================================================
csv_path = 'mc_toy_diag.csv';
header30 = ['mc_idx,SNR_dB,B_hz,K_s,' ...
            'theta_true_deg,r_true_m,' ...
            'theta_B1_deg,r_B1_m,' ...
            'theta_B2_deg,r_B2_m,' ...
            'theta_B4_deg,r_B4_m,' ...
            'theta_B5_deg,r_B5_m,' ...
            'err_theta_B1_deg,err_r_B1_m,' ...
            'err_theta_B2_deg,err_r_B2_m,' ...
            'err_theta_B4_deg,err_r_B4_m,' ...
            'err_theta_B5_deg,err_r_B5_m,' ...
            'crb_r_m,crb_theta_deg,' ...
            'N_RF,M,N,d,r_hi_fac,Q_theta'];

fid = fopen(csv_path, 'w');
if fid == -1
    error('mc_toy_test: cannot create %s', csv_path);
end
fprintf(fid, '%s\n', header30);

for mc = 1 : N_MC_toy
    th_t  = th_true_all(1, mc) * 180/pi;
    r_t   = r_true_all(1, mc);

    th_b1 = th_B1_all(1, mc) * 180/pi;
    r_b1  = r_B1_all(1, mc);
    th_b2 = th_B2_all(1, mc) * 180/pi;
    r_b2  = r_B2_all(1, mc);
    th_b4 = th_B4_all(1, mc) * 180/pi;
    r_b4  = r_B4_all(1, mc);
    th_b5 = th_B5_all(1, mc) * 180/pi;
    r_b5  = r_B5_all(1, mc);

    e_th_b1 = (th_B1_all(1,mc) - th_true_all(1,mc)) * 180/pi;
    e_r_b1  = r_B1_all(1,mc) - r_true_all(1,mc);
    e_th_b2 = (th_B2_all(1,mc) - th_true_all(1,mc)) * 180/pi;
    e_r_b2  = r_B2_all(1,mc) - r_true_all(1,mc);
    e_th_b4 = (th_B4_all(1,mc) - th_true_all(1,mc)) * 180/pi;
    e_r_b4  = r_B4_all(1,mc) - r_true_all(1,mc);
    e_th_b5 = (th_B5_all(1,mc) - th_true_all(1,mc)) * 180/pi;
    e_r_b5  = r_B5_all(1,mc) - r_true_all(1,mc);

    fprintf(fid, '%d,%.2f,%.6g,%d,', mc, SNR_dB, P.B, P.K_s);
    fprintf(fid, '%.6g,%.6g,', th_t, r_t);
    fprintf(fid, '%.6g,%.6g,', th_b1, r_b1);
    fprintf(fid, '%.6g,%.6g,', th_b2, r_b2);
    fprintf(fid, '%.6g,%.6g,', th_b4, r_b4);
    fprintf(fid, '%.6g,%.6g,', th_b5, r_b5);
    fprintf(fid, '%.6g,%.6g,', e_th_b1, e_r_b1);
    fprintf(fid, '%.6g,%.6g,', e_th_b2, e_r_b2);
    fprintf(fid, '%.6g,%.6g,', e_th_b4, e_r_b4);
    fprintf(fid, '%.6g,%.6g,', e_th_b5, e_r_b5);
    fprintf(fid, '%.6g,%.6g,', crb_r_val, crb_theta_val);
    fprintf(fid, '%d,%d,%d,%d,%.4f,%d\n', ...
        P.N_RF, P.M, P.N, P.d, P.r_hi_fac, P.Q_theta);
end
fclose(fid);
fprintf('[CSV] Wrote %s (%d rows, 30 columns)\n\n', csv_path, N_MC_toy);

% =========================================================================
%  T1 -- ALL METHODS EXECUTE WITHOUT ERROR
% =========================================================================
fprintf('--- T1: Method execution ---\n');
T1_pass = true;
check_flag('B1 WB-BPD',    ~err_B1_flag);
check_flag('B2 WB-P-SOMP', ~err_B2_flag);
check_flag('B4 WB-CL-KL',  ~err_B4_flag);
check_flag('B5 WB-DL-OMP', ~err_B5_flag);
check_flag('B6 CRB',       ~err_B6_flag);
T1_pass = ~(err_B1_flag || err_B2_flag || err_B4_flag || err_B5_flag || err_B6_flag);
print_result('T1', T1_pass);

% =========================================================================
%  T2 -- OUTPUT DIMENSIONS CORRECT
% =========================================================================
fprintf('\n--- T2: Output dimensions ---\n');
check_flag('B1 theta/r: [d x 1]', ~dim_fail_B1);
check_flag('B2 theta/r: [d x 1]', ~dim_fail_B2);
check_flag('B4 theta/r: [d x 1]', ~dim_fail_B4);
check_flag('B5 theta/r: [d x 1]', ~dim_fail_B5);
check_flag('B6 crb_r:   [d x 1]', ~dim_fail_crb);
T2_pass = ~(dim_fail_B1 || dim_fail_B2 || dim_fail_B4 || dim_fail_B5 || dim_fail_crb);
print_result('T2', T2_pass);

% =========================================================================
%  T3 -- CSV WRITTEN CORRECTLY
% =========================================================================
fprintf('\n--- T3: CSV structure ---\n');
T3_pass = false;

csv_ok = isfile(csv_path);
fprintf('  CSV exists:          %s\n', yn(csv_ok));
if csv_ok
    % Check header
    fid2 = fopen(csv_path, 'r');
    hdr  = fgetl(fid2);
    expected_cols = 30;
    actual_cols   = numel(strsplit(hdr, ','));
    fprintf('  Header column count: %d (expected %d)  %s\n', ...
        actual_cols, expected_cols, yn(actual_cols == expected_cols));

    % Count data rows
    n_rows = 0;
    while ~feof(fid2)
        line = fgetl(fid2);
        if ischar(line) && ~isempty(strtrim(line))
            n_rows = n_rows + 1;
        end
    end
    fclose(fid2);
    fprintf('  Data rows:           %d (expected %d)  %s\n', ...
        n_rows, N_MC_toy, yn(n_rows == N_MC_toy));

    % Check for NaN in estimate columns (B1, B2, B4, B5)
    % Columns: theta_B1_deg (7), r_B1_m (8), theta_B2_deg (9), r_B2_m (10),
    %          theta_B4_deg (11), r_B4_m (12), theta_B5_deg (13), r_B5_m (14)
    data_mat = readmatrix(csv_path, 'NumHeaderLines', 1);
    est_cols = data_mat(:, 7:14);
    has_nan  = any(isnan(est_cols(:)));
    fprintf('  No NaN in B1-B5 est cols: %s\n', yn(~has_nan));

    T3_pass = csv_ok && (actual_cols == expected_cols) && ...
              (n_rows == N_MC_toy) && ~has_nan;
end
print_result('T3', T3_pass);

% =========================================================================
%  T4 -- B1 SANITY CHECK (BPD accuracy at SNR=10 dB, r_hi_fac=0.30)
% =========================================================================
fprintf('\n--- T4: B1 accuracy sanity (BPD at SNR=10 dB, r_hi_fac=0.30) ---\n');
T4_pass = false;

valid_mask = ~isnan(th_B1_all(1,:)) & ~isnan(r_B1_all(1,:));
if sum(valid_mask) < 1
    fprintf('  No valid B1 realisations -- T4 FAIL.\n');
else
    abs_th_err_deg = abs(th_B1_all(1, valid_mask) - th_true_all(1, valid_mask)) * 180/pi;
    rel_r_err      = abs(r_B1_all(1, valid_mask) - r_true_all(1, valid_mask)) ...
                     ./ abs(r_true_all(1, valid_mask));

    med_th = median(abs_th_err_deg);
    med_r  = median(rel_r_err);

    thr_th = 5.0;    % deg
    thr_r  = 0.50;   % relative

    fprintf('  Median |theta_err| = %.3f deg  (threshold: < %.1f deg)  %s\n', ...
        med_th, thr_th, yn(med_th < thr_th));
    fprintf('  Median |r_err|/r   = %.4f     (threshold: < %.2f)       %s\n', ...
        med_r, thr_r, yn(med_r < thr_r));

    T4_pass = (med_th < thr_th) && (med_r < thr_r);
end
print_result('T4', T4_pass);

% =========================================================================
%  FINAL VERDICT
% =========================================================================
fprintf('\n=========================================\n');
all_pass = T1_pass && T2_pass && T3_pass && T4_pass;
if all_pass
    fprintf('  RESULT: ALL TESTS PASSED -- Phase 2 gate satisfied.\n');
else
    fprintf('  RESULT: ONE OR MORE TESTS FAILED.\n');
    fprintf('  T1=%s  T2=%s  T3=%s  T4=%s\n', ...
        pf(T1_pass), pf(T2_pass), pf(T3_pass), pf(T4_pass));
end
fprintf('=========================================\n');

% =========================================================================
%  LOCAL HELPERS
% =========================================================================
function check_flag(label, ok)
fprintf('  %-30s %s\n', label, yn(ok));
end

function s = yn(ok)
if ok; s = 'OK'; else; s = 'FAIL'; end
end

function s = pf(ok)
if ok; s = 'PASS'; else; s = 'FAIL'; end
end

function print_result(tag, ok)
if ok
    fprintf('  %s: PASS\n', tag);
else
    fprintf('  %s: FAIL\n', tag);
end
end
