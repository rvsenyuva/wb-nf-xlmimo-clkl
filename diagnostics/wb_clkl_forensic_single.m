function wb_clkl_forensic_single
%WB_CLKL_FORENSIC_SINGLE  Single-realisation forensic for Paper C Task 11.3 rev3.
%
%  Runs two target realisations at SNR=10 dB in convergence mode:
%    mc=1   -- typical hit-cap realisation (Phase B drifts; H1/H3 trajectory test).
%    mc=30  -- fast-converge realisation (Phase B stationary at BPD; H2/H3 test).
%
%  For each target this script:
%    1. Generates observations identically to run_one_realisation
%       (rng(1000+mc,'twister') -> wb_channel_gen_ofdm_nf_fixed).
%    2. Runs BPD warm-start.
%    3. Reproduces v2 driver Phase A initialisation (eta_init).
%    4. Reproduces v2 driver Phase B inner loop step-by-step,
%       logging (iter, r, theta, L, ||grad_kappa||, ||grad_omega||, alpha)
%       at every iteration up to iter=50.
%    5. Reproduces v2 driver Phase C 4-pass alternating MF scan
%       (seeded from BPD per D5; for d=1, deflation is trivial).
%    6. Reproduces Phase D KL arbitration; additionally evaluates
%       L_C = KL(theta_BPD, u_BPD) -- the BPD-anchor candidate that
%       current Phase D does NOT consider (this is the H3 test).
%
%  Output files:
%    forensic_single_mc1_traj.csv   -- iter-by-iter Phase B trajectory mc=1
%    forensic_single_mc30_traj.csv  -- iter-by-iter Phase B trajectory mc=30
%    forensic_single_diag.csv       -- summary row per mc target
%
%  Discrimination (read after running):
%    H1 confirmed if r in trajectory drifts monotonically away from r_true.
%    H3 confirmed if L_C < min(L_A, L_B), especially in mc=30 where
%      Phase B is stationary -- meaning corruption is purely in Phase D.
%
%  Author : R. V. Senyuva (Maltepe University)
%  Date   : May 2026
%  Ref    : Paper C Phase 3, Task 11.3 (rev3) Forensic Diagnostic
%  Reqs   : setup_production_P_v4.m, wb_clkl_estimator.m (v2-canonical),
%           wb_nf_fresnel_steer.m, bpd_baseline.m

% =========================================================================
%  CONFIGURATION
% =========================================================================
mc_targets = [1, 30];
SNR_dB     = 10;

% =========================================================================
%  PARAMETER SETUP (identical to production convergence sweep)
% =========================================================================
P                  = setup_production_P_v4('convergence', 'restricted');
P.use_multi_start  = false;   % single-start for clean Phase B trajectory

d         = P.d;
M         = P.M;
N_RF      = P.N_RF;
K_s       = P.K_s;
c_lin     = 2*pi*P.d_ant   / P.lambda_c;
c_quad    =     pi*P.d_ant^2 / P.lambda_c;
n_log_iter = 50;   % first 50 Phase B iterations logged in trajectory CSV
abl_ponly  = false;

fprintf('\n================================================================\n');
fprintf(' Paper C -- Task 11.3 rev3 Forensic Diagnostic (single-realisation)\n');
fprintf(' SNR = %d dB,  mc targets = [%d, %d]\n', SNR_dB, mc_targets(1), mc_targets(2));
fprintf(' Fixed scene: theta_true = %.2f deg, r_true = %.2f m\n', ...
        P.conv_theta*180/pi, P.conv_r);
fprintf('================================================================\n');

% Pre-allocate summary storage (one row per mc target)
n_mc = numel(mc_targets);
summary = nan(n_mc, 22);   % columns documented at CSV write below

% =========================================================================
%  PER-TARGET FORENSIC LOOP
% =========================================================================
for it = 1:n_mc
    mc = mc_targets(it);

    fprintf('\n----------------------------------------------------------------\n');
    fprintf(' mc_idx = %d   (rng seed = %d)\n', mc, 1000+mc);
    fprintf('----------------------------------------------------------------\n');

    % --- Generate observations identically to run_one_realisation ---
    rng(1000+mc, 'twister');
    P_fixed                       = P;
    P_fixed.convergence_fixed_scene = false;   % prevent recursion guard
    [X_full, Y_full, ~, theta_true_v, r_true_v, ~, ~, W_comb] = ...
        wb_channel_gen_ofdm_nf_fixed(P_fixed, SNR_dB, P.conv_theta, P.conv_r);

    theta_true_deg = theta_true_v(1) * 180/pi;
    r_true_m       = r_true_v(1);

    % --- Compressed sample covariance per subcarrier ---
    R_hat_cell = cell(K_s, 1);
    for k = 1:K_s
        Yk             = Y_full(:, :, k);
        R_hat_cell{k}  = (Yk * Yk') / P.N;
        R_hat_cell{k}  = (R_hat_cell{k} + R_hat_cell{k}') / 2;
    end

    % --- BPD warm-start ---
    [th_bpd, r_bpd_vec, ~, ~] = bpd_baseline(X_full, P);
    theta_bpd     = th_bpd(1);
    r_bpd         = r_bpd_vec(1);
    theta_bpd_deg = theta_bpd * 180/pi;

    fprintf('  r_true     = %10.4f m   theta_true     = %8.4f deg\n', ...
            r_true_m, theta_true_deg);
    fprintf('  r_bpd      = %10.4f m   theta_bpd      = %8.4f deg   (err: %+.4f m, %+.4f deg)\n', ...
            r_bpd, theta_bpd_deg, r_bpd - r_true_m, theta_bpd_deg - theta_true_deg);

    % =====================================================================
    %  PHASE A -- INITIALISATION (mirror wb_clkl_driver Phase A)
    % =====================================================================
    R_bar = zeros(N_RF, N_RF);
    for k = 1:K_s, R_bar = R_bar + R_hat_cell{k}; end
    R_bar      = (R_bar + R_bar') / (2*K_s);
    ev_sorted  = sort(real(eig(R_bar)), 'ascend');
    n_noise_ev = max(1, N_RF - d);
    N0_init    = max(mean(ev_sorted(1:n_noise_ev)), 1e-12);

    p_init_val = ones(d, 1) / d;
    omega_init = c_lin  * cos(theta_bpd);
    u_init     = min(P.u_max, max(P.u_min, 1./r_bpd));
    kappa_init = c_quad * sin(theta_bpd).^2 .* u_init;

    eta = [omega_init; kappa_init; p_init_val; N0_init];

    % Initial estimator call -> warm-start KL and preconditioner
    [L_curr, grad_L, ~] = wb_clkl_estimator(eta, R_hat_cell, W_comb, P);
    L_warm   = L_curr;
    eps_prec = 1e-6;
    D_prec   = build_preconditioner(grad_L, d, eps_prec, abl_ponly);
    alpha0   = P.alpha_p / N_RF;

    % =====================================================================
    %  PHASE B -- BLOCK-PRECONDITIONED ARMIJO (mirror v2 driver lines ~201-272)
    % =====================================================================
    % Trajectory log buffer (only first n_log_iter+1 iterations)
    traj = nan(n_log_iter+1, 7);   % cols: iter, r, theta_deg, L, ||grad_k||, ||grad_o||, alpha

    [theta_c, u_c]   = eta_to_physical(eta, d, c_lin, c_quad, P.u_min, P.u_max);
    traj(1, :) = [0, 1./u_c, theta_c*180/pi, L_curr, ...
                  norm(grad_L(d+1:2*d)), norm(grad_L(1:d)), NaN];

    fprintf('\n  Phase B trajectory (first %d iters logged; max_iter=%d):\n', ...
            n_log_iter, P.max_iter);
    fprintf('   iter | r_curr (m) | theta (deg) |    L_curr    |  ||grad_k||  |  ||grad_o||  |   alpha\n');
    fprintf('   -----|------------|-------------|--------------|--------------|--------------|----------\n');
    fprintf('     0  | %10.4f | %11.4f |  %10.4f  |  %.3e  |  %.3e  |    -\n', ...
            1./u_c, theta_c*180/pi, L_curr, ...
            norm(grad_L(d+1:2*d)), norm(grad_L(1:d)));

    converged = false;
    n_iter    = 0;
    min_iter  = 3;
    consec    = 0;

    for t = 1:P.max_iter

        if abl_ponly
            grad_step             = zeros(3*d+1, 1);
            grad_step(2*d+1:3*d)  = grad_L(2*d+1:3*d);
        else
            grad_step = grad_L;
        end

        dir      = D_prec .* grad_step;
        suff_dec = grad_step' * dir;

        if suff_dec < 1e-15
            converged = true;
            n_iter    = t;
            break;
        end

        alpha  = alpha0;
        L_prev = L_curr;
        for ls = 1:25
            eta_try                = eta - alpha * dir;
            eta_try(2*d+1:3*d)     = max(0,     eta_try(2*d+1:3*d));
            eta_try(3*d+1)         = max(1e-12, eta_try(3*d+1));
            [L_try, ~, ~]          = wb_clkl_estimator(eta_try, R_hat_cell, W_comb, P);
            if L_try <= L_prev - P.ls_sigma * alpha * suff_dec
                break;
            end
            alpha = alpha * P.ls_beta;
        end

        eta    = eta_try;
        L_prev = L_curr;
        [L_curr, grad_L, ~] = wb_clkl_estimator(eta, R_hat_cell, W_comb, P);
        n_iter = t;

        % Log first n_log_iter iterations
        if t <= n_log_iter
            [theta_c, u_c] = eta_to_physical(eta, d, c_lin, c_quad, P.u_min, P.u_max);
            traj(t+1, :) = [t, 1./u_c, theta_c*180/pi, L_curr, ...
                            norm(grad_L(d+1:2*d)), norm(grad_L(1:d)), alpha];
            fprintf('   %3d  | %10.4f | %11.4f |  %10.4f  |  %.3e  |  %.3e  |  %.2e\n', ...
                    t, 1./u_c, theta_c*180/pi, L_curr, ...
                    norm(grad_L(d+1:2*d)), norm(grad_L(1:d)), alpha);
        end

        rel_change = abs(L_curr - L_prev) / (abs(L_prev) + 1e-15);
        if rel_change < P.tol_clkl && t > 5
            consec = consec + 1;
        else
            consec = 0;
        end
        if consec >= min_iter
            converged = true;
            break;
        end
    end

    if n_iter == n_log_iter && ~converged
        fprintf('   ... (Phase B continues; trajectory log truncated at iter %d)\n', n_log_iter);
    end

    % --- Phase B end state ---
    [theta_PhaseB, u_PhaseB] = eta_to_physical(eta, d, c_lin, c_quad, P.u_min, P.u_max);
    r_PhaseB     = 1 ./ u_PhaseB;
    L_PhaseB     = L_curr;
    p_ref        = max(0, eta(2*d+1:3*d));
    N0_converged = max(1e-12, eta(3*d+1));

    fprintf('\n  Phase B end (n_iter=%d, converged=%d):\n', n_iter, converged);
    fprintf('    r_PhaseB     = %10.4f m   theta_PhaseB     = %8.4f deg\n', ...
            r_PhaseB, theta_PhaseB*180/pi);
    fprintf('    drift from BPD: dr = %+.4f m,  dtheta = %+.4f deg\n', ...
            r_PhaseB - r_bpd, (theta_PhaseB - theta_bpd)*180/pi);
    fprintf('    delta_L: L_warm=%.4f -> L_PhaseB=%.4f  (Delta=%.4f)\n', ...
            L_warm, L_PhaseB, L_warm - L_PhaseB);

    % --- Save trajectory CSV ---
    traj_name = sprintf('forensic_single_mc%d_traj.csv', mc);
    fid = fopen(traj_name, 'w');
    fprintf(fid, 'iter,r_curr_m,theta_curr_deg,L_curr,grad_kappa_norm,grad_omega_norm,alpha_step\n');
    for i = 1:size(traj, 1)
        if isnan(traj(i, 1)) && i > 1, break; end
        if isnan(traj(i, 7))
            fprintf(fid, '%d,%.6f,%.6f,%.6f,%.6e,%.6e,\n', ...
                    traj(i,1), traj(i,2), traj(i,3), traj(i,4), traj(i,5), traj(i,6));
        else
            fprintf(fid, '%d,%.6f,%.6f,%.6f,%.6e,%.6e,%.6e\n', traj(i, :));
        end
    end
    fclose(fid);
    fprintf('  Wrote: %s\n', traj_name);

    % =====================================================================
    %  PHASE C -- 4-PASS ALTERNATING SCAN (mirror driver Phase C)
    %  Seeded from BPD per D5; for d=1, deflation is a no-op.
    % =====================================================================
    th_scan  = linspace(P.theta_lo, P.theta_hi, P.Q_theta);
    Q_scan_u = P.G_r;    u_scan   = linspace(P.u_min, P.u_max, Q_scan_u);

    theta_ref     = theta_bpd;        % seeded from BPD per D5
    u_ref         = u_init;           % seeded from BPD
    u_ref_pass    = nan(1, 4);

    fprintf('\n  Phase C -- 4-pass alternating MF scan (BPD-seeded):\n');

    for pass = 1:4
        if mod(pass, 2) == 1
            % Theta scan, u fixed
            sc = zeros(1, P.Q_theta);
            for k = 1:K_s
                Ds = scan_atoms_theta(th_scan, u_ref, P.alpha_k_vec(k), W_comb, M, P);
                sc = sc + real(sum(conj(Ds) .* (R_hat_cell{k} * Ds), 1));
            end
            [~, bi]   = max(sc);
            theta_ref = th_scan(bi);
            fprintf('    Pass %d (theta-scan, u_ref=%.4f, r=%.4f m): theta_ref = %8.4f deg\n', ...
                    pass, u_ref, 1./u_ref, theta_ref*180/pi);
        else
            % u scan, theta fixed
            sc = zeros(1, Q_scan_u);
            for k = 1:K_s
                Ds = scan_atoms_u(u_scan, theta_ref, P.alpha_k_vec(k), W_comb, M, P);
                sc = sc + real(sum(conj(Ds) .* (R_hat_cell{k} * Ds), 1));
            end
            [~, bi] = max(sc);
            u_ref   = u_scan(bi);
            fprintf('    Pass %d (u-scan, theta_ref=%.4f deg): u_ref = %.4f -> r = %8.4f m\n', ...
                    pass, theta_ref*180/pi, u_ref, 1./u_ref);

            % Pass 2 -- flatness diagnostic at u_true and u_bpd
            if pass == 2
                u_true_val = 1 ./ r_true_m;
                idx_max    = bi;
                idx_true   = max(1, min(Q_scan_u, ...
                    round((u_true_val - P.u_min) / (P.u_max - P.u_min) * (Q_scan_u-1)) + 1));
                idx_bpd    = max(1, min(Q_scan_u, ...
                    round((u_init - P.u_min) / (P.u_max - P.u_min) * (Q_scan_u-1)) + 1));
                sc_max  = sc(idx_max);
                sc_true = sc(idx_true);
                sc_bpd  = sc(idx_bpd);
                fprintf('      u-scan flatness diagnostic:\n');
                fprintf('        argmax: idx=%4d  u=%.4f  r=%8.4f m  score=%12.4f  (peak)\n', ...
                        idx_max,  u_scan(idx_max),  1./u_scan(idx_max),  sc_max);
                fprintf('        u_true: idx=%4d  u=%.4f  r=%8.4f m  score=%12.4f  (%.4f x peak)\n', ...
                        idx_true, u_scan(idx_true), 1./u_scan(idx_true), sc_true, sc_true/sc_max);
                fprintf('        u_bpd : idx=%4d  u=%.4f  r=%8.4f m  score=%12.4f  (%.4f x peak)\n', ...
                        idx_bpd,  u_scan(idx_bpd),  1./u_scan(idx_bpd),  sc_bpd,  sc_bpd /sc_max);
            end
        end
        u_ref_pass(pass) = u_ref;
    end

    % =====================================================================
    %  PHASE D -- KL ARBITRATION (with diagnostic Cand C at BPD anchor)
    % =====================================================================
    fprintf('\n  Phase D -- KL arbitration:\n');

    % Cand A: (theta_ref, u_ref from Phase C scan)
    omega_A    = c_lin  * cos(theta_ref);
    kappa_A    = c_quad * sin(theta_ref).^2 .* u_ref;
    eta_A      = [omega_A; kappa_A; p_ref; N0_converged];
    [L_A, ~, ~] = wb_clkl_estimator(eta_A, R_hat_cell, W_comb, P);

    % Cand B: (theta_ref, u_init from BPD warm-start)
    kappa_B    = c_quad * sin(theta_ref).^2 .* u_init;
    eta_B      = [omega_A; kappa_B; p_ref; N0_converged];
    [L_B, ~, ~] = wb_clkl_estimator(eta_B, R_hat_cell, W_comb, P);

    % Cand C: (theta_BPD, u_BPD) -- BPD-anchor (NOT in current Phase D logic)
    omega_C    = c_lin  * cos(theta_bpd);
    kappa_C    = c_quad * sin(theta_bpd).^2 .* u_init;
    eta_C      = [omega_C; kappa_C; p_ref; N0_converged];
    [L_C, ~, ~] = wb_clkl_estimator(eta_C, R_hat_cell, W_comb, P);

    fprintf('    L_A (theta_ref, u_ref) = %12.4f   r=%.4f m, theta=%.4f deg\n', ...
            L_A, 1./u_ref, theta_ref*180/pi);
    fprintf('    L_B (theta_ref, u_bpd) = %12.4f   r=%.4f m, theta=%.4f deg\n', ...
            L_B, 1./u_init, theta_ref*180/pi);
    fprintf('    L_C (theta_bpd, u_bpd) = %12.4f   r=%.4f m, theta=%.4f deg   [DIAGNOSTIC: not in current Phase D]\n', ...
            L_C, 1./u_init, theta_bpd*180/pi);

    % Current Phase D selects min(L_A, L_B)
    if L_A <= L_B
        sel_label   = 'A';
        r_final     = 1 ./ u_ref;
        theta_final = theta_ref;
    else
        sel_label   = 'B';
        r_final     = 1 ./ u_init;
        theta_final = theta_ref;
    end
    fprintf('    -> Phase D (current) selects Cand %s: r_hat = %.4f m, theta_hat = %.4f deg\n', ...
            sel_label, r_final, theta_final*180/pi);

    % H3 diagnostic
    if L_C < min(L_A, L_B)
        fprintf('    *** H3 CONFIRMED: L_C < min(L_A, L_B). BPD anchor would be best. ***\n');
        fprintf('        L_A - L_C = %.4f,  L_B - L_C = %.4f\n', L_A - L_C, L_B - L_C);
        fprintf('        Adding Cand C to Phase D would yield r_hat = %.4f m (truth %.4f m).\n', ...
                1./u_init, r_true_m);
    elseif L_C >= min(L_A, L_B)
        fprintf('    H3 NOT confirmed at this realisation: L_C is not the minimum.\n');
        fprintf('        Failure mechanism is more subtle than D-OUT-2 alone.\n');
    end

    % --- Populate summary row ---
    summary(it, :) = [mc, SNR_dB, r_bpd, theta_bpd_deg, ...
                      r_PhaseB, theta_PhaseB*180/pi, ...
                      r_final, theta_final*180/pi, ...
                      r_true_m, theta_true_deg, ...
                      L_warm, L_PhaseB, L_A, L_B, L_C, ...
                      double(L_A <= L_B), ...
                      double(abs(r_final - r_bpd) < 0.5), ...
                      u_ref_pass(1), u_ref_pass(2), u_ref_pass(3), u_ref_pass(4), ...
                      n_iter];
end  % mc target loop

% =========================================================================
%  WRITE SUMMARY CSV
% =========================================================================
sum_name = 'forensic_single_diag.csv';
fid = fopen(sum_name, 'w');
fprintf(fid, ['mc_idx,SNR_dB,r_bpd_m,theta_bpd_deg,', ...
              'r_PhaseB_m,theta_PhaseB_deg,r_final_m,theta_final_deg,', ...
              'r_true_m,theta_true_deg,', ...
              'L_warm,L_PhaseB,L_A,L_B,L_C,', ...
              'phase_D_select_AleqB,proxy_CandB_picked,', ...
              'u_ref_pass1,u_ref_pass2,u_ref_pass3,u_ref_pass4,', ...
              'n_iter\n']);
for i = 1:n_mc
    fprintf(fid, '%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%.6f,%.6f,%.6f,%.6f,%d\n', ...
            summary(i, 1), summary(i, 2), summary(i, 3), summary(i, 4), ...
            summary(i, 5), summary(i, 6), summary(i, 7), summary(i, 8), ...
            summary(i, 9), summary(i, 10), summary(i, 11), summary(i, 12), ...
            summary(i, 13), summary(i, 14), summary(i, 15), summary(i, 16), ...
            summary(i, 17), summary(i, 18), summary(i, 19), summary(i, 20), ...
            summary(i, 21), summary(i, 22));
end
fclose(fid);

fprintf('\n================================================================\n');
fprintf(' Forensic complete. Summary: %s\n', sum_name);
fprintf(' Per-trajectory CSVs:\n');
for i = 1:n_mc
    fprintf('   forensic_single_mc%d_traj.csv\n', mc_targets(i));
end
fprintf('================================================================\n');

end  % main function


% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function D = build_preconditioner(grad_L, d, eps_prec, abl_ponly)
%BUILD_PRECONDITIONER  Diagonal scaling for Phase B preconditioned descent.
%  Mirrors v2 driver's build_preconditioner local function.
g_omega = abs(grad_L(1:d))         + eps_prec;
g_kappa = abs(grad_L(d+1:2*d))     + eps_prec;
g_p     = abs(grad_L(2*d+1:3*d))   + eps_prec;
g_N0    = abs(grad_L(3*d+1))       + eps_prec;

D                  = zeros(3*d+1, 1);
D(1:d)             = 1 ./ g_omega;
D(d+1:2*d)         = 1 ./ g_kappa;
D(2*d+1:3*d)       = 1 ./ g_p;
D(3*d+1)           = 1 ./ g_N0;

if abl_ponly
    D(1:2*d)       = 0;
    D(3*d+1)       = 0;
end
end

function [theta, u] = eta_to_physical(eta, d, c_lin, c_quad, u_min, u_max)
%ETA_TO_PHYSICAL  Decompose eta = [omega; kappa; p; N0] into (theta, u).
%  Mirrors v2 driver's eta_to_physical local function.
omega = eta(1:d);
kappa = eta(d+1:2*d);

% Recover theta from omega via cos^{-1}(omega/c_lin); clamp argument
arg_c   = max(-1, min(1, omega ./ c_lin));
theta   = acos(arg_c);

% Recover u from kappa via kappa = c_quad * sin(theta)^2 * u
sin2    = sin(theta).^2;
sin2    = max(sin2, 1e-12);
u_raw   = kappa ./ (c_quad .* sin2);
u       = min(u_max, max(u_min, u_raw));
end

function Ds = scan_atoms_theta(th_scan, u_fix, alpha_k, W_comb, M, P)
%SCAN_ATOMS_THETA  Build compressed steering atoms across theta_scan grid.
%  u fixed; theta varies. Returns N_RF x Q_theta complex matrix.
m_bar  = ((0:M-1) - (M-1)/2).';
Q      = numel(th_scan);
Ds     = zeros(P.N_RF, Q);
for q = 1:Q
    th_q   = th_scan(q);
    omega  = (2*pi*P.d_ant/P.lambda_c) * cos(th_q);
    kappa  = (pi*P.d_ant^2/P.lambda_c) * sin(th_q)^2 * u_fix;
    a_th   = exp(1j*alpha_k*(omega*m_bar - kappa*m_bar.^2)) / sqrt(M);
    Ds(:,q) = W_comb' * a_th;
end
% Normalise atoms to unit norm for fair score across grid
nrm    = sqrt(sum(abs(Ds).^2, 1));
nrm(nrm < 1e-15) = 1;
Ds     = Ds ./ nrm;
end

function Ds = scan_atoms_u(u_scan, theta_fix, alpha_k, W_comb, M, P)
%SCAN_ATOMS_U  Build compressed steering atoms across u_scan grid.
%  theta fixed; u varies. Returns N_RF x Q_u complex matrix.
m_bar  = ((0:M-1) - (M-1)/2).';
Q      = numel(u_scan);
Ds     = zeros(P.N_RF, Q);
omega  = (2*pi*P.d_ant/P.lambda_c) * cos(theta_fix);
sin2   = sin(theta_fix)^2;
for q = 1:Q
    kappa  = (pi*P.d_ant^2/P.lambda_c) * sin2 * u_scan(q);
    a_u    = exp(1j*alpha_k*(omega*m_bar - kappa*m_bar.^2)) / sqrt(M);
    Ds(:,q) = W_comb' * a_u;
end
nrm    = sqrt(sum(abs(Ds).^2, 1));
nrm(nrm < 1e-15) = 1;
Ds     = Ds ./ nrm;
end

function [X_full, Y_full, H_true, theta_true, r_true, p_true, N0, W_comb] = ...
        wb_channel_gen_ofdm_nf_fixed(P, SNR_dB, theta_fix, r_fix)
%WB_CHANNEL_GEN_OFDM_NF_FIXED  Channel generator with pinned geometry.
%  VERBATIM copy from run_monte_carlo_paperC.m local function (per B3
%  confirmation). Local copy makes forensic self-contained.
P_tmp           = P;
P_tmp.r_lo_fac  = max(P.r_lo_fac, r_fix / P.r_RD * 0.99);
P_tmp.r_hi_fac  = min(P.r_hi_fac, r_fix / P.r_RD * 1.01);

[~, ~, ~, ~, ~, p_true, N0, W_comb] = wb_channel_gen_ofdm_nf(P_tmp, SNR_dB);

theta_true = repmat(theta_fix, P.d, 1);
r_true     = repmat(r_fix,     P.d, 1);
p_true     = ones(P.d, 1) / P.d;

M    = P.M;
N    = P.N;
K_s  = P.K_s;
m_bar = ((0:M-1) - (M-1)/2).';

X_full = zeros(M, N, K_s);
H_true = zeros(M, P.d, K_s);

for k = 1:K_s
    alpha_k = P.alpha_k_vec(k);
    A_k     = zeros(M, P.d);
    for l = 1:P.d
        omega_l   = (2*pi*P.d_ant/P.lambda_c) * cos(theta_true(l));
        kappa_l   = (pi*P.d_ant^2/P.lambda_c) * sin(theta_true(l))^2 / r_true(l);
        A_k(:,l)  = exp(1j*alpha_k*(omega_l*m_bar - kappa_l*m_bar.^2)) * sqrt(M);
    end
    H_true(:,:,k) = A_k;
    S_k = sqrt(p_true(1)) * (randn(P.d, N) + 1j*randn(P.d, N)) / sqrt(2);
    W_k = sqrt(N0)        * (randn(M,   N) + 1j*randn(M,   N)) / sqrt(2);
    X_full(:,:,k) = A_k * S_k + W_k;
end

Y_full = zeros(P.N_RF, N, K_s);
for k = 1:K_s
    Y_full(:,:,k) = W_comb' * X_full(:,:,k);
end
end
