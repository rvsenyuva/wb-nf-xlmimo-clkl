function [theta_hat, r_hat, p_hat, N0_hat, info] = ...
    wb_clkl_driver(R_hat_cell, W_comb, theta_init, r_init, p_init, P)
%WB_CLKL_DRIVER  4-phase WB-CL-KL optimisation driver (v2, Component 2).
%
%  Paper C Phase 3, Task 11.3 (rev3) Component 2.  Orchestrates the full
%  WB-CL-KL estimation pipeline by wrapping wb_clkl_estimator.m (Task 11.3)
%  in a block-preconditioned Armijo projected-gradient-descent loop with
%  multi-start BPD warm-start and post-loop 4-pass residual matched-filter
%  refinement.
%
%  Called by wb_clkl_driver_toy_test_v2.m (Task 11.4 toy test, v2) and by
%  run_monte_carlo_paperC.m (Task 11.8).  The driver is the ONLY entry
%  point into the WB-CL-KL pipeline; wb_clkl_estimator.m is NEVER called
%  directly by Monte Carlo.
%
%  v2 CHANGES (Component 2)
%  ------------------------
%  - Multi-start warm-start replaces single-start default (Design D2).
%    n_ms_starts candidate (omega, kappa) inits are screened via
%    short_iter_ms iterations of Phase B each; the candidate with the
%    lowest final KL objective seeds the full Phase B run.
%  - New optional P fields (safe defaults):
%      P.n_ms_starts     = 4      (number of candidates)
%      P.short_iter_ms   = 20     (Phase B iterations per short trial)
%      P.use_multi_start = true   (opt-out via false; legacy flag aliased)
%  - info struct extended with .best_start (1..n_ms_starts) and
%    .L_hist_ms (short_iter_ms x n_ms_starts L histories).
%  - Driver signature unchanged (6-arg; backward compatible).
%
%  FOUR-PHASE ARCHITECTURE
%  -----------------------
%  Phase A -- Initialisation
%    N0_init from averaged eigenvalue estimate across K_s subcarriers
%    (locked formula, Task 11.3 header).  Convert physical (theta,r) ->
%    (omega,kappa).  Build eta_init (3d+1 real).  Multi-start screen
%    (if P.use_multi_start) selects best of n_ms_starts candidates.
%
%  Phase B -- Block-preconditioned Armijo projected gradient descent
%    Joint update of all 3d+1 components of eta: (omega, kappa, p, N0).
%    Gradient from wb_clkl_estimator; block-diagonal preconditioner
%    (D_prec, see Design Decision D4) equalises step magnitude across
%    parameter groups before Armijo backtracking.
%    Projections: p >= 0, N0 >= 1e-12.
%    L_hist recorded; convergence requires tol_clkl for min_iter=3
%    consecutive iterations (Design Decision D6).
%
%  Phase C -- Post-loop 4-pass alternating 2D residual MF scan
%    Seeds theta_ref/u_ref from the BPD warm-start, NOT from Phase B eta
%    (Design Decision D5).  Score accumulated over all K_s subcarriers.
%    Inherits deflation structure from Paper B nf_clkl.m lines 220-268.
%    Controlled by P.ablation_skip_scan (default false).
%
%  Phase D -- Parameter extraction + info struct assembly.
%
%  DESIGN DECISIONS (documented for Paper C Phase 4 / Sec. III)
%  -------------------------------------------------------------
%  D1 (joint gradient): Phase B updates all 3d+1 components jointly.
%  D2 (multi-start warm-start, v2): n_ms_starts candidate inits screened
%     via short Phase B trial; lowest-L candidate seeds the full run.
%     Default ON.  Opt-out via P.use_multi_start = false.  Candidates:
%     (1) BPD (theta_BPD, r_BPD); (2) (theta_BPD, r_ring) Hussain ring;
%     (3) (theta_BPD, r=r_max); (4) (theta_BPD, r=0.5*(r_lo_fac+r_hi_fac)*r_RD).
%  D3 (Phase C wideband): scan score over all K_s subcarriers.
%  D4 (block preconditioner): (omega,kappa,p,N0) span orders of magnitude
%     in gradient scale.  Per-group RMS normalisation at phase start
%     ensures each group receives comparably-sized steps.  Built per
%     candidate inside multi_start_select and once for the main Phase B.
%  D5 (Phase C from BPD): Phase C seeded from warm-start, not Phase B eta
%     and not multi-start best.  Robust to Phase B kappa stagnation;
%     scan corrects range independently.
%  D6 (consecutive convergence): tol_clkl must hold for min_iter=3
%     consecutive iterations AND t > 5 before declaring convergence.
%  D7 (Phase D KL arbitration): two candidate (theta_ref, u) pairs
%     evaluated; lower-KL pair selected for r_hat.
%
%  INPUTS
%    R_hat_cell : K_s x 1 cell; {k} is N_RF x N_RF raw sample covariance
%                 NOT whitened.
%    W_comb     : M x N_RF frequency-flat constant-modulus combiner
%    theta_init : d x 1 warm-start angles [rad]
%    r_init     : d x 1 warm-start ranges [m]
%    p_init     : d x 1 warm-start powers, or [] for equal-power 1/d
%    P          : parameter struct
%
%  REQUIRED P FIELDS
%    P.M, P.N_RF, P.d, P.K_s, P.lambda_c, P.lambda, P.d_ant
%    P.alpha_k_vec, P.lambda_reg, P.max_iter, P.tol_clkl
%    P.alpha_p, P.ls_beta, P.ls_sigma
%    P.u_min, P.u_max, P.theta_lo, P.theta_hi
%    P.beta_delta, P.c0 (or P.c)
%
%  OPTIONAL P FIELDS (multi-start, Component 2; safe defaults)
%    P.n_ms_starts     -- default 4   (clamped to <=4 canonical candidates)
%    P.short_iter_ms   -- default 20  (Phase B iters per candidate)
%    P.use_multi_start -- default true; legacy alias ablation_multi_start
%
%  OPTIONAL ABLATION FLAGS (default false; do NOT set in production)
%    P.ablation_skip_scan      -- skip Phase C
%    P.ablation_joint_update   -- (default true) set false for p-only update
%    P.ablation_multi_start    -- legacy alias: false => use_multi_start=false
%                                 (otherwise multi-start is ON by default)
%
%  OUTPUTS
%    theta_hat : d x 1  estimated angles [rad]
%    r_hat     : d x 1  estimated ranges [m]
%    p_hat     : d x 1  estimated path powers
%    N0_hat    : scalar estimated noise variance
%    info.L_hist     -- n_iter x 1 main Phase B L history
%    info.n_iter     -- main Phase B iteration count
%    info.N0_init    -- initial N0 estimate (eigenvalue floor)
%    info.converged  -- main Phase B convergence flag
%    info.best_start -- multi-start winning candidate (1..n_ms_starts);
%                       == 1 when use_multi_start = false
%    info.L_hist_ms  -- short_iter_ms x n_ms_starts short-trial histories
%                       ([] when use_multi_start = false)
%
%  Author : R. V. Senyuva (Maltepe University)
%  Date   : May 2026 (v2, Task 11.3 rev3 Component 2)

% =========================================================================
%  INPUT VALIDATION
% =========================================================================
assert(numel(theta_init) == P.d, ...
    'wb_clkl_driver: numel(theta_init)=%d but P.d=%d.', numel(theta_init), P.d);
assert(numel(r_init) == P.d, ...
    'wb_clkl_driver: numel(r_init)=%d but P.d=%d.', numel(r_init), P.d);
assert(numel(R_hat_cell) == P.K_s, ...
    'wb_clkl_driver: numel(R_hat_cell)=%d but P.K_s=%d.', numel(R_hat_cell), P.K_s);

theta_init = theta_init(:);
r_init     = r_init(:);

% =========================================================================
%  MULTI-START DEFAULTS (Component 2)
% =========================================================================
if ~isfield(P, 'n_ms_starts'),   P.n_ms_starts   = 4;  end
if ~isfield(P, 'short_iter_ms'), P.short_iter_ms = 20; end
if ~isfield(P, 'use_multi_start')
    % Default: multi-start active unless legacy ablation_multi_start = false
    if isfield(P, 'ablation_multi_start') && ~P.ablation_multi_start
        P.use_multi_start = false;
    else
        P.use_multi_start = true;
    end
end

abl_noscan  = isfield(P, 'ablation_skip_scan')    && P.ablation_skip_scan;
abl_ponly   = isfield(P, 'ablation_joint_update') && ~P.ablation_joint_update;

% =========================================================================
%  PHASE A -- INITIALISATION
% =========================================================================
d    = P.d;
M    = P.M;
N_RF = P.N_RF;
K_s  = P.K_s;

% N0 from K_s-averaged eigenvalue estimate
R_bar = zeros(N_RF, N_RF);
for k = 1:K_s
    R_bar = R_bar + R_hat_cell{k};
end
R_bar      = (R_bar + R_bar') / (2 * K_s);
ev_sorted  = sort(real(eig(R_bar)), 'ascend');
n_noise_ev = max(1, N_RF - d);
N0_init    = max(mean(ev_sorted(1:n_noise_ev)), 1e-12);

% Power initialisation
if isempty(p_init)
    p_init_val = ones(d, 1) / d;
else
    p_init_val = max(p_init(:), 0);
    if sum(p_init_val) < 1e-15
        p_init_val = ones(d, 1) / d;
    end
end

% Physical (theta, r) -> parametric (omega, kappa)
c_lin  = 2 * pi * P.d_ant / P.lambda_c;
c_quad =     pi * P.d_ant^2 / P.lambda_c;

omega_init = c_lin  * cos(theta_init);
u_init     = min(P.u_max, max(P.u_min, 1 ./ r_init));
kappa_init = c_quad * sin(theta_init).^2 .* u_init;

eta = [omega_init; kappa_init; p_init_val; N0_init];

% --- Multi-start screen (Component 2) ---
best_start_idx = 1;
L_hist_ms      = [];
if P.use_multi_start
    [eta, N0_init, best_start_idx, L_hist_ms] = multi_start_select( ...
        R_hat_cell, W_comb, P, eta, N0_init, omega_init, kappa_init, ...
        d, c_lin, c_quad, theta_init, r_init, abl_ponly);
end

% =========================================================================
%  PHASE B -- BLOCK-PRECONDITIONED ARMIJO PROJECTED GRADIENT DESCENT
% =========================================================================

% Evaluate at initial eta (post multi-start) to build preconditioner (D4)
[L_curr, grad_L, ~] = wb_clkl_estimator(eta, R_hat_cell, W_comb, P);

% Build fixed block preconditioner from initial gradient magnitudes
eps_prec = 1e-6;
D_prec   = build_preconditioner(grad_L, d, eps_prec, abl_ponly);

alpha0    = P.alpha_p / N_RF;
L_hist    = nan(P.max_iter, 1);
converged = false;
n_iter    = 0;
min_iter  = 3;    % D6: consecutive convergence guard
consec    = 0;

for t = 1:P.max_iter

    % Gradient selection (D1 / ablation)
    if abl_ponly
        grad_step           = zeros(3*d+1, 1);
        grad_step(2*d+1:3*d) = grad_L(2*d+1:3*d);
    else
        grad_step = grad_L;
    end

    % Preconditioned direction
    dir      = D_prec .* grad_step;
    suff_dec = grad_step' * dir;

    if suff_dec < 1e-15
        % Gradient effectively zero
        converged = true;
        L_hist(t) = L_curr;
        n_iter    = t;
        break;
    end

    % Armijo backtracking
    alpha  = alpha0;
    L_prev = L_curr;

    for ls = 1:25
        eta_try                = eta - alpha * dir;
        eta_try(2*d+1 : 3*d)  = max(0,     eta_try(2*d+1 : 3*d));
        eta_try(3*d+1)         = max(1e-12, eta_try(3*d+1));

        [L_try, ~, ~] = wb_clkl_estimator(eta_try, R_hat_cell, W_comb, P);

        if L_try <= L_prev - P.ls_sigma * alpha * suff_dec
            break;
        end
        alpha = alpha * P.ls_beta;
    end

    eta    = eta_try;
    L_prev = L_curr;

    [L_curr, grad_L, ~] = wb_clkl_estimator(eta, R_hat_cell, W_comb, P);

    L_hist(t) = L_curr;
    n_iter    = t;

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

L_hist = L_hist(1:n_iter);

% Phase B physical result (used if Phase C is skipped)
info.eta_PhaseB    = eta;
[theta_phB, u_phB] = eta_to_physical(eta, d, c_lin, c_quad, P.u_min, P.u_max);
p_ref              = max(0, eta(2*d+1 : 3*d));
N0_converged       = max(1e-12, eta(3*d+1));

% =========================================================================
%  PHASE C -- 4-PASS ALTERNATING 2D RESIDUAL MF SCAN
%  Seeded from BPD warm-start (D5), not Phase B eta or multi-start best.
% =========================================================================
Q_scan_th = 192;
Q_scan_u  = 256;
th_scan   = linspace(P.theta_lo, P.theta_hi, Q_scan_th);
u_scan    = linspace(P.u_min,    P.u_max,    Q_scan_u);

% Seed from BPD warm-start (D5)
theta_ref = theta_init;
u_ref     = u_init;

if ~abl_noscan

    % Per-subcarrier reference atoms for deflation:
    %   d_ref_cell{k}(:,ell) = W' * a_{ell,k}  at current (theta_ref, u_ref)
    % We need per-subcarrier atoms because the scan score uses R_hat_cell{k}
    % matched against atoms at subcarrier k (Design Decision D3 corrected).
    d_ref_cell = cell(K_s, 1);
    for k = 1:K_s
        Dk = zeros(N_RF, d);
        for ell = 1:d
            a        = wb_nf_fresnel_steer(theta_ref(ell), u_ref(ell), ...
                                           P.alpha_k_vec(k), P) * sqrt(M);
            Dk(:,ell) = W_comb' * a;
        end
        d_ref_cell{k} = Dk;
    end

    for pass = 1:4
        for ell = 1:d

            if mod(pass, 2) == 1
                % ---------------------------------------------------------
                %  Odd pass: scan theta, u fixed at u_ref(ell)
                %  Correct score: sum_k d_test_k^H * R_hat_deflated_k * d_test_k
                %  where R_hat_deflated_k = R_hat_cell{k} - sum_{j!=ell} p_j*d_{j,k}*d_{j,k}^H
                % ---------------------------------------------------------
                sc = zeros(1, Q_scan_th);
                for k = 1:K_s
                    % Per-subcarrier deflated covariance
                    R_def_k = R_hat_cell{k};
                    for j = 1:d
                        if j == ell; continue; end
                        dj      = d_ref_cell{k}(:, j);
                        R_def_k = R_def_k - p_ref(j) * (dj * dj');
                    end
                    R_def_k = (R_def_k + R_def_k') / 2;

                    % Scan atoms at subcarrier k
                    Ds  = scan_atoms_theta(th_scan, u_ref(ell), ...
                                           P.alpha_k_vec(k), W_comb, M, P);
                    sc = sc + real(sum(conj(Ds) .* (R_def_k * Ds), 1));
                end
                [~, bi]        = max(sc);
                theta_ref(ell) = th_scan(bi);

            else
                % ---------------------------------------------------------
                %  Even pass: scan u, theta fixed at theta_ref(ell)
                % ---------------------------------------------------------
                sc = zeros(1, Q_scan_u);
                for k = 1:K_s
                    R_def_k = R_hat_cell{k};
                    for j = 1:d
                        if j == ell; continue; end
                        dj      = d_ref_cell{k}(:, j);
                        R_def_k = R_def_k - p_ref(j) * (dj * dj');
                    end
                    R_def_k = (R_def_k + R_def_k') / 2;

                    Ds  = scan_atoms_u(u_scan, theta_ref(ell), ...
                                       P.alpha_k_vec(k), W_comb, M, P);
                    sc = sc + real(sum(conj(Ds) .* (R_def_k * Ds), 1));
                end
                [~, bi]    = max(sc);
                u_ref(ell) = u_scan(bi);
            end

            % Update per-subcarrier reference atoms for source ell
            for k = 1:K_s
                a = wb_nf_fresnel_steer(theta_ref(ell), u_ref(ell), ...
                                        P.alpha_k_vec(k), P) * sqrt(M);
                d_ref_cell{k}(:, ell) = W_comb' * a;
            end

        end
    end

end  % ~abl_noscan

% =========================================================================
%  PHASE D -- PARAMETER EXTRACTION
%
%  Range selection via KL arbitration (Design Decision D7):
%    The Phase C u-scan MF score is unreliable for range in the compressed
%    domain: the spatial Fresnel quadratic phase is small at the operating
%    point (r ~ r_RD, M=64, 28 GHz), making the coherence peak in u very
%    broad (3dB width >> u_scan range). The scan score is dominated by the
%    noise floor and can peak at an arbitrary u.
%
%    The KL objective IS sensitive to kappa (5x kappa difference between
%    u_scan and u_bpd corresponds to a measurable KL gap). After Phase C
%    has refined theta, we therefore evaluate the KL at two candidate
%    (theta_ref, u) pairs:
%      Candidate A: (theta_ref, u_ref)  -- Phase C scan u
%      Candidate B: (theta_ref, u_init) -- BPD warm-start u (preserved)
%    and select whichever achieves the lower KL objective.
%    Theta is always taken from the Phase C scan (theta refinement is
%    reliable because the linear phase omega is large and well-discriminated).
% =========================================================================
if abl_noscan
    theta_hat = theta_phB(:);
    r_hat     = (1 ./ u_phB(:));
else
    % Theta always from Phase C scan
    theta_out = theta_ref(:);

    % Build eta at candidate A: (theta_ref, u_ref from scan)
    omega_out = c_lin  * cos(theta_out);
    kappa_A   = c_quad * sin(theta_out).^2 .* u_ref(:);
    eta_A     = [omega_out; kappa_A; p_ref(:); N0_converged];

    % Build eta at candidate B: (theta_ref, u_init from BPD)
    kappa_B   = c_quad * sin(theta_out).^2 .* u_init(:);
    eta_B     = [omega_out; kappa_B; p_ref(:); N0_converged];

    % Evaluate KL objective at both candidates (single call each)
    %[L_A, ~, ~] = wb_clkl_estimator(eta_A, R_hat_cell, W_comb, P);
    %[L_B, ~, ~] = wb_clkl_estimator(eta_B, R_hat_cell, W_comb, P);

    %if L_A <= L_B
    %    theta_hat = theta_out;
    %    r_hat     = (1 ./ u_ref(:));
    %else
    %    theta_hat = theta_out;
    %    r_hat     = (1 ./ u_init(:));
    %end

    % === Component 3 (Task 11.3 rev3): BPD-anchor candidate C ===
    % Phase C theta_ref drifts ~0.1 deg from theta_init; that arc-minute
    % shift biases L_A and L_B against BPD's correct geometry by 4-12 nats.
    % Cand C re-anchors at the original BPD warm-start (theta_init, u_init)
    % to provide a Phase D safeguard. Forensic confirmed at mc=1 and mc=30,
    % SNR=10 dB convergence-mode. See Roadmap v18.1 (lessons L48' / L49).
    omega_C   = c_lin  * cos(theta_init(:));
    kappa_C   = c_quad * sin(theta_init(:)).^2 .* u_init(:);
    eta_C     = [omega_C; kappa_C; p_ref(:); N0_converged];
    % Evaluate KL objective at all three candidates (single call each)
    [L_A, ~, ~] = wb_clkl_estimator(eta_A, R_hat_cell, W_comb, P);
    [L_B, ~, ~] = wb_clkl_estimator(eta_B, R_hat_cell, W_comb, P);
    [L_C, ~, ~] = wb_clkl_estimator(eta_C, R_hat_cell, W_comb, P);
    % Three-way argmin selection (replaces 2-way L_A vs L_B)
    [~, sel_idx] = min([L_A, L_B, L_C]);
    switch sel_idx
        case 1   % Cand A: Phase-C theta + Phase-C u
            theta_hat = theta_out;
            r_hat     = (1 ./ u_ref(:));
        case 2   % Cand B: Phase-C theta + BPD u
            theta_hat = theta_out;
            r_hat     = (1 ./ u_init(:));
        case 3   % Cand C: BPD theta + BPD u (anchor)
            theta_hat = theta_init(:);
            r_hat     = (1 ./ u_init(:));
    end
    info.L_PhaseD      = [L_A, L_B, L_C];
    PhaseD_labels      = 'ABC';
    info.PhaseD_select = PhaseD_labels(sel_idx);
end

p_hat  = p_ref(:);
N0_hat = N0_converged;

info.L_hist     = L_hist;
info.n_iter     = n_iter;
info.N0_init    = N0_init;
info.converged  = converged;
info.best_start = best_start_idx;
info.L_hist_ms  = L_hist_ms;

end  % wb_clkl_driver


% =========================================================================
%  LOCAL: build_preconditioner
%  D_prec(i) = 1 / max(rms_of_block(i), eps_prec).
% =========================================================================
function D = build_preconditioner(grad_L, d, eps_prec, abl_ponly)
D = ones(3*d+1, 1);
if abl_ponly
    D(1      : d)    = 0;
    D(d+1    : 2*d)  = 0;
    D(3*d+1)         = 0;
    D(2*d+1  : 3*d)  = 1 / max(rms(grad_L(2*d+1:3*d)), eps_prec);
    return;
end
D(1      : d)    = 1 / max(rms(grad_L(1      : d)),    eps_prec);
D(d+1    : 2*d)  = 1 / max(rms(grad_L(d+1    : 2*d)),  eps_prec);
D(2*d+1  : 3*d)  = 1 / max(rms(grad_L(2*d+1  : 3*d)),  eps_prec);
D(3*d+1)         = 1 / max(abs(grad_L(3*d+1)),          eps_prec);
end


% =========================================================================
%  LOCAL: eta_to_physical
% =========================================================================
function [theta, u] = eta_to_physical(eta, d, c_lin, c_quad, u_min, u_max)
omega = eta(1   : d);
kappa = eta(d+1 : 2*d);
arg   = min(1 - 1e-6, max(1e-6, omega / c_lin));
theta = acos(arg);
sin2  = sin(theta).^2;
u     = kappa ./ (c_quad * sin2);
u     = min(u_max, max(u_min, u));
end


% =========================================================================
%  LOCAL: scan_atoms_theta   (N_RF x Q theta-sweep atoms at subcarrier k)
% =========================================================================
function Ds = scan_atoms_theta(th_scan, u_fix, alpha_k, W_comb, M, P)
Q  = numel(th_scan);
Ds = zeros(size(W_comb, 2), Q);
for qi = 1:Q
    a       = wb_nf_fresnel_steer(th_scan(qi), u_fix, alpha_k, P) * sqrt(M);
    Ds(:,qi) = W_comb' * a;
end
end


% =========================================================================
%  LOCAL: scan_atoms_u   (N_RF x Q u-sweep atoms at subcarrier k)
% =========================================================================
function Ds = scan_atoms_u(u_scan, theta_fix, alpha_k, W_comb, M, P)
Q  = numel(u_scan);
Ds = zeros(size(W_comb, 2), Q);
for qi = 1:Q
    a       = wb_nf_fresnel_steer(theta_fix, u_scan(qi), alpha_k, P) * sqrt(M);
    Ds(:,qi) = W_comb' * a;
end
end


% =========================================================================
%  LOCAL: multi_start_select   (Component 2; replaces v1 zero-iter stub)
%
%  Builds n_ms_starts canonical candidate eta vectors at the BPD
%  theta_init (angular warm-start retained; range varied), runs
%  short_iter_ms iterations of the same Armijo-preconditioned Phase B
%  loop from each, and returns the candidate eta with the lowest final
%  KL objective.
%
%  Candidates (in order):
%    1. BPD: (omega_init, kappa_init)                -- baseline; always 1st
%    2. Ring: u_ring = 1 / (Z_delta * sin^2 theta), Hussain beam-depth
%    3. r=r_max: u = u_min for each path             (far-field end)
%    4. r-mid: r = 0.5*(r_lo_fac+r_hi_fac)*r_RD, clamped to [u_min,u_max]
%             (r-domain midpoint; Component 2-rev)
%
%  P.n_ms_starts is silently clamped to <=4 (canonical set size).
% =========================================================================
function [eta_best, N0_out, best_idx, L_hist_ms] = multi_start_select( ...
    R_hat_cell, W_comb, P, eta_bpd, N0_in, omega_init, kappa_init, ...
    d, c_lin, c_quad, theta_init, r_init, abl_ponly) %#ok<INUSD>

n_short      = P.short_iter_ms;
eps_prec     = 1e-6;
alpha0       = P.alpha_p / P.N_RF;
min_iter_loc = 3;   % D6 consec guard inside short trial

% --- Build canonical candidate kappa vectors ---------------------------
% Use theta_init directly (sin^2) -- avoids the acos numerical inversion
% of omega_init/c_lin near the +/-1 endpoints.
sin2_bpd = sin(theta_init(:)).^2;

D_ap    = (P.M - 1) * P.d_ant;
Z_delta = D_ap^2 / (2 * P.beta_delta^2 * P.lambda_c);

u_ring = zeros(d, 1);
for ii = 1:d
    if sin2_bpd(ii) > 1e-6
        u_ring(ii) = min(P.u_max, max(P.u_min, ...
            1 / (Z_delta * sin2_bpd(ii))));
    else
        u_ring(ii) = P.u_min;
    end
end

% --- Validate required P fields for r-midpoint candidate (Component 2-rev) --
if ~isfield(P, 'r_lo_fac') || ~isfield(P, 'r_hi_fac') || ~isfield(P, 'r_RD')
    error('multi_start_select: P missing one of r_lo_fac, r_hi_fac, r_RD required for r-midpoint candidate (Component 2-rev).');
end

u_rmax    = P.u_min * ones(d, 1);            % r = r_max  (u = u_min)
r_lo_grid = P.r_lo_fac * P.r_RD;            % physical lower grid edge
r_hi_grid = P.r_hi_fac * P.r_RD;            % physical upper grid edge
r_mid     = 0.5 * (r_lo_grid + r_hi_grid);  % r-domain midpoint
u_rmid    = (1 / r_mid) * ones(d, 1);       % convert to u
u_rmid    = min(P.u_max, max(P.u_min, u_rmid)); % clamp to legal u-range

kappa_bpd  = kappa_init;                  % Cand. 1
kappa_ring = c_quad * sin2_bpd .* u_ring; % Cand. 2
kappa_rmax = c_quad * sin2_bpd .* u_rmax; % Cand. 3
kappa_rmid = c_quad * sin2_bpd .* u_rmid; % Cand. 4

kappa_set = {kappa_bpd, kappa_ring, kappa_rmax, kappa_rmid};
n_starts  = min(P.n_ms_starts, numel(kappa_set));   % clamp <=4

% Common power and N0 for all candidates (start from BPD)
p_init0 = eta_bpd(2*d+1 : 3*d);

% Storage
L_hist_ms = nan(n_short, n_starts);
eta_after = cell(n_starts, 1);
L_final   = inf(n_starts, 1);

% --- Run short Phase B trial from each candidate -----------------------
for si = 1:n_starts

    % Build initial eta for candidate si
    eta_si = [omega_init; kappa_set{si}; p_init0; N0_in];

    % Initial estimator evaluation (gradient for preconditioner)
    [L_si, grad_si, ~] = wb_clkl_estimator(eta_si, R_hat_cell, W_comb, P);
    D_si = build_preconditioner(grad_si, d, eps_prec, abl_ponly);

    % Short Armijo loop -- mirrors main Phase B (max_iter -> n_short)
    consec    = 0;
    L_curr_si = L_si;
    grad_curr = grad_si;

    for tt = 1:n_short

        % Gradient selection (matches main loop semantics)
        if abl_ponly
            grad_step              = zeros(3*d+1, 1);
            grad_step(2*d+1 : 3*d) = grad_curr(2*d+1 : 3*d);
        else
            grad_step = grad_curr;
        end

        dir      = D_si .* grad_step;
        suff_dec = grad_step' * dir;

        if suff_dec < 1e-15
            % Effective stationary point reached; freeze rest of L_hist row
            L_hist_ms(tt:end, si) = L_curr_si;
            break;
        end

        % Armijo backtracking
        alpha   = alpha0;
        L_ref   = L_curr_si;
        eta_try = eta_si;            % init for safety if line search fails

        for ls = 1:25
            eta_try                = eta_si - alpha * dir;
            eta_try(2*d+1 : 3*d)  = max(0,     eta_try(2*d+1 : 3*d));
            eta_try(3*d+1)         = max(1e-12, eta_try(3*d+1));

            [L_try, ~, ~] = wb_clkl_estimator(eta_try, R_hat_cell, W_comb, P);

            if L_try <= L_ref - P.ls_sigma * alpha * suff_dec
                break;
            end
            alpha = alpha * P.ls_beta;
        end

        eta_si    = eta_try;
        L_prev_si = L_curr_si;

        [L_curr_si, grad_curr, ~] = wb_clkl_estimator(eta_si, ...
                                                      R_hat_cell, W_comb, P);

        L_hist_ms(tt, si) = L_curr_si;

        rel_change = abs(L_curr_si - L_prev_si) / (abs(L_prev_si) + 1e-15);
        if rel_change < P.tol_clkl && tt > 5
            consec = consec + 1;
        else
            consec = 0;
        end
        if consec >= min_iter_loc
            % Early termination: fill remaining L_hist row with current L
            if tt + 1 <= n_short
                L_hist_ms(tt+1:end, si) = L_curr_si;
            end
            break;
        end
    end

    eta_after{si} = eta_si;
    L_final(si)   = L_curr_si;
end

% --- Select best candidate ---------------------------------------------
[~, best_idx] = min(L_final);
eta_best      = eta_after{best_idx};
N0_out        = max(1e-12, eta_best(3*d+1));

end
