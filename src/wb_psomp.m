function [theta_hat, r_hat, info] = wb_psomp(R_hat_cell, W_comb, P)
%WB_PSOMP  Wideband Polar-grid Simultaneous OMP (WB-P-SOMP) -- Baseline B2.
%
%  Paper C Phase 2, Task 11.6.
%  Wideband extension of nf_psomp.m (Paper B, Task 2 v3).
%
%  ARCHITECTURE OVERVIEW
%  ---------------------
%  The narrowband P-SOMP builds ONE dictionary D_polar (N_RF x S) and
%  computes a single covariance-domain score per atom.  This wideband
%  version decomposes the extension into two independently reused blocks:
%
%  Block 1 -- POLAR GRID (VERBATIM from nf_psomp lines 63-128, v3):
%    Per-angle Hussain beam-depth range sampling.  The polar grid is
%    FREQUENCY-INDEPENDENT.  theta_atoms (1 x S) and r_atoms (1 x S)
%    are built ONCE from the geometry parameters in P.
%
%  Block 2 -- WIDEBAND COVARIANCE-DOMAIN OMP:
%    For the score computation, per-subcarrier dictionaries D_k (N_RF x S)
%    are built ON-THE-FLY inside the OMP loop using wb_nf_fresnel_steer
%    at each alpha_k.  The wideband score aggregates over all K_s
%    subcarriers:
%
%      score(s) = sum_{k=1}^{K_s}  real( d_{s,k}^H * R_hat_cell{k} * d_{s,k} )
%
%    where d_{s,k} = W_comb' * (sqrt(M) * wb_nf_fresnel_steer(theta_s, 1/r_s, alpha_k, P))
%
%  DEFLATION:
%    After each OMP iteration ell, the orthogonal-projection deflation is
%    applied PER SUBCARRIER:
%
%      R_res_cell{k} = P_perp_k * R_hat_cell{k} * P_perp_k'
%
%    where P_perp_k = I - D_sel_k * pinv(D_sel_k) uses the ell x k-specific
%    selected compressed atoms D_sel_k (N_RF x ell) at subcarrier k.
%
%  SIGNAL MODEL (Branch B locked, Phase 1 C3 decision memo):
%    Phase-only USW; NO element-amplitude factor.
%    Atom:  [a_{l,k}]_m = exp( j*alpha_k*(omega_l*m_bar - kappa_l*m_bar^2) )
%    omega_l = (2*pi*d_ant/lambda_c)*cos(theta)    [positive cosine]
%    kappa_l = (pi*d_ant^2/lambda_c)*sin(theta)^2/r
%    alpha_k = f_k/fc = P.alpha_k_vec(k)
%
%  KEY DIFFERENCES vs. nf_psomp (Paper B):
%    nf_psomp:  single R_hat (N_RF x N_RF), single D_polar, alpha_k = 1.
%    wb_psomp:  R_hat_cell (K_s x 1 cell), D_k built per-subcarrier,
%               score = sum_k d_{s,k}^H * R_hat_cell{k} * d_{s,k}.
%    'lambda'  -> 'lambda_c'    (Paper C notation lock).
%    variable 'L' (Paper B)  -> 'd'  (Paper C notation lock).
%
%  WHAT NOT TO DO:
%    - Never use R_mean = (1/K_s)*sum_k R_hat_cell{k} for scoring (Bug B2).
%    - Never rebuild theta_atoms / r_atoms per subcarrier; the grid is
%      frequency-independent.
%    - Never add element-amplitude factor (Branch B locked).
%    - Never call wb_channel_gen_ofdm_nf from inside this function.
%    - Never whiten R_hat_cell here; caller provides raw sample covariance.
%    - Do not use P.lambda alone; always use P.lambda_c.
%
%  INPUTS
%  ------
%  R_hat_cell : K_s x 1 cell array; R_hat_cell{k} is (N_RF x N_RF)
%               compressed sample covariance at subcarrier k.
%               Formed by caller as (1/N)*Yk*Yk' (Hermitian-symmetrised).
%  W_comb     : M x N_RF hybrid combiner matrix (constant-modulus,
%               frequency-flat, same for all subcarriers).
%  P          : parameter struct.  Required fields:
%                 Core:    M, N_RF, d, lambda_c, d_ant
%                 Grid:    Q_theta, theta_lo, theta_hi, r_lo_fac, r_hi_fac
%                          (r_RD is derived internally if absent; see below)
%                          beta_delta (default 1.2, unused in Hussain grid
%                          but kept for struct compatibility)
%                 Wideband: K_s, alpha_k_vec (K_s x 1)
%
%  OUTPUTS
%  -------
%  theta_hat : d x 1  estimated angles [rad]
%  r_hat     : d x 1  estimated ranges [m]
%  info      : struct with fields:
%                .selected   -- 1 x d atom indices in the polar dictionary
%                .p_est      -- d x 1 non-negative power estimates
%                .Q_total    -- total dictionary size S
%                .Q_r_used   -- 4 (nominal minimum per angle, for compat.)
%                .theta_all  -- 1 x S atom angles [rad]
%                .r_all      -- 1 x S atom ranges [m]
%                .sampling   -- 'beam_depth_rBD'
%                .K_s_used   -- P.K_s (wideband subcarrier count used)
%
%  DEPENDENCIES
%  ------------
%  wb_nf_fresnel_steer.m  (Phase 2 utility, Task 11.2 suite)
%
%  PAPER C POSITIONING
%  -------------------
%  Baseline B2 in the six-baseline comparison (Paper C Sec. V).
%  Fairness class: compressed-domain (same data access as proposed CL-KL,
%  Baseline B4); see wb_clkl_driver.m for B4.
%
%  REFERENCE
%  ---------
%  Hussain et al., "Near-Field Channel Estimation in Hybrid MIMO Systems:
%  Polar-Domain Sparsity and Beam-Depth Sampling," IEEE TWC 2025.
%  (Algorithm 1: per-angle beam-depth range grid.)
%
%  Author : R. V. Senyuva (Maltepe University)
%  Date   : May 2026

% =========================================================================
% 0. Unpack parameters
% =========================================================================
M        = P.M;
N_RF     = P.N_RF;
d_path   = P.d;             % number of paths ('d' per Paper C notation lock)
lambda_c = P.lambda_c;      % carrier wavelength (never bare P.lambda)
d_ant    = P.d_ant;
Q_theta  = P.Q_theta;
K_s      = P.K_s;
alpha_k_vec = P.alpha_k_vec(:);   % K_s x 1

% Derived geometry
D_ap  = (M - 1) * d_ant;
if isfield(P, 'r_RD')
    r_RD = P.r_RD;
else
    r_RD = 2 * D_ap^2 / lambda_c;
end
r_min = P.r_lo_fac * r_RD;
r_max = P.r_hi_fac * r_RD;

theta_grid = linspace(P.theta_lo, P.theta_hi, Q_theta);   % 1 x Q_theta

% =========================================================================
% 1. Build per-angle beam-depth polar grid (Hussain TWC 2025 Algorithm 1)
%    VERBATIM TRANSFER from nf_psomp lines 63-128.
%    Only changes: 'lambda' -> 'lambda_c',  loop variable 'k' -> 'ni' or
%    'ri' (kept),  'L' does not appear in this block.
%    The grid is FREQUENCY-INDEPENDENT; built once here.
% =========================================================================
theta_atoms = [];   % 1 x S -- angle for each atom
r_atoms     = [];   % 1 x S -- range for each atom

for ni = 1:Q_theta
    theta_n  = theta_grid(ni);
    cos2_n   = cos(theta_n)^2;
    r_EBRD_n = r_RD / 10 * cos2_n;   % EBRD at this angle

    % --- Generate r_BD-based sample points for this angle ---
    r_samps = r_min;
    r_F     = r_min;
    n_samp  = 1;
    while r_F < r_max && n_samp < 8
        % r_BD formula (Hussain Theorem 1) -- valid only inside EBRD
        if r_F >= r_EBRD_n * 0.99
            break;   % at or beyond EBRD: beam depth -> inf, stop stepping
        end
        denom = r_RD * cos2_n * (1 - 100*r_F^2 / (r_RD^2 * cos2_n^2));
        if denom <= 1e-12
            break;
        end
        r_BD_n = 20 * r_F^2 / denom;
        r_next = min(r_F + r_BD_n, r_max);
        % Avoid duplicate samples
        if abs(r_next - r_samps(end)) > 0.01
            r_samps(end+1) = r_next; %#ok<AGROW>
            n_samp = n_samp + 1;
        end
        if abs(r_next - r_max) < 1e-4; break; end
        r_F = r_next;
    end
    % Always include r_max as endpoint
    if abs(r_samps(end) - r_max) > 0.01
        r_samps(end+1) = r_max;
    end

    % --- Enforce minimum 4 samples (pad with uniform if needed) ---
    Qr_n = numel(r_samps);
    if Qr_n < 4
        n_add = 4 - Qr_n;
        for ai = 1:n_add
            extra = r_min + (r_max - r_min) * ai / (n_add + 1);
            if ~any(abs(extra - r_samps) < 0.1)
                r_samps(end+1) = extra; %#ok<AGROW>
            end
        end
        r_samps = sort(r_samps);
    end

    Qr_n = numel(r_samps);
    theta_atoms = [theta_atoms, repmat(theta_n, 1, Qr_n)]; %#ok<AGROW>
    r_atoms     = [r_atoms,     r_samps];                   %#ok<AGROW>
end

S       = numel(theta_atoms);   % total dictionary size
u_atoms = 1 ./ r_atoms;         % 1 x S inverse ranges

% =========================================================================
% 2. Wideband covariance-domain OMP with deflation
%
%  The score for each atom s is accumulated over all K_s subcarriers.
%  Per-subcarrier dictionaries D_k (N_RF x S) are built ON-THE-FLY using
%  wb_nf_fresnel_steer to avoid storing K_s full matrices simultaneously.
%
%  Deflation is applied PER SUBCARRIER: R_res_cell{k} is maintained
%  independently for each k, projecting out the subspace spanned by the
%  selected compressed atoms at that subcarrier.
% =========================================================================

% Initialise residual covariance cell (Hermitian-symmetrised)
R_res_cell = cell(K_s, 1);
for k = 1:K_s
    Rk = R_hat_cell{k};
    R_res_cell{k} = (Rk + Rk') / 2;
end

selected   = zeros(1, d_path);    % selected atom indices
D_sel_cell = cell(K_s, 1);        % per-subcarrier selected atoms (N_RF x ell)
for k = 1:K_s
    D_sel_cell{k} = zeros(N_RF, 0);
end

for ell = 1:d_path
    % ---- Accumulate wideband matched-filter score over subcarriers ----
    % score(s) = sum_{k=1}^{K_s}  d_{s,k}^H * R_res_cell{k} * d_{s,k}
    % Build D_k on-the-fly at each alpha_k; never use R_mean.
    score_total = zeros(1, S);

    for k = 1:K_s
        alpha_k = alpha_k_vec(k);
        Rk_res  = R_res_cell{k};

        % Build compressed dictionary D_k (N_RF x S) at subcarrier k
        % Atom:  d_{s,k} = W_comb' * (sqrt(M) * wb_nf_fresnel_steer(theta_s, u_s, alpha_k, P))
        A_k = zeros(M, S);
        for s = 1:S
            A_k(:,s) = wb_nf_fresnel_steer(theta_atoms(s), u_atoms(s), alpha_k, P);
        end
        A_k   = sqrt(M) * A_k;     % M x S, unnormalised (norm = sqrt(M))
        D_k   = W_comb' * A_k;     % N_RF x S, compressed atoms

        % Covariance-domain score for this subcarrier
        RD_k        = Rk_res * D_k;
        score_k     = real(sum(conj(D_k) .* RD_k, 1));   % 1 x S
        score_total = score_total + score_k;
    end

    % ---- Select best atom ----
    [~, best]    = max(score_total);
    selected(ell) = best;

    % ---- Per-subcarrier orthogonal-projection deflation ----
    for k = 1:K_s
        alpha_k = alpha_k_vec(k);

        % Rebuild compressed atom for selected index at this subcarrier
        a_new = wb_nf_fresnel_steer(theta_atoms(best), u_atoms(best), alpha_k, P);
        a_new = sqrt(M) * a_new;
        d_new = W_comb' * a_new;   % N_RF x 1

        % Append to selected set for this subcarrier
        D_sel_cell{k} = [D_sel_cell{k}, d_new];   % N_RF x ell

        % Orthogonal projector and deflation
        P_proj_k = D_sel_cell{k} * pinv(D_sel_cell{k});
        P_perp_k = eye(N_RF) - P_proj_k;

        Rk_new = P_perp_k * R_hat_cell{k} * P_perp_k';
        R_res_cell{k} = (Rk_new + Rk_new') / 2;
    end
end

% =========================================================================
% 3. Power estimation
%    Use the carrier-frequency atom (alpha_k=1) and the original
%    (un-deflated) covariance R_hat_cell{k_c} for each selected path.
%    Find the subcarrier index closest to alpha_k = 1.
% =========================================================================
[~, k_c] = min(abs(alpha_k_vec - 1.0));   % carrier-subcarrier index

p_hat = zeros(d_path, 1);
for ell = 1:d_path
    s_idx = selected(ell);
    a_el  = wb_nf_fresnel_steer(theta_atoms(s_idx), u_atoms(s_idx), alpha_k_vec(k_c), P);
    a_el  = sqrt(M) * a_el;
    d_el  = W_comb' * a_el;               % N_RF x 1
    nd2   = max(real(d_el' * d_el), 1e-15);
    p_hat(ell) = max(0, real(d_el' * R_hat_cell{k_c} * d_el) / nd2^2);
end

% =========================================================================
% 4. Assemble outputs
% =========================================================================
theta_hat = theta_atoms(selected).';   % d x 1
r_hat     = r_atoms(selected).';       % d x 1

info.selected   = selected;
info.p_est      = p_hat;
info.Q_r_used   = 4;                   % nominal minimum per angle (compat.)
info.Q_total    = S;
info.theta_all  = theta_atoms;
info.r_all      = r_atoms;
info.sampling   = 'beam_depth_rBD';
info.K_s_used   = K_s;

end
