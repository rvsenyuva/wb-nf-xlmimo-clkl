function [crb_r, crb_theta, info] = wb_crb_compressed(theta_true, r_true, ...
                                                        p_true, N0, P, opts)
%WB_CRB_COMPRESSED  Wideband compressed-domain CRB for near-field XL-MIMO.
%
%  [crb_r, crb_theta, info] = wb_crb_compressed(theta_true, r_true, ...
%                                                p_true, N0, P)
%  [crb_r, crb_theta, info] = wb_crb_compressed(theta_true, r_true, ...
%                                                p_true, N0, P, opts)
%
%  Extends the GLOBECOM 2026 CRB codebase (wb_crb_globecom2026.m) to:
%    (C4) General d >= 1 paths  (multi-path Slepian-Bangs FIM).
%    (C1) Random-W combiner envelope  (N_seed realisations; mean/std output).
%
%  Signal model (compressed domain, per subcarrier k):
%    y_k(n) = W^H * (A_k * diag(p)^{1/2} * s_k(n) + w(n))
%    R_{y,k} = W^H * (A_k * diag(p) * A_k^H + N0 * I_M) * W
%
%  Parameter vector:
%    eta = [omega_1,...,omega_d, kappa_1,...,kappa_d, p_1,...,p_d, N0]
%    (3d+1) x 1 real.
%
%  Wideband FIM (Assumption W1: cross-subcarrier independence):
%    J_WB = sum_{k=1}^{K_s} J_k,   J_k is (3d+1) x (3d+1).
%
%  CRB is extracted via SVD pseudoinverse with relative threshold 1e-6.
%  Error propagation (Chain rule):
%    CRB_r_ell     = CRB_kappa_ell * (r_ell / kappa_ell)^2
%    CRB_theta_ell = CRB_omega_ell / (d_omega/d_theta_ell)^2
%
%  Random-W envelope (C1):
%    N_seed independent constant-modulus combiners are drawn.
%    crb_r and crb_theta are the seed-averaged CRB vectors (per path).
%    info contains per-seed raw arrays and summary statistics.
%
%  Branch B (locked, Phase 1 C3 decision memo, Apr 29 2026):
%    Phase-only USW steering vector -- no element-amplitude factor.
%    FIM derivatives use:
%      da/d_omega_ell  = +j * alpha_k * m_bar   .* a_{ell,k}  (alpha_k mandatory)
%      da/d_kappa_ell  = -j * alpha_k * m_bar.^2 .* a_{ell,k} (alpha_k mandatory)
%
%  Noise term (locked): dR/dN0 = W^H * W = WtW  (NOT eye(N_RF)).
%
%  INPUTS
%  ------
%  theta_true : d x 1  angles of arrival [rad]
%  r_true     : d x 1  ranges [m]
%  p_true     : d x 1  path powers (linear scale)
%  N0         : scalar noise variance
%  P          : parameter struct with fields:
%                 .M          -- number of ULA elements
%                 .N          -- number of snapshots
%                 .N_RF       -- number of RF chains
%                 .lambda_c   -- carrier wavelength [m]  (= P.c / P.fc)
%                 .d_ant      -- element spacing [m]     (= P.lambda_c / 2)
%                 .K_s        -- number of subcarriers used in CRB
%                 .alpha_k_vec -- K_s x 1 frequency ratios f_k/f_c
%               Optional:
%                 .N_seed     -- number of W realisations (default: 50)
%                 .rng_seed_W -- base RNG seed for W draws (default: 0)
%
%  opts (optional struct):
%    .return_J    -- if true, return J_WB in info.J_WB (default: false)
%    .verbose     -- if true, print per-seed diagnostics (default: false)
%
%  OUTPUTS
%  -------
%  crb_r     : d x 1  seed-averaged sqrt-CRB for range [m]
%               (= sqrt(E_W[ CRB_r_ell ]) element-wise, averaged over seeds)
%  crb_theta : d x 1  seed-averaged sqrt-CRB for angle [deg]
%  info      : struct with fields:
%               .crb_r_all    -- d x N_seed  sqrt(CRB_r) per path per seed [m]
%               .crb_theta_all -- d x N_seed sqrt(CRB_theta) per path per seed [deg]
%               .crb_r_mean   -- d x 1  mean over seeds [m]
%               .crb_r_std    -- d x 1  std  over seeds [m]
%               .crb_r_p10    -- d x 1  10th-percentile [m]
%               .crb_r_p90    -- d x 1  90th-percentile [m]
%               .crb_theta_mean -- d x 1 [deg]
%               .crb_theta_std  -- d x 1 [deg]
%               .cond_J_all   -- 1 x N_seed  condition number of J_WB
%               .n_singular_all -- 1 x N_seed  number of near-zero singular values
%               .J_WB         -- (3d+1)x(3d+1) wideband FIM (last seed; only
%                                if opts.return_J = true)
%
%  NOTATION (Paper C locked)
%  -------------------------
%  d        -- number of paths (use 'd', never 'L')
%  omega_ell = (2*pi*d_ant/lambda_c)*cos(theta_ell)   [linear phase slope]
%  kappa_ell = (pi*d_ant^2/lambda_c)*sin^2(theta_ell)/r_ell  [curvature]
%  alpha_k   = f_k / f_c  (P.alpha_k_vec)
%
%  INTERFACE CONTRACT
%  ------------------
%  Called by: run_monte_carlo_paperC.m (Task 11.8) -- ONCE per sweep point.
%  The random-W envelope is computed INSIDE this function (not in outer MC).
%
%  COMPARISON WITH GLOBECOM SCRIPT (wb_crb_globecom2026.m)
%  --------------------------------------------------------
%  GLOBECOM: d=1, fixed W, standalone script with sweeps and CSV export.
%  Task 11.5: d>=1, N_seed random W realisations, pure function (no sweep).
%  Sweep logic (over SNR, B, N_RF, r) is the CALLER's responsibility.
%
%  REFERENCES
%  ----------
%  [1] P. Stoica and A. Nehorai, "Performance study of conditional and
%      unconditional direction-of-arrival estimation," IEEE Trans. ASSP,
%      vol.38, no.10, pp.1783-1795, Oct. 1990.  (Slepian-Bangs FIM)
%  [2] R. V. Senyuva, GLOBECOM 2026, arXiv:2604.08531.
%
%  Author:  R. V. Senyuva (Maltepe University)
%  Date:    May 2026

% =========================================================================
%  0. Parse inputs and set defaults
% =========================================================================
theta_true = theta_true(:);   % d x 1
r_true     = r_true(:);       % d x 1
p_true     = p_true(:);       % d x 1
d          = numel(theta_true);

if nargin < 6 || isempty(opts)
    opts = struct();
end
return_J = isfield(opts, 'return_J')  && opts.return_J;
verbose  = isfield(opts, 'verbose')   && opts.verbose;

N_seed     = 50;
rng_seed_W = 0;
if isfield(P, 'N_seed'),     N_seed     = P.N_seed;     end
if isfield(P, 'rng_seed_W'), rng_seed_W = P.rng_seed_W; end

% =========================================================================
%  1. Derived parameters
% =========================================================================
M        = P.M;
N        = P.N;            % number of snapshots
N_RF     = P.N_RF;
lam_c    = P.lambda_c;
d_ant    = P.d_ant;
K_s      = P.K_s;
alpha_k  = P.alpha_k_vec(:);   % K_s x 1

n_par    = 3*d + 1;           % FIM dimension

% Centred element indices (M x 1)
m_bar    = ((0:M-1).' - (M-1)/2);
m_bar2   = m_bar.^2;

% Pre-compute omega and kappa at carrier frequency (alpha=1)
omega_true = (2*pi*d_ant/lam_c) * cos(theta_true);    % d x 1
c_coef     = (pi*d_ant^2/lam_c) * sin(theta_true).^2; % d x 1
kappa_true = c_coef ./ r_true;                         % d x 1

% Jacobian denominators for error propagation
domega_dtheta = -(2*pi*d_ant/lam_c) * sin(theta_true);  % d x 1
dkappa_dr     = -c_coef ./ r_true.^2;                    % d x 1

% =========================================================================
%  2. Allocate output arrays
% =========================================================================
crb_r_all     = zeros(d, N_seed);   % sqrt(CRB_r) per path per seed [m]
crb_theta_all = zeros(d, N_seed);   % sqrt(CRB_theta) per path per seed [deg]
cond_J_all    = zeros(1, N_seed);
n_sing_all    = zeros(1, N_seed);

% =========================================================================
%  3. Loop over N_seed random W realisations
% =========================================================================
rng_state_save = rng;   % save caller's RNG state

for iseed = 1:N_seed

    % Draw a fresh constant-modulus combiner for this seed
    rng(rng_seed_W + iseed - 1, 'twister');
    W = (1/sqrt(M)) * exp(1j * 2*pi * rand(M, N_RF));   % M x N_RF

    % -----------------------------------------------------------------
    %  Accumulate wideband FIM: J_WB = sum_k J_k
    % -----------------------------------------------------------------
    J_WB = zeros(n_par, n_par);

    for ks = 1:K_s
        ak = alpha_k(ks);   % scalar frequency ratio for this subcarrier

        J_k = loc_compute_fim_multipath( ...
            omega_true, kappa_true, p_true, N0, ak, m_bar, m_bar2, W, N);

        J_WB = J_WB + J_k;
    end

    % -----------------------------------------------------------------
    %  Invert J_WB via SVD pseudoinverse (relative threshold 1e-6)
    % -----------------------------------------------------------------
    [CRB_mat, n_sing, cond_num] = loc_pinv_svd(J_WB);

    cond_J_all(iseed)  = cond_num;
    n_sing_all(iseed)  = n_sing;

    % -----------------------------------------------------------------
    %  Extract marginal CRBs and apply chain rule (per path)
    % -----------------------------------------------------------------
    for ell = 1:d
        var_omega = max(CRB_mat(ell,     ell),   0);  % omega_ell block
        var_kappa = max(CRB_mat(d+ell, d+ell),   0);  % kappa_ell block

        % CRB_r: via kappa chain rule  CRB_r = CRB_kappa * (r/kappa)^2
        crb_r_all(ell, iseed)     = sqrt(var_kappa) / abs(dkappa_dr(ell));

        % CRB_theta: via omega chain rule  CRB_theta = CRB_omega / (domega/dtheta)^2
        crb_theta_all(ell, iseed) = sqrt(var_omega) / abs(domega_dtheta(ell)) ...
                                    * (180/pi);   % radians -> degrees
    end

    if verbose
        fprintf('  seed %3d | cond(J)=%.2e | n_sing=%d | sqrt_crb_r(1)=%.4f m\n', ...
            iseed, cond_num, n_sing, crb_r_all(1, iseed));
    end
end

% Restore caller's RNG state
rng(rng_state_save);

% =========================================================================
%  4. Summary statistics across seeds
% =========================================================================
crb_r     = mean(crb_r_all,     2);   % d x 1 (seed mean)
crb_theta = mean(crb_theta_all, 2);   % d x 1 (seed mean)

info.crb_r_all       = crb_r_all;
info.crb_theta_all   = crb_theta_all;
info.crb_r_mean      = crb_r;
info.crb_r_std       = std(crb_r_all,  0, 2);
info.crb_r_p10       = prctile(crb_r_all,  10, 2);
info.crb_r_p90       = prctile(crb_r_all,  90, 2);
info.crb_theta_mean  = crb_theta;
info.crb_theta_std   = std(crb_theta_all, 0, 2);
info.cond_J_all      = cond_J_all;
info.n_singular_all  = n_sing_all;

if return_J
    % Return J_WB from the LAST seed (deterministic reference)
    rng(rng_seed_W + N_seed - 1, 'twister');
    W_last = (1/sqrt(M)) * exp(1j * 2*pi * rand(M, N_RF));
    rng(rng_state_save);
    J_last = zeros(n_par, n_par);
    for ks = 1:K_s
        J_last = J_last + loc_compute_fim_multipath( ...
            omega_true, kappa_true, p_true, N0, alpha_k(ks), ...
            m_bar, m_bar2, W_last, N);
    end
    info.J_WB = J_last;
end

end   % wb_crb_compressed


% =========================================================================
%  LOCAL FUNCTION: Per-subcarrier Slepian-Bangs FIM (multi-path)
% =========================================================================
function J = loc_compute_fim_multipath(omega, kappa, p, N0, alpha_k, ...
                                        m_bar, m_bar2, W, N_snap)
%LOC_COMPUTE_FIM_MULTIPATH  Per-subcarrier (3d+1) x (3d+1) Slepian-Bangs FIM.
%
%  Extends loc_compute_fim_single from wb_crb_globecom2026.m to d >= 1 paths.
%  The function signature and FIM assembly are identical; only d is general.
%
%  Inputs:
%    omega   - d x 1  spatial frequencies (at carrier, scaled by alpha_k here)
%    kappa   - d x 1  Fresnel curvatures  (at carrier, scaled by alpha_k here)
%    p       - d x 1  path powers
%    N0      - scalar noise variance
%    alpha_k - scalar frequency ratio f_k / f_c
%    m_bar   - M x 1  centred element indices
%    m_bar2  - M x 1  centred element indices squared (pre-computed)
%    W       - M x N_RF  analog combiner
%    N_snap  - scalar  number of snapshots
%
%  Output:
%    J       - (3d+1) x (3d+1)  per-subcarrier FIM
%
%  Branch B (locked): phase-only steering vector; NO amplitude factor.
%  Derivatives (alpha_k MANDATORY in both):
%    da/d_omega_ell  = +j * alpha_k * m_bar   .* a_{ell,k}
%    da/d_kappa_ell  = -j * alpha_k * m_bar.^2 .* a_{ell,k}
%  Noise: dR/dN0 = W^H * W = WtW  (NOT eye(N_RF)).

d_loc    = numel(omega);
N_RF_loc = size(W, 2);
n_par    = 3*d_loc + 1;

WtW = W' * W;

% Build R_{y,k} and cache compressed steering vectors d_{ell,k} = W^H a_{ell,k}
Ry   = N0 * WtW;
d_vec = zeros(N_RF_loc, d_loc);   % d_{ell,k} = W^H a_{ell,k}
a_mat = zeros(numel(m_bar), d_loc);  % unnormalized steering vectors (norm = sqrt(M))

for ell = 1:d_loc
    % Unnormalized near-field steering vector at subcarrier k (Branch B)
    a_ell = exp( 1j * alpha_k * omega(ell) * m_bar ...
               - 1j * alpha_k * kappa(ell) * m_bar2);   % M x 1
    a_mat(:, ell)  = a_ell;
    d_vec(:, ell)  = W' * a_ell;                         % N_RF x 1
    Ry = Ry + p(ell) * (d_vec(:,ell) * d_vec(:,ell)');
end

% Regularise R_{y,k} (Hermitian symmetrisation + floor)
Ry    = (Ry + Ry') / 2 + 1e-12 * eye(N_RF_loc);
Ry_inv = inv(Ry);  %#ok<MINV>

% Build covariance derivatives dR_{y,k}/d_eta_i  (N_RF x N_RF each)
% Parameter ordering: [omega_1..omega_d, kappa_1..kappa_d, p_1..p_d, N0]
dR = cell(n_par, 1);

for ell = 1:d_loc
    a_ell    = a_mat(:, ell);
    Wd_ell   = d_vec(:, ell);   % W^H a_ell

    % d/d(omega_ell): da/d_omega = +j*alpha_k*m_bar .* a_ell
    da_domega     = 1j * alpha_k * (m_bar .* a_ell);
    Wd_ell_omega  = W' * da_domega;
    dR{ell} = p(ell) * (Wd_ell_omega * Wd_ell' + Wd_ell * Wd_ell_omega');

    % d/d(kappa_ell): da/d_kappa = -j*alpha_k*m_bar^2 .* a_ell
    da_dkappa     = -1j * alpha_k * (m_bar2 .* a_ell);
    Wd_ell_kappa  = W' * da_dkappa;
    dR{d_loc + ell} = p(ell) * (Wd_ell_kappa * Wd_ell' + Wd_ell * Wd_ell_kappa');

    % d/d(p_ell): dR/dp_ell = d_ell * d_ell^H
    dR{2*d_loc + ell} = Wd_ell * Wd_ell';
end

% d/d(N0): dR/dN0 = W^H * W = WtW  (NOT I_{N_RF})
dR{n_par} = WtW;

% Assemble FIM (Slepian-Bangs formula, exploit Hermitian symmetry)
J = zeros(n_par, n_par);
for ii = 1:n_par
    Ry_inv_dRi = Ry_inv * dR{ii};
    for jj = ii:n_par
        val = N_snap * real(trace(Ry_inv_dRi * Ry_inv * dR{jj}));
        J(ii, jj) = val;
        J(jj, ii) = val;
    end
end

end   % loc_compute_fim_multipath


% =========================================================================
%  LOCAL FUNCTION: SVD pseudoinverse with diagnostics
% =========================================================================
function [J_inv, n_sing, cond_num] = loc_pinv_svd(J)
%LOC_PINV_SVD  SVD pseudoinverse with relative threshold 1e-6.
%  Matches the implementation in wb_crb_globecom2026.m loc_pinv_svd.
%
%  Additional outputs: n_sing (number of near-zero singular values),
%  cond_num (condition number, ratio of max to min nonzero singular value).

[U, S, V] = svd(J);
s       = diag(S);
tol     = 1e-6 * max(s);
active  = s > tol;
s_inv   = zeros(size(s));
s_inv(active) = 1 ./ s(active);
J_inv   = V * diag(s_inv) * U';

n_sing  = sum(~active);
if any(active)
    cond_num = max(s(active)) / min(s(active));
else
    cond_num = Inf;
end

end   % loc_pinv_svd
