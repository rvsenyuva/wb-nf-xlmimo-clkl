function [L, grad_L, Ry_cell] = wb_clkl_estimator(eta, R_hat_cell, W_comb, P)
%WB_CLKL_ESTIMATOR  Wideband KL covariance fitting objective + analytic gradient.
%
%  Paper C Phase 2, Task 11.3 (v2 -- SNR-adaptive diagonal loading).
%  Core computation for the WB-CL-KL estimator.
%  Called by wb_clkl_driver.m (Task 11.4) at each iteration.
%
%  CHANGE FROM v1:  Fixed 1e-12 diagonal regularisation replaced by
%  SNR-adaptive loading (Park & Gerstoft OJSP 2024):
%    eps_load = P.eps_reg * trace(R_hat_cell{k}) / N_RF
%    Ry_reg   = R_y_k + eps_load * I_NRF
%  This ties the floor to the observed data scale, eliminating the
%  RCOND=1e-16 singularity at high SNR (root cause of Phase 3 B4 failure).
%  See Paper_C_Roadmap_v14.md Sec.9, Task 11.3 (rev2).
%
%  OBJECTIVE (cross-subcarrier KL sum, Paper C Sec. III):
%
%    L(eta) = sum_{k=1}^{K_s} [ log det R_{y,k}  +  trace(R_{y,k}^{-1} R_hat_k) ]
%             + lambda_reg * sum(p)
%
%  MODEL COVARIANCE per subcarrier k:
%
%    R_{y,k}(eta) = W' * A_k * diag(p) * A_k' * W  +  N0 * (W' * W)
%    A_k(:,l) = wb_nf_fresnel_steer(theta_l, u_l, alpha_k, P) * sqrt(M)
%
%  WHITENING NOTE (Design Decision D3, Task 11.2):
%    R_hat_cell{k} is the RAW compressed sample covariance
%    (1/N)*Y_full(:,:,k)*Y_full(:,:,k)',  NOT whitened.
%    Noise term is N0*(W'*W),  NOT N0*I_{N_RF}.
%    This differs from Paper B nf_clkl.m, which receives a whitened
%    covariance (from nf_hybrid_combiner.m) and uses N0*I.
%    See Paper_C_Roadmap_v7.md Sec. 9, Task 11.2 Design Decision D3.
%
%  PARAMETER VECTOR layout (3d+1 real elements):
%    eta(1     : d)    = omega_l = (2*pi*d_ant/lambda_c)*cos(theta_l)
%    eta(d+1   : 2d)   = kappa_l = (pi*d_ant^2/lambda_c)*sin^2(theta_l)/r_l
%    eta(2d+1  : 3d)   = p_l     (path powers, >= 0)
%    eta(3d+1)         = N0      (noise variance, > 0)
%
%  Physical-to-parametric mapping (locked Paper C notation):
%    omega_l = (2*pi*d_ant/lambda_c) * cos(theta_l)
%    kappa_l = (pi*d_ant^2/lambda_c) * sin^2(theta_l) / r_l
%    alpha_k = f_k / fc  (= P.alpha_k_vec(k))
%    u_l     = 1/r_l
%
%  GRADIENT layout (real, same as eta):
%    grad_L(1     : d)   = dL/d_omega_l   (Stoica-Nehorai sum over k)
%    grad_L(d+1   : 2d)  = dL/d_kappa_l   (alpha_k factor included)
%    grad_L(2d+1  : 3d)  = dL/d_p_l + lambda_reg
%    grad_L(3d+1)        = dL/d_N0
%
%  SIGN CONVENTIONS (chain rule through lambda_eff = lambda_c/alpha_k):
%    da_{l,k}/d_omega_l = +j * alpha_k * m_bar .* a_{l,k}   [RULE B]
%    da_{l,k}/d_kappa_l = -j * alpha_k * m_bar.^2 .* a_{l,k} [RULE C]
%    Both derivatives require alpha_k because wb_nf_fresnel_steer uses
%    lambda_eff = lambda_c/alpha_k.
%
%  N0 INITIALISATION HINT (for wb_clkl_driver.m, Task 11.4):
%    R_bar      = (1/K_s) * sum_k R_hat_cell{k}
%    ev_sorted  = sort(real(eig(R_bar)), 'ascend')
%    n_noise_ev = max(1, N_RF - d)
%    N0_init    = max(mean(ev_sorted(1:n_noise_ev)), 1e-12)
%
%  INPUTS
%    eta        : (3d+1) x 1  real parameter vector
%    R_hat_cell : K_s x 1 cell; {k} is N_RF x N_RF Hermitian (raw, NOT whitened)
%    W_comb     : M x N_RF  frequency-flat constant-modulus combiner
%    P          : parameter struct.  Required fields:
%                   P.M          -- ULA element count
%                   P.N_RF       -- RF chain count
%                   P.d          -- path count
%                   P.lambda_c   -- carrier wavelength c0/fc  [m]
%                   P.lambda     -- alias: must equal P.lambda_c
%                   P.d_ant      -- element spacing lambda_c/2  [m]
%                   P.K_s        -- number of used subcarriers
%                   P.alpha_k_vec-- K_s x 1 frequency scaling f_k/fc
%                   P.lambda_reg -- L1 regularisation weight (default 1e-3)
%                   P.u_min      -- lower inverse-range bound  [1/m]
%                   P.u_max      -- upper inverse-range bound  [1/m]
%                   P.c0         -- speed of light 3e8  [m/s]
%                   P.eps_reg    -- SNR-adaptive loading coeff (default 1e-3)
%                                   Ry_reg = R_y_k + eps_reg*trace(R_hat_k)/N_RF*I
%
%  OUTPUTS
%    L          : scalar  wideband KL objective value (real, finite)
%    grad_L     : (3d+1) x 1  real analytic gradient w.r.t. eta
%    Ry_cell    : K_s x 1 cell; {k} is N_RF x N_RF model covariance at eta
%
%  CALLED BY
%    wb_clkl_driver.m (Task 11.4)  -- optimisation loop
%    wb_clkl_toy_test.m             -- finite-difference gradient verification
%
%  DEPENDENCIES
%    wb_nf_fresnel_steer.m  (Task 11.1 shared Phase 2 utility)
%
%  Author: R. V. Senyuva (Maltepe University)
%  Date  : May 2026 (v2: SNR-adaptive diagonal loading, Task 11.3 rev2)

% =========================================================================
%  PHASE A -- Unpack eta, apply floors, recover physical parameters
% =========================================================================

d      = P.d;
M      = P.M;
N_RF   = P.N_RF;
K_s    = P.K_s;

% Unpack parameter vector
omega  = eta(1     : d);           % d x 1  linear-phase coefficients
kappa  = eta(d+1   : 2*d);         % d x 1  quadratic-phase coefficients
p      = eta(2*d+1 : 3*d);         % d x 1  path powers
N0_raw = eta(3*d+1);               % scalar noise variance

% Apply physical floors  (before any matrix inversion)
p      = max(p,  0);               % powers are non-negative
N0     = max(N0_raw, 1e-12);       % noise floor for numerical safety

% Physical parameter constants (Paper C locked notation)
c_lin  = 2 * pi * P.d_ant / P.lambda_c;   % (2*pi*d_ant/lambda_c)
c_quad =     pi * P.d_ant^2 / P.lambda_c; % (pi*d_ant^2/lambda_c)

% Recover theta_l from omega_l  (clamp argument to stay in (0, pi/2))
%   omega_l = c_lin * cos(theta_l)  =>  cos(theta_l) = omega_l / c_lin
arg_theta   = omega / c_lin;
arg_theta   = min(1 - 1e-6, max(1e-6, arg_theta));   % clamp to (0,1) exclusive
theta_l     = acos(arg_theta);                         % d x 1  [rad]

% Recover u_l = 1/r_l from kappa_l
%   kappa_l = c_quad * sin^2(theta_l) * u_l
sin2_l = sin(theta_l).^2;                              % d x 1
u_l    = kappa ./ (c_quad * sin2_l);                   % d x 1
u_l    = min(P.u_max, max(P.u_min, u_l));              % clamp to valid range

% Precompute quantities used in every iteration of the k-loop
m_bar  = ((0:M-1).' - (M-1)/2);           % M x 1  centred element indices
m_bar2 = m_bar.^2;                         % M x 1  element-index squares
WtW    = W_comb' * W_comb;                 % N_RF x N_RF  (frequency-flat, computed once)
I_NRF  = eye(N_RF);

% =========================================================================
%  Initialise accumulators
% =========================================================================
L           = 0;
grad_omega  = zeros(d, 1);
grad_kappa  = zeros(d, 1);
grad_p      = zeros(d, 1);
grad_N0_acc = 0;
Ry_cell     = cell(K_s, 1);

% =========================================================================
%  PHASE B -- k-loop: accumulate objective and gradient over K_s subcarriers
% =========================================================================

for k = 1:K_s
    alpha_k = P.alpha_k_vec(k);    % scalar frequency scaling f_k / fc

    % ----------------------------------------------------------------
    %  Build per-subcarrier steering matrix A_k (M x d)
    %  Columns are UNNORMALISED atoms: norm(A_k(:,l)) = sqrt(M).
    %  wb_nf_fresnel_steer returns unit-norm; multiply by sqrt(M).
    % ----------------------------------------------------------------
    A_k = zeros(M, d);
    for l = 1:d
        A_k(:,l) = wb_nf_fresnel_steer(theta_l(l), u_l(l), alpha_k, P) * sqrt(M);
    end

    % Compressed atoms  D_k = W' * A_k  (N_RF x d)
    D_k = W_comb' * A_k;

    % ----------------------------------------------------------------
    %  Model covariance  [RULE A: noise term is N0*(W'W), NOT N0*I]
    % ----------------------------------------------------------------
    R_y_k = D_k * diag(p) * D_k'  +  N0 * WtW;
    R_y_k = (R_y_k + R_y_k') / 2;          % enforce Hermitian symmetry

    % Regularised version for stable inversion (SNR-adaptive, Park & Gerstoft 2024)
    %   eps_load scales with trace(R_hat_k)/N_RF = average observed eigenvalue,
    %   which grows with SNR.  This maintains a non-trivial regularisation
    %   floor relative to matrix scale across the full SNR range.
    %   [v1 used fixed 1e-12, which failed at high SNR: RCOND=1e-16 at SNR=20 dB]
    eps_load = P.eps_reg * trace(R_hat_cell{k}) / N_RF;
    Ry_reg   = R_y_k + eps_load * I_NRF;

    % KL objective contribution: log det + trace term
    %   Use try-catch mirroring Paper B nf_clkl kl_obj() pattern
    try
        L = L + real( log(det(Ry_reg)) + trace(Ry_reg \ R_hat_cell{k}) );
    catch
        L = Inf;
        grad_L  = zeros(3*d+1, 1);
        Ry_cell = cell(K_s, 1);
        return
    end

    % ----------------------------------------------------------------
    %  G-matrix (Stoica-Nehorai 1990 pattern)
    %    G_k = R_{y,k}^{-1}  -  R_{y,k}^{-1} * R_hat_k * R_{y,k}^{-1}
    % ----------------------------------------------------------------
    Ry_inv = inv(Ry_reg);                    %#ok<MINV>  (N_RF <= 32: direct inv OK)
    G_k    = Ry_inv - Ry_inv * R_hat_cell{k} * Ry_inv;

    % ----------------------------------------------------------------
    %  Per-path gradient accumulation
    % ----------------------------------------------------------------
    for l = 1:d
        a_lk = A_k(:, l);      % M x 1 unnormalised steering vector at subcarrier k
        d_lk = D_k(:, l);      % N_RF x 1 compressed atom

        % -- omega_l gradient ------------------------------------------
        %   da_{l,k}/d_omega_l = +j * alpha_k * m_bar .* a_{l,k}  [RULE B]
        da_domega = 1j * alpha_k * m_bar .* a_lk;   % M x 1  [alpha_k mandatory]
        Wd_om = W_comb' * da_domega;                 % N_RF x 1
        grad_omega(l) = grad_omega(l) + ...
            2 * real( p(l) * (d_lk' * G_k * Wd_om) );

        % -- kappa_l gradient ------------------------------------------
        %   da_{l,k}/d_kappa_l = -j * alpha_k * m_bar.^2 .* a_{l,k}  [RULE C]
        da_dkappa = -1j * alpha_k * (m_bar2 .* a_lk);   % M x 1
        Wd_kp = W_comb' * da_dkappa;                     % N_RF x 1
        grad_kappa(l) = grad_kappa(l) + ...
            2 * real( p(l) * (d_lk' * G_k * Wd_kp) );

        % -- p_l gradient ----------------------------------------------
        %   dR_{y,k}/d_p_l = d_lk * d_lk'
        %   (lambda_reg term added in Phase C)
        grad_p(l) = grad_p(l) + real(trace((d_lk * d_lk') * G_k));
    end

    % -- N0 gradient -----------------------------------------------
    %   dR_{y,k}/d_N0 = W' * W  (= WtW, precomputed)
    %   dL/d_N0 = sum_k trace(WtW * G_k)
    grad_N0_acc = grad_N0_acc + real(trace(WtW * G_k));

    % Store model covariance for output
    Ry_cell{k} = R_y_k;

end  % k-loop

% =========================================================================
%  PHASE C -- Finalise: add regularisation and assemble gradient
% =========================================================================

L         = L + P.lambda_reg * sum(p);     % L1 regularisation term
grad_p    = grad_p + P.lambda_reg;         % regularisation gradient

% Assemble full gradient vector (same layout as eta)
grad_L = [grad_omega; grad_kappa; grad_p; grad_N0_acc];

end  % wb_clkl_estimator
