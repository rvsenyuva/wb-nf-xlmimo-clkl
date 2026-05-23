function [theta_hat, r_hat, alpha_hat, diag_info] = bpd_baseline(X_full, params)
%BPD_BASELINE  Cui & Dai 2023 Bilinear Pattern Detection for wideband
%              near-field XL-MIMO channel estimation. Paper C baseline B4.
%
%  [theta_hat, r_hat, alpha_hat, diag_info] = bpd_baseline(X_full, params)
%
%  Implements Algorithm 1 of:
%    M. Cui and L. Dai, "Near-field wideband channel estimation for
%    extremely large-scale MIMO," Sci. China Inf. Sci., vol. 66,
%    p. 172303, Jul. 2023. https://doi.org/10.1007/s11432-022-3654-y
%
%  PAPER C ROLE
%  ------------
%  B4 is the strongest published full-array wideband near-field baseline.
%  It operates on FULL-ARRAY observations X_full (no hybrid combiner),
%  giving it a structural advantage over Paper C's WB-CL-KL which works
%  on compressed observations {R_hat_y_k}. This asymmetry is documented
%  in Paper C Sec. V fairness framing.
%
%  SIGNAL MODEL (Paper C Sec. II / Phase 1 locked notation)
%  ---------------------------------------------------------
%  Per-subcarrier snapshot matrix: X_full(:,:,k) in C^{M x N}
%    x_k(n) = sum_{l=1}^{d} s_{l,k}(n)*a_{l,k}(theta_l,r_l) + w_k(n)
%  Wideband Fresnel steering vector (GLOBECOM eq. a_k):
%    [a_{l,k}]_m = exp(j*alpha_k*(m_bar*omega_l + m_bar^2*kappa_l))
%  where
%    omega_l  = (2*pi*d_ant/lambda_c)*cos(theta_l)  [Paper B sign convention]
%    kappa_l  = (pi*d_ant^2/lambda_c)*sin(theta_l)^2/r_l
%    alpha_k  = f_k/fc   (wideband frequency scaling factor)
%    s_{l,k}  ~ CN(0,p_l), cross-frequency independent (Assumption W1)
%    w_k(n)   ~ CN(0,N0*I_M), white per subcarrier
%
%  SIGN CONVENTION NOTE
%  --------------------
%  nf_fresnel_steer (Paper B) uses omega = +(2*pi*d_ant/lambda_c)*cos(theta).
%  wb_nf_fresnel_steer (called below) inherits this convention.
%  The channel generator (wb_channel_gen_ofdm_nf, Task 11.2) MUST use
%  the same sign convention -- dictionary and data must match.
%
%  DEVIATIONS FROM CUI & DAI 2023 (documented for Paper C Sec. V)
%  ---------------------------------------------------------------
%  D1: Pre-whitening skipped. Cui & Dai pre-whiten because compressed noise
%      C = sigma^2*diag{A_p*A_p^H} is colored. Here A_p = I_M (full array),
%      so C = sigma^2*I -- scalar scaling, cancels in all power ratios.
%  D2: Rectangular polar grid (Day-1). Cui & Dai use a ring-based range grid
%      (Lemma 1) with angle-dependent range samples. Rectangular grid used
%      for Day 1; ring-based grid is a planned Day 2 investigation (NOTE_D2
%      marks). Investigation DEFERRED: Day-1 toy test already passed
%      (|theta_err| < 0.02 deg, |r_err| < 0.11 m at SNR=10 dB) so accuracy
%      gap motivating the ring grid is not yet demonstrated. Revisit after
%      Monte Carlo sweep in Phase 2.
%  D3: Gain estimated via per-subcarrier LS on sample mean (greedy pass),
%      then joint LS polishing across all d_max paths (Section 5b).
%      Standard OMP convention; Cui & Dai give no explicit formula.
%  D4: Residual is M x N x K_s tensor; Cui & Dai stack pilot slots. The
%      bilinear gather and residual update are equivalent in both cases.
%  D5 (Day-2): Joint LS polishing added (Section 5b). After greedy
%      selection fixes the support (na_picked, nd_picked), a single batched
%      LS pass re-estimates all d_max gains jointly per subcarrier. Support
%      is unchanged; only gains are updated. Controlled by bpd_do_ls_polish.
%
%  INPUTS
%  ------
%  X_full  : M x N x K_s complex tensor of full-array OFDM snapshots.
%            X_full(:,:,k) is the M x N snapshot matrix at subcarrier k.
%  params  : struct with fields listed below.
%
%    -- Inherited from Paper B nf_params (REQUIRED) --
%    fc        : carrier frequency [Hz]              (default 28e9)
%    c0        : speed of light [m/s]                (default 3e8)
%    M         : ULA element count                   (default 256)
%    N         : snapshots per subcarrier             (default 64)
%    d_ant     : element spacing [m]                 (= lambda_c/2)
%    lambda_c  : carrier wavelength [m]              (= c0/fc)
%    theta_lo  : min angle [rad], must be >= 15*pi/180
%    theta_hi  : max angle [rad]
%    r_lo_fac  : r_min = r_lo_fac * r_RD
%    r_hi_fac  : r_max = r_hi_fac * r_RD
%    r_RD      : Rayleigh distance [m]               (= 2*D_ap^2/lambda_c)
%    beta_delta: Cui & Dai coherence parameter        (default 1.2)
%    -- Wideband additions (REQUIRED) --
%    K_s       : number of used subcarriers           (default min(K,512))
%    k_indices : K_s x 1 integer vector of subcarrier indices (1-based)
%    Delta_f   : subcarrier spacing [Hz]              (default 120e3)
%    alpha_k_vec : K_s x 1 vector of f_k/fc          (precomputed outside)
%    -- BPD-specific (REQUIRED) --
%    G_theta   : number of polar grid angle bins      (default 128)
%    G_r       : number of polar grid range bins      (default 64)
%    d_max     : max paths to detect (= L_hat in paper)
%    -- Optional --
%    bpd_write_csv   : logical, write CSV diagnostic  (default true)
%    bpd_csv_path    : char, path for diagnostic CSV
%                      (default 'bpd_baseline_diag.csv')
%    bpd_do_ls_polish: logical, run joint LS polishing after greedy
%                      selection (Section 5b)          (default true)
%                      When true, alpha_hat uses gain_polished instead of
%                      gain_store; res_norm_polished column added to CSV.
%
%  OUTPUTS
%  -------
%  theta_hat  : d_max x 1  estimated angles [rad]
%  r_hat      : d_max x 1  estimated ranges [m]
%  alpha_hat  : d_max x 1  complex path gains (carrier-referenced)
%  diag_info  : struct with fields
%    .T_power          d_max x G_theta x G_r  BPD power maps per iteration
%    .res_norm         d_max x 1  residual Frobenius norm after each greedy step
%    .res_norm_polished scalar     residual Frobenius norm after joint LS polish
%                                 (set to NaN if bpd_do_ls_polish == false)
%    .na_picked        d_max x 1  angle grid index selected per iteration
%    .nd_picked        d_max x 1  range grid index selected per iteration
%    .T_peak           d_max x 1  detection statistic at peak per iteration
%    .theta_grid       G_theta x 1  angle grid [rad]
%    .r_grid           G_r x 1    range grid [m]
%    .alpha_k_vec      K_s x 1    frequency scaling factors used
%
%  Author : R. V. Senyuva (Maltepe University)
%  Version: Day-2 (rectangular grid, full-array, no pre-whitening,
%           joint LS polishing)
%  Date   : May 2026

% =========================================================================
%  SECTION 0: INPUT VALIDATION AND PARAMETER UNPACKING
% =========================================================================

[M_in, N_in, Ks_in] = size(X_full);

req = {'fc','c0','M','N','d_ant','lambda_c','theta_lo','theta_hi', ...
       'r_lo_fac','r_hi_fac','r_RD','K_s','k_indices','Delta_f', ...
       'alpha_k_vec','G_theta','G_r','d_max'};
for fi = 1:numel(req)
    assert(isfield(params, req{fi}), ...
        'bpd_baseline: missing required params field: %s', req{fi});
end

M        = params.M;
N_snaps  = params.N;  %#ok<NASGU> -- stored for documentation; N_in used for actual tensor dim
lambda_c = params.lambda_c;
theta_lo = params.theta_lo;
theta_hi = params.theta_hi;
r_lo_fac = params.r_lo_fac;
r_hi_fac = params.r_hi_fac;
r_RD     = params.r_RD;
K_s      = params.K_s;
alpha_k  = params.alpha_k_vec(:);   % K_s x 1, column
G_theta  = params.G_theta;
G_r      = params.G_r;
d_max    = params.d_max;

if isfield(params,'bpd_write_csv');    do_csv      = params.bpd_write_csv;
else;                                  do_csv      = true;           end
if isfield(params,'bpd_csv_path');     csv_path    = params.bpd_csv_path;
else;                                  csv_path    = 'bpd_baseline_diag.csv'; end
if isfield(params,'bpd_do_ls_polish'); do_ls_polish = params.bpd_do_ls_polish;
else;                                  do_ls_polish = true;          end

assert(M_in == M, ...
    'bpd_baseline: X_full dim 1 = %d but params.M = %d', M_in, M);
assert(Ks_in == K_s, ...
    'bpd_baseline: X_full dim 3 = %d but params.K_s = %d', Ks_in, K_s);
assert(numel(alpha_k) == K_s, ...
    'bpd_baseline: alpha_k_vec length %d != K_s = %d', numel(alpha_k), K_s);
assert(theta_lo >= 15*pi/180, ...
    'bpd_baseline: theta_lo = %.1f deg < 15 deg (physical regime guard)', ...
    theta_lo*180/pi);

% =========================================================================
%  SECTION 1: POLAR GRID CONSTRUCTION (RECTANGULAR, DAY-1)
% =========================================================================
% Angle grid: uniform in sin-domain over [theta_lo, theta_hi].
% Range grid: uniform over [r_lo_fac*r_RD, r_hi_fac*r_RD].
%
% NOTE_D2: Day-2 ring-based upgrade DEFERRED (see D2 in header).
%   Cui & Dai (2023) Lemma 1 formula for reference:
%     r_{na,nd} = D_ap^2*cos^2(theta_grid(na))/(2*beta_delta^2*lambda_c*nd)
%   Simplification using r_RD = 2*D_ap^2/lambda_c (eliminates lambda_c):
%     r_{na,nd} = r_RD*cos^2(theta_grid(na)) / (4*beta_delta^2*nd)
%   Implementation requires 2D u_grid (G_theta x G_r) and 3D Lambda map
%   (G_theta x G_r x K_s). Deferred until rectangular-grid NMSE motivates it.

sin_lo   = sin(theta_lo);
sin_hi   = sin(theta_hi);
sin_grid = sin_lo + (0:G_theta-1).' * (sin_hi - sin_lo) / (G_theta - 1);
theta_grid = asin(sin_grid);        % G_theta x 1  [rad]

r_min   = r_lo_fac * r_RD;
r_max   = r_hi_fac * r_RD;
r_grid  = r_min + (0:G_r-1).' * (r_max - r_min) / (G_r - 1); % G_r x 1 [m]
u_grid  = 1 ./ r_grid;             % G_r x 1  [1/m]

% =========================================================================
%  SECTION 2: CARRIER-FREQUENCY DICTIONARY
% =========================================================================
% Psi_mat(:, (nd-1)*G_theta + na) = unit-norm Fresnel atom at (theta_grid(na),
%   u_grid(nd)) evaluated at f_c (alpha_k = 1).
% Built via 3D broadcasting: M x G_theta x G_r, then reshaped.

m_bar        = ((0:M-1).' - (M-1)/2);          % M x 1

omega_grid   = (2*pi*params.d_ant/lambda_c) .* cos(theta_grid);   % G_theta x 1
c_theta_grid = (pi*params.d_ant^2/lambda_c) .* sin(theta_grid).^2;% G_theta x 1

m3  = reshape(m_bar,        M,       1,       1);
om3 = reshape(omega_grid,   1, G_theta,       1);
ct3 = reshape(c_theta_grid, 1, G_theta,       1);
u3  = reshape(u_grid,       1,       1, G_r);

% Phase: m_bar*omega - m_bar^2*c_theta*u  (carrier frequency, alpha_k=1)
% Sign: positive linear, negative quadratic -- matches wb_nf_fresnel_steer
% at alpha_k=1 and the Paper B nf_fresnel_steer convention.
phase_tensor = m3 .* om3 - (m3.^2) .* (ct3 .* u3);  % M x G_theta x G_r
Psi_tensor   = exp(1j * phase_tensor);                % M x G_theta x G_r

norms3     = sqrt(sum(abs(Psi_tensor).^2, 1));        % 1 x G_theta x G_r
Psi_tensor = Psi_tensor ./ norms3;

% Column order: angle-major  j = (nd-1)*G_theta + na
Psi_mat = reshape(Psi_tensor, M, G_theta * G_r);      % M x (G_theta*G_r)

% =========================================================================
%  SECTION 3: BILINEAR PATTERN INDEX MAPS
% =========================================================================
% Gamma(na, k)  = nearest angle-grid index for alpha_k(k)*sin_grid(na)
% Lambda(nd, k) = nearest range-grid index for alpha_k(k)*u_grid(nd)
% Both precomputed as G_theta x K_s and G_r x K_s integer matrices.
% This implements Cui & Dai (2023) Eqs. 15-16.

sin_scaled = sin_grid * alpha_k.';   % G_theta x K_s  (outer product)

Gamma = zeros(G_theta, K_s);
for k = 1:K_s
    sc  = sin_scaled(:, k);
    sc  = max(sin_lo, min(sin_hi, sc));   % clamp to grid bounds
    idx = round(1 + (sc - sin_lo) / (sin_hi - sin_lo) * (G_theta - 1));
    Gamma(:, k) = max(1, min(G_theta, idx));
end

% Lambda(nd, k) = nearest r_grid index to r_grid(nd) / alpha_k(k).
% Derivation: the bilinear-shifted range at subcarrier k for a path at
% r_grid(nd) is r_shifted = r_grid(nd) / alpha_k(k), because
%   u_shifted = alpha_k * u_grid(nd) = alpha_k / r_grid(nd)
%   r_shifted = 1 / u_shifted = r_grid(nd) / alpha_k(k).
% Since r_grid is linearly spaced, the nearest-index formula is linear in r.
%
% NOTE: the earlier u-domain formula
%   idx = round(1 + (u_hi - u_shifted)/(u_hi - u_lo) * (G_r-1))
% is WRONG for a linear r_grid because u_grid = 1/r_grid is NOT linearly
% spaced in u.  Using it maps nd=22 (r=6.1 m) to idx=30 (r=8.0 m) at
% alpha_k~1.0001 -- a 1.9 m error -- causing the power map peak to appear
% at the wrong range index.  The r-domain formula below is exact.

r_scaled = r_grid * (1 ./ alpha_k.');   % G_r x K_s: r_shifted(nd,k) = r_grid(nd)/alpha_k(k)

Lambda = zeros(G_r, K_s);
for k = 1:K_s
    rs  = r_scaled(:, k);                           % G_r x 1
    rs  = max(r_min, min(r_max, rs));               % clamp to grid bounds
    idx = round(1 + (rs - r_min) / (r_max - r_min) * (G_r - 1));
    Lambda(:, k) = max(1, min(G_r, idx));
end

% =========================================================================
%  SECTION 4: INITIALISE RESIDUAL TENSOR
% =========================================================================
% R_tensor(:,:,k) tracks the per-subcarrier residual, initialised to X_full.

R_tensor = X_full;   % M x N_in x K_s

% =========================================================================
%  SECTION 5: GREEDY DETECTION LOOP
% =========================================================================
% Implements Cui & Dai (2023) Algorithm 1, Steps 4-16.
% Outer loop: l = 1, ..., d_max  (fixed path count stopping rule).
% Per iteration:
%   (a) Compute BPD power map T(na,nd) using bilinear index gather.
%   (b) Pick peak atom (na*, nd*).
%   (c) Estimate complex gain for picked atom at each subcarrier.
%   (d) Subtract estimated component from R_tensor.

na_picked   = zeros(d_max, 1);
nd_picked   = zeros(d_max, 1);
T_peak_vec  = zeros(d_max, 1);
res_norm    = zeros(d_max, 1);
T_power_all = zeros(d_max, G_theta, G_r);
gain_store  = zeros(d_max, K_s);   % complex, d_max x K_s

for l = 1:d_max

    % ------------------------------------------------------------------
    %  (a) BPD power accumulation map T(na, nd)
    % ------------------------------------------------------------------
    % T(na,nd) = sum_{k=1}^{K_s} || psi_{j_k(na,nd)}^H * R_k ||_F^2
    % where j_k(na,nd) = (Lambda(nd,k)-1)*G_theta + Gamma(na,k).
    %
    % Per subcarrier k:
    %   1. Form atom power map P_mat_k(na,nd) = power of carrier-freq
    %      atom (na,nd) acting on R_k, as G_theta x G_r matrix.
    %   2. Gather: T(na,nd) += P_mat_k(Gamma(na,k), Lambda(nd,k)).
    %
    % The gather step shifts attention from the carrier-freq atom (na,nd)
    % to the bilinear-predicted atom position at subcarrier k -- this is
    % the BPD contribution over plain SOMP.

    T = zeros(G_theta, G_r);

    for k = 1:K_s
        Rk = R_tensor(:, :, k);                    % M x N_in

        % Batched matched filter: (G_theta*G_r) x N_in
        ProjRk = Psi_mat' * Rk;

        % Per-atom squared Frobenius norm: (G_theta*G_r) x 1
        atom_power_k = sum(abs(ProjRk).^2, 2);

        % Reshape to G_theta x G_r (angle-major)
        P_mat_k = reshape(atom_power_k, G_theta, G_r);

        % Bilinear gather:
        %   T(na,nd) += P_mat_k(Gamma(na,k), Lambda(nd,k))
        % Linear index: lin(na,nd) = (Lambda(nd,k)-1)*G_theta + Gamma(na,k)
        % Gamma(:,k) is G_theta x 1; Lambda(:,k) is G_r x 1.
        Gamma_k  = Gamma(:, k);    % G_theta x 1
        Lambda_k = Lambda(:, k);   % G_r x 1

        % Broadcast to G_theta x G_r index matrix
        lin_idx = (Lambda_k.' - 1) * G_theta + Gamma_k;  % G_theta x G_r

        T = T + P_mat_k(lin_idx);
    end

    T_power_all(l, :, :) = T;

    % ------------------------------------------------------------------
    %  (b) Pick peak atom
    % ------------------------------------------------------------------
    [T_peak_val, flat_idx] = max(T(:));
    [na_star, nd_star]     = ind2sub([G_theta, G_r], flat_idx);

    na_picked(l)  = na_star;
    nd_picked(l)  = nd_star;
    T_peak_vec(l) = T_peak_val;

    % ------------------------------------------------------------------
    %  (c)+(d) Gain estimation and residual update (combined, one k-loop)
    % ------------------------------------------------------------------
    % For stochastic sources s_{l,k}(n) ~ CN(0,p_l), each snapshot n has
    % a DIFFERENT realisation.  The correct OMP residual update is the
    % per-snapshot orthogonal projection onto span{a_lk}^perp:
    %
    %   R_tensor(:,:,k) -= a_lk * (a_lk' * R_tensor(:,:,k))
    %
    % where (a_lk' * R_tensor(:,:,k)) is 1 x N_in -- a different scalar
    % per snapshot.  This correctly removes the picked path from every
    % snapshot column of the residual.
    %
    % BUG FIX (Day-1 v2): the previous code subtracted
    %   a_lk * gain_scalar * ones(1,N_in)
    % where gain_scalar = a_lk' * mean(R,2) was a single scalar applied
    % uniformly to all snapshots.  For stochastic sources this left a
    % large residual and caused the greedy loop to re-pick the same atom
    % on every iteration without updating.

    theta_star = theta_grid(na_star);
    u_star     = u_grid(nd_star);

    for k = 1:K_s
        a_lk   = wb_nf_fresnel_steer(theta_star, u_star, alpha_k(k), params);
        proj_k = a_lk' * R_tensor(:, :, k);          % 1 x N_in per-snapshot coeff
        gain_store(l, k) = mean(proj_k);              % scalar mean for output
        R_tensor(:, :, k) = R_tensor(:, :, k) - a_lk * proj_k;  % project out
    end

    res_norm(l) = sqrt(sum(abs(R_tensor(:)).^2));

end  % greedy loop

% =========================================================================
%  SECTION 5b: JOINT LS POLISHING  (runs if bpd_do_ls_polish == true)
% =========================================================================
% After the greedy loop fixes the support {(na_picked(l), nd_picked(l))},
% re-estimate all d_max path gains jointly per subcarrier using a single
% overdetermined LS solve.  The support is NOT changed; only gains are
% updated.
%
% For each subcarrier k, build the support matrix:
%   A_support_k  in C^{M x d_max}  (unnormalised atoms, each has norm sqrt(M))
% Then solve the overdetermined system:
%   A_support_k * s_hat_k = X_full(:,:,k)   (M x N_in RHS)
% MATLAB backslash A_support_k \ X_full(:,:,k) calls LAPACK dgels and
% returns the d_max x N_in matrix of LS solutions -- one solution vector
% per snapshot column simultaneously.
%
% The per-subcarrier mean gain across snapshots replaces gain_store in the
% alpha_hat computation in Section 6.
%
% res_norm_polished reports the Frobenius norm of the full-support residual
% after all d_max paths are simultaneously removed by the LS solution.

if do_ls_polish && (d_max >= 1)

    gain_polished = zeros(d_max, K_s);   % complex, replaces gain_store in Sec 6

    % Build support matrices once per subcarrier (cache for residual calc).
    % A_all_k{k} = M x d_max unnormalised support matrix at subcarrier k.
    A_all_k = cell(K_s, 1);
    for k = 1:K_s
        Ak = zeros(M, d_max);
        for l = 1:d_max
            % Unnormalised atom: unit-norm * sqrt(M) so columns have power M.
            Ak(:, l) = wb_nf_fresnel_steer( ...
                theta_grid(na_picked(l)), u_grid(nd_picked(l)), ...
                alpha_k(k), params) * sqrt(M);
        end
        A_all_k{k} = Ak;
    end

    % LS solve and gain extraction.
    for k = 1:K_s
        % Overdetermined system: (M x d_max) \ (M x N_in) --> d_max x N_in.
        % Each column of X_full(:,:,k) is an independent LS problem;
        % MATLAB handles all N_in columns in one LAPACK call.
        s_hat_k = A_all_k{k} \ X_full(:, :, k);   % d_max x N_in
        gain_polished(:, k) = mean(s_hat_k, 2);    % d_max x 1 mean over snapshots
    end

    % Compute polished full-support residual.
    % R_polished_k = X_full(:,:,k) - A_support_k * s_hat_k  for all k,
    % then total Frobenius norm.
    res_polished_sq = 0;
    for k = 1:K_s
        s_hat_k        = A_all_k{k} \ X_full(:, :, k);   % d_max x N_in
        R_pol_k        = X_full(:, :, k) - A_all_k{k} * s_hat_k;
        res_polished_sq = res_polished_sq + sum(abs(R_pol_k(:)).^2);
    end
    res_norm_polished = sqrt(res_polished_sq);

else
    % LS polish disabled -- fall back to greedy gains and report NaN.
    gain_polished     = gain_store;
    res_norm_polished = NaN;
end

% =========================================================================
%  SECTION 6: PHYSICAL PARAMETER CONVERSION AND GAIN REFERENCING
% =========================================================================

theta_hat = theta_grid(na_picked);   % d_max x 1  [rad]
r_hat     = r_grid(nd_picked);       % d_max x 1  [m]

% Carrier-referenced gain: remove alpha_k scaling from each subcarrier.
% Uses gain_polished when LS polish is enabled, gain_store otherwise.
alpha_hat = zeros(d_max, 1);
for l = 1:d_max
    alpha_hat(l) = mean(gain_polished(l, :).' ./ alpha_k);
end

diag_info.T_power          = T_power_all;
diag_info.res_norm         = res_norm;
diag_info.res_norm_polished = res_norm_polished;
diag_info.na_picked        = na_picked;
diag_info.nd_picked        = nd_picked;
diag_info.T_peak           = T_peak_vec;
diag_info.theta_grid       = theta_grid;
diag_info.r_grid           = r_grid;
diag_info.alpha_k_vec      = alpha_k;

% =========================================================================
%  CSV DIAGNOSTIC BLOCK
% =========================================================================
if do_csv
    fid = fopen(csv_path, 'w');
    if fid == -1
        warning('bpd_baseline: could not open %s for writing', csv_path);
    else
        fprintf(fid, ...
            'iteration,theta_idx,r_idx,theta_hat_deg,r_hat_m,T_peak,res_norm,res_norm_polished\n');
        for l = 1:d_max
            fprintf(fid, '%d,%d,%d,%.6f,%.6f,%.6e,%.6e,%.6e\n', ...
                l, na_picked(l), nd_picked(l), ...
                theta_hat(l)*180/pi, r_hat(l), ...
                T_peak_vec(l), res_norm(l), res_norm_polished);
        end
        fclose(fid);
        fprintf('bpd_baseline: diagnostics written to %s\n', csv_path);
    end
end

end  % end bpd_baseline
