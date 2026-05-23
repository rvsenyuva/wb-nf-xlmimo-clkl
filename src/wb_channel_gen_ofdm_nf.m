function [X_full, Y_full, H_true, theta_true, r_true, p_true, N0, W_comb] = ...
        wb_channel_gen_ofdm_nf(P, SNR_dB, use_exact)
%WB_CHANNEL_GEN_OFDM_NF  Wideband near-field OFDM channel generator.
%
%  Paper C Phase 2 shared utility.  Generates per-subcarrier snapshot
%  tensors for all Paper C estimators.  Wideband extension of Paper B's
%  nf_gen_channel.m via alpha_k = f_k/fc frequency scaling.
%
%  CRITICAL RULES
%  --------------
%  Rule 1: geometry (theta_l, r_l) drawn ONCE, outside k-loop.
%          Enforces 3GPP TR 38.901 Sec. 7.6.5 cross-frequency consistency.
%          Violating this destroys the cross-frequency score-function
%          structure underlying CRB Proposition 1 (data diversity).
%  Rule 2: path gains s_{l,k}(n) ~ CN(0,p_l) drawn fresh per subcarrier k.
%          Enforces Assumption W1 cross-frequency independence.
%          Violating this introduces spurious gain correlation that
%          invalidates FIM additivity J_WB = sum_k J_k.
%
%  SIGNAL MODEL (Paper C Sec. II, Phase 1 locked notation)
%  -------------------------------------------------------
%  Per-subcarrier full-array snapshot (M x N at subcarrier k):
%    X_full(:,:,k) = A_k * S_k + W_k
%  where
%    A_k(:,l) = a_{l,k}  -- wideband Fresnel/USW steering vector
%    S_k(l,n) = s_{l,k}(n) ~ CN(0,p_l)  (fresh per k, Rule 2)
%    W_k      ~ CN(0, N0*I_M)
%
%  Wideband Fresnel steering vector (GLOBECOM eq. a_k):
%    [a_{l,k}]_m = exp(j*alpha_k*(m_bar*omega_l + m_bar^2*kappa_l))
%    omega_l = (2*pi*d_ant/lambda_c)*cos(theta_l)  [POSITIVE cosine]
%    kappa_l = (pi*d_ant^2/lambda_c)*sin(theta_l)^2 / r_l
%    alpha_k = f_k/fc   (= P.alpha_k_vec(k))
%
%  Implemented by calling wb_nf_fresnel_steer (unit-norm) * sqrt(M).
%
%  Wideband USW truth model (use_exact = true, Branch B locked):
%    Phase-only; no element-amplitude factor (1/r_m suppressed).
%    Wideband form substitutes lambda_c -> lambda_c/alpha_k in
%    nf_usw_steer.m (i.e., lambda_eff = lambda_c/alpha_k).
%
%  SNR CONVENTION
%  --------------
%    sig_pow = mean over k of trace(A_k * diag(p_true) * A_k')
%    N0      = sig_pow / (M * 10^(SNR_dB/10))
%  N0 is per-element, per-snapshot, per-subcarrier noise variance.
%  Total noise power = M * K_s * N * N0.
%  Matches bpd_toy_test.m toy_scene() helper exactly.
%
%  COMBINER
%  --------
%  W_comb (M x N_RF) is random constant-modulus, frequency-flat (same
%  for all subcarriers, per Paper C Sec. II).  Built ONCE in Phase B.
%  Y_full(:,:,k) = W_comb' * X_full(:,:,k)  (no whitening here;
%  whitening is the estimator's responsibility).
%
%  INPUTS
%  ------
%  P         : parameter struct.  Required fields (Paper B + wideband):
%                Paper B core: M, N, N_RF, d, lambda, lambda_c, d_ant,
%                  theta_lo, theta_hi, r_lo_fac, r_hi_fac, r_RD
%                Wideband additions: K, K_s, Delta_f, B,
%                  alpha_k_vec (K_s x 1), k_indices (K_s x 1)
%                Optional: wb_gen_write_csv (default true),
%                           wb_gen_csv_path (default 'wb_channel_gen_diag.csv')
%  SNR_dB    : scalar per-subcarrier SNR [dB]
%  use_exact : logical (default true)
%                true  -- phase-only USW truth model (Branch B locked)
%                false -- Fresnel approximation
%
%  OUTPUTS
%  -------
%  X_full     : M x N x K_s  -- full-array snapshot tensor
%  Y_full     : N_RF x N x K_s  -- compressed snapshots (W_comb' * X_full)
%  H_true     : M x d x K_s  -- per-subcarrier steering matrices
%                 H_true(:,:,k) has columns a_{l,k} (unnormalised, norm=sqrt(M))
%  theta_true : d x 1 [rad]  -- true angles (same for all k, Rule 1)
%  r_true     : d x 1 [m]    -- true ranges (same for all k, Rule 1)
%  p_true     : d x 1        -- true path powers (= 1/d, equal)
%  N0         : scalar       -- noise variance per element per subcarrier
%  W_comb     : M x N_RF     -- hybrid combiner (return so downstream
%                               estimators use the identical W that
%                               generated Y_full)
%
%  CALLED BY
%  ---------
%    wb_channel_gen_toy_test.m  -- Task 11.2 toy verification script
%    wb_clkl_driver.m           -- Task 11.4
%    run_monte_carlo_paperC.m   -- Task 11.8
%
%  DEPENDENCIES
%  ------------
%    wb_nf_fresnel_steer.m      -- Task 11.1, shared Phase 2 steering utility
%    nf_usw_steer.m             -- Paper B exact USW (use_exact = true path)
%
%  Author : R. V. Senyuva (Maltepe University)
%  Date   : May 2026

% =========================================================================
%  PHASE A -- Parameter validation
% =========================================================================

if nargin < 3; use_exact = true; end

assert(P.theta_lo >= 15*pi/180, ...
    ['wb_channel_gen_ofdm_nf: theta_lo = %.1f deg is below 15 deg.\n' ...
     'Inherited from nf_gen_channel.m Step 1 audit: theta_lo >= 15 deg\n' ...
     'required for identifiable near-field curvature. Use 20 deg.'], ...
    P.theta_lo * 180/pi);
assert(P.theta_lo < P.theta_hi, ...
    'wb_channel_gen_ofdm_nf: theta_lo must be strictly less than theta_hi.');
assert(isfield(P,'K_s') && isfield(P,'alpha_k_vec'), ...
    'wb_channel_gen_ofdm_nf: P must contain wideband fields K_s and alpha_k_vec.');
assert(P.K_s <= P.K, ...
    'wb_channel_gen_ofdm_nf: K_s must be <= K (used <= total subcarriers).');
assert(numel(P.alpha_k_vec) == P.K_s, ...
    'wb_channel_gen_ofdm_nf: numel(alpha_k_vec) must equal K_s.');
assert(isfield(P,'lambda_c'), ...
    'wb_channel_gen_ofdm_nf: P.lambda_c required (= c0/fc). Add to param struct.');

% Unpack scalars used in multiple phases
M        = P.M;
N        = P.N;
N_RF     = P.N_RF;
d        = P.d;
K_s      = P.K_s;
alpha_k  = P.alpha_k_vec;   % K_s x 1

% Rayleigh distance and range bounds (inherited from nf_gen_channel.m)
% r_RD may already be in P (set by nf_update_derived_pub); recompute
% defensively to ensure consistency with current P.M and P.lambda.
D        = (M - 1) * P.d_ant;
r_RD     = 2 * D^2 / P.lambda_c;   % use lambda_c (Paper C name)
r_min    = P.r_lo_fac * r_RD;
r_max    = P.r_hi_fac * r_RD;

% CSV diagnostic defaults
if ~isfield(P,'wb_gen_write_csv');  P.wb_gen_write_csv = true;          end
if ~isfield(P,'wb_gen_csv_path');   P.wb_gen_csv_path  = ...
        'wb_channel_gen_diag.csv';                                       end

% =========================================================================
%  PHASE B -- Geometry draw (ONCE, outside k-loop) [RULE 1]
%             + combiner build (frequency-flat)
% =========================================================================

% --- Random geometry: uniform on [theta_lo, theta_hi] and [r_min, r_max]
% These draws happen exactly ONCE per call.  Every downstream estimator
% that calls this function in a Monte Carlo loop gets an independent
% realization but sees the same geometry across all K_s subcarriers.
theta_true = P.theta_lo + (P.theta_hi - P.theta_lo) * rand(d, 1);  % d x 1 [rad]
r_true     = r_min      + (r_max - r_min)            * rand(d, 1);  % d x 1 [m]
p_true     = ones(d, 1) / d;                                        % equal power

% --- Frequency-flat constant-modulus combiner (Paper C Sec. II)
% [W_comb]_{m,j} = (1/sqrt(M)) * exp(j*phi_{mj}),  phi ~ U[0, 2*pi)
% Same W for all subcarriers -- this is the frequency-flat assumption.
% Returned as 8th output so downstream estimators use the identical W.
phi_mat = 2*pi * rand(M, N_RF);
W_comb  = (1/sqrt(M)) * exp(1j * phi_mat);   % M x N_RF

% =========================================================================
%  PHASE C -- Per-subcarrier loop
%             Step C1: compute average signal power for N0
%             Step C2: generate snapshots
% =========================================================================
% Two-pass design: C1 accumulates signal power, C2 generates data.
% Cost is 2*K_s steering-vector builds per path -- acceptable for K_s<=512.

% --- C1: Signal power (averaged over K_s subcarriers) for N0 definition
sig_pow = 0;
for k = 1:K_s
    Ak       = build_A_k(theta_true, r_true, alpha_k(k), use_exact, M, d, P);
    sig_pow  = sig_pow + real(trace(Ak * diag(p_true) * Ak'));
end
sig_pow = sig_pow / K_s;
N0      = sig_pow / (M * 10^(SNR_dB / 10));

% --- C2: Snapshot generation
X_full = zeros(M,   N,   K_s, 'like', 1+1j);
Y_full = zeros(N_RF, N,  K_s, 'like', 1+1j);
H_true = zeros(M,   d,   K_s, 'like', 1+1j);

for k = 1:K_s
    % Build per-subcarrier steering matrix A_k (M x d)
    A_k = build_A_k(theta_true, r_true, alpha_k(k), use_exact, M, d, P);

    % Draw path gains S_k (d x N): FRESH per subcarrier k [RULE 2]
    % s_{l,k}(n) ~ CN(0, p_l)  independently for each k
    S_k = zeros(d, N);
    for ell = 1:d
        S_k(ell, :) = sqrt(p_true(ell) / 2) * ...
            (randn(1, N) + 1j * randn(1, N));
    end

    % Additive noise W_k (M x N): CN(0, N0*I_M)
    W_k = sqrt(N0 / 2) * (randn(M, N) + 1j * randn(M, N));

    % Assemble snapshots and compress
    X_full(:, :, k) = A_k * S_k + W_k;
    Y_full(:, :, k) = W_comb' * X_full(:, :, k);
    H_true(:, :, k) = A_k;
end

% =========================================================================
%  PHASE D -- Diagnostics / CSV output
% =========================================================================
if P.wb_gen_write_csv
    % 6-column CSV: k_idx, alpha_k, theta_true_deg (path 1), r_true_m (path 1),
    %               sig_power_k, N0
    % For multi-path (d>1) reports path-1 geometry only; full geometry is
    % in the returned theta_true/r_true vectors.
    fid = fopen(P.wb_gen_csv_path, 'w');
    fprintf(fid, 'k_idx,alpha_k,theta_hat_deg,r_hat_m,sig_power_k,N0\n');
    for k = 1:K_s
        Ak_k     = H_true(:, :, k);
        spow_k   = real(trace(Ak_k * diag(p_true) * Ak_k'));
        fprintf(fid, '%d,%.10f,%.6f,%.6f,%.10f,%.10e\n', ...
            k, alpha_k(k), theta_true(1)*180/pi, r_true(1), spow_k, N0);
    end
    fclose(fid);
end

end  % end main function


% =========================================================================
%  LOCAL FUNCTION: build_A_k
%  Build the M x d steering matrix at subcarrier k.
%  Called twice per k (power pass + snapshot pass) to avoid storing A_k.
% =========================================================================
function A_k = build_A_k(theta_true, r_true, alpha_k_val, use_exact, M, d, P)
%BUILD_A_K  Per-subcarrier steering matrix (M x d).
%
%  Columns are unnormalised atoms: norm(A_k(:,l)) = sqrt(M).
%
%  use_exact = true  : wideband phase-only USW
%                        lambda_eff = lambda_c / alpha_k_val
%                        [a]_m = exp(-j*(2*pi/lambda_eff)*(dist_m - r))
%  use_exact = false : wideband Fresnel (wb_nf_fresnel_steer * sqrt(M))

A_k = zeros(M, d);
for ell = 1:d
    if use_exact
        % Wideband USW: substitute lambda_c -> lambda_c/alpha_k in nf_usw_steer.
        % Create a local copy of P with P_local.lambda = lambda_eff.
        % nf_usw_steer reads only P.lambda, P.M, P.d_ant.
        % P.lambda_c is NOT overridden; only P.lambda changes for this call.
        P_local        = P;
        P_local.lambda = P.lambda_c / alpha_k_val;   % lambda_eff at subcarrier k
        % nf_usw_steer returns unnormalised: ||a||^2 = M (Branch B: phase-only)
        A_k(:, ell) = nf_usw_steer(theta_true(ell), r_true(ell), P_local);
    else
        % Fresnel: wb_nf_fresnel_steer returns unit-norm; multiply by sqrt(M)
        A_k(:, ell) = wb_nf_fresnel_steer( ...
            theta_true(ell), 1/r_true(ell), alpha_k_val, P) * sqrt(M);
    end
end
end
