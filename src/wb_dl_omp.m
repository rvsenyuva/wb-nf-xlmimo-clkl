function [theta_hat, r_hat, info] = wb_dl_omp(X_full, P)
%WB_DL_OMP  Wideband Dictionary-Learning OMP (DL-OMP) near-field estimator.
%
%  [theta_hat, r_hat, info] = wb_dl_omp(X_full, P)
%
%  Wideband extension of the narrowband DL-OMP estimator from Paper B
%  (nf_zhang.m).  Implements Zhang, Zhang & Eldar (IEEE Trans. Commun. 2024)
%  extended to K_s OFDM subcarriers using the Paper C locked wideband
%  steering vector model:
%
%    [a_{l,k}]_m = exp(j*alpha_k*omega_l*m_bar - j*alpha_k*kappa_l*m_bar^2)
%
%  where omega_l = (2*pi*d_ant/lambda_c)*cos(theta_l),
%        kappa_l = (pi*d_ant^2/lambda_c)*sin^2(theta_l)/r_l,
%        alpha_k = f_k/f_c  (per-subcarrier frequency ratio).
%
%  This is Baseline B5 in the Paper C six-baseline comparison.
%  DL-OMP is a FULL-ARRAY method: it receives the uncompressed snapshot
%  tensor X_full (M x N x K_s).  No hybrid combiner W_comb is used.
%
%  WIDEBAND EXTENSIONS vs. nf_zhang (narrowband):
%  -----------------------------------------------
%  1. Subarray covariances: K_s cell arrays R1_cell{k}, R2_cell{k} instead
%     of single matrices R1_full, R2_full.
%  2. OMP score: cross-subcarrier sum  sum_k real(w_{s,k,j}^H z_j_cell{k} w_{s,k,j})
%     instead of single-covariance score.
%  3. Dictionary correction: alpha_k factor in correction exponent.
%  4. Deflation: per-subcarrier z1_cell{k}, z2_cell{k} updated independently.
%
%  FREQUENCY-INDEPENDENT BLOCKS (transferred verbatim from nf_zhang):
%  - Subarray geometry (M_s, idx1, idx2, delta)
%  - Angle grid (theta_grid, m_s, m_s2, omega_g)
%  - Base dictionary at r_max (W0, carrier freq alpha_k=1)
%  - Law-of-sines range formula (geometric, alpha_k-independent)
%
%  INPUTS
%  ------
%  X_full : M x N x K_s  full-array wideband snapshot tensor
%           (1st output of wb_channel_gen_ofdm_nf)
%  P      : parameter struct with fields:
%             .M          -- ULA element count
%             .N          -- snapshots per subcarrier
%             .d          -- number of paths to estimate
%             .lambda_c   -- carrier wavelength [m]
%             .d_ant      -- element spacing [m]
%             .r_RD       -- Rayleigh distance [m]
%             .r_lo_fac   -- r_min = r_lo_fac * r_RD
%             .r_hi_fac   -- r_max = r_hi_fac * r_RD
%             .theta_lo   -- minimum angle [rad]
%             .theta_hi   -- maximum angle [rad]
%             .Q_theta    -- angle dictionary size (default 256)
%             .K_s        -- number of active subcarriers
%             .alpha_k_vec-- K_s x 1 vector of alpha_k = f_k/f_c values
%
%  OUTPUTS
%  -------
%  theta_hat : d x 1  estimated angles [rad]
%  r_hat     : d x 1  estimated ranges [m]
%  info      : struct with diagnostic fields:
%                .K_iter_used  -- K_iter constant used (always 3)
%                .d_hat        -- number of paths returned (= P.d)
%                .delta_m      -- subarray centroid separation [m]
%                .r_los        -- raw law-of-sines range for last path [m]
%                               (before clamping; diagnostic only)
%                .r_per_path   -- d x 1 estimated ranges per path [m]
%                .theta_per_path -- d x 1 estimated angles per path [rad]
%
%  NOTE ON K_iter:
%  K_iter = 3 is a local constant (Paper B default, sufficient per Zhang
%  Fig. 9). It is NOT a P struct field; the caller does not set it.
%
%  Author  : R. V. Senyuva (Maltepe University)
%  Date    : May 2026
%  Ref     : Zhang, Zhang & Eldar, IEEE Trans. Commun. 2024 (Alg. 1-2)
%            Paper C Phase 2, Task 11.7

% ---- Extract parameters ------------------------------------------------
M        = P.M;
N        = P.N;
d        = P.d;
lambda_c = P.lambda_c;
d_ant    = P.d_ant;
K_s      = P.K_s;
K_iter   = 3;              % dictionary update iterations (local constant)

r_min = P.r_lo_fac * P.r_RD;
r_max = P.r_hi_fac * P.r_RD;
u_max = 1 / r_max;         % inverse range at r_max (dictionary init point)

% ---- Subarray geometry (verbatim from nf_zhang lines 58-68) -----------
% Only change: lambda -> lambda_c (notation lock)
M_s       = floor(M / 2);
idx1      = 1:M_s;
idx2      = (M - M_s + 1):M;

cent1_pos = mean(idx1 - 1) * d_ant;
cent2_pos = mean(idx2 - 1) * d_ant;
delta     = cent2_pos - cent1_pos;   % subarray separation [m]

% ---- Angle grid and initial Fresnel dictionary at r_max ---------------
% Verbatim from nf_zhang lines 70-80; lambda -> lambda_c.
% W0 is built ONCE at carrier freq (alpha_k = 1); corrections applied later.
Q_th       = P.Q_theta;
theta_grid = linspace(P.theta_lo, P.theta_hi, Q_th).';  % Q_th x 1
m_s        = ((0:M_s-1).' - (M_s-1)/2);                  % M_s x 1
m_s2       = m_s .^ 2;                                   % M_s x 1
omega_g    = (2*pi*d_ant/lambda_c) * cos(theta_grid).';  % 1 x Q_th (carrier)

% Initial Fresnel dictionary at r_max, carrier freq (alpha_k = 1)
kappa_max = (pi*d_ant^2/lambda_c) * sin(theta_grid).^2 / r_max; % Q_th x 1
PHI_max   = m_s * omega_g - m_s2 * kappa_max.';          % M_s x Q_th
W0        = exp(1j * PHI_max);                            % M_s x Q_th

% Precompute carrier-frequency curvature scale per grid atom (for correction)
% c_q(s) = (pi*d_ant^2/lambda_c) * sin(theta_grid(s))^2
c_q = (pi*d_ant^2/lambda_c) * sin(theta_grid).^2;        % Q_th x 1

% ---- Wideband: per-subcarrier subarray covariances --------------------
% R1_cell{k}, R2_cell{k}: M_s x M_s covariance matrices at each subcarrier.
% z1_cell{k}, z2_cell{k}: working residuals (initialised to R1/R2_cell{k}).
R1_cell = cell(K_s, 1);
R2_cell = cell(K_s, 1);
z1_cell = cell(K_s, 1);
z2_cell = cell(K_s, 1);

for k = 1:K_s
    X1_k         = X_full(idx1, :, k);                % M_s x N
    X2_k         = X_full(idx2, :, k);
    R1_k         = (1/N) * (X1_k * X1_k');            % M_s x M_s
    R2_k         = (1/N) * (X2_k * X2_k');
    R1_cell{k}   = (R1_k + R1_k') / 2;               % Hermitian-symmetrise
    R2_cell{k}   = (R2_k + R2_k') / 2;
    z1_cell{k}   = R1_cell{k};
    z2_cell{k}   = R2_cell{k};
end

% ---- Multi-path OMP loop ----------------------------------------------
active_theta = zeros(d, 1);
active_r     = zeros(d, 1);
r_los_last   = NaN;           % raw law-of-sines for last path (diagnostic)

for l = 1:d

    % Per-path working dictionaries (reset to W0 each path, then updated)
    W1 = W0;   % M_s x Q_th  (subarray 1 dictionary, carrier freq init)
    W2 = W0;   % M_s x Q_th  (subarray 2 dictionary, carrier freq init)

    % Initial estimates for K_iter loop
    theta1_hat  = theta_grid(1);
    theta2_hat  = theta_grid(end);
    r_hat_l     = 0.5 * (r_min + r_max);
    idx1_best   = 1;
    idx2_best   = Q_th;

    for ki = 1:K_iter

        % ---- Step 1: wideband OMP score for subarray 1 ----------------
        % score_1(s) = sum_k real( w_{s,k,1}^H * z1_cell{k} * w_{s,k,1} )
        % Per-subcarrier dict W1_k is W1 with alpha_k correction on
        % selected column; unselected columns remain at W0 (Option B).
        score_1 = zeros(Q_th, 1);
        score_2 = zeros(Q_th, 1);

        for k = 1:K_s
            alpha_k = P.alpha_k_vec(k);

            % Build W1_k: copy W1, then replace selected col with alpha_k version
            W1_k = W1;                              % M_s x Q_th (copy)
            if idx1_best >= 1
                % Correction for previously selected atom at subcarrier k
                corr1_k = exp(-1j * m_s2 * alpha_k * c_q(idx1_best) * ...
                              (1/r_hat_l - u_max));
                W1_k(:, idx1_best) = corr1_k .* W0(:, idx1_best);
            end

            % Build W2_k similarly for subarray 2
            W2_k = W2;
            if idx2_best >= 1
                corr2_k = exp(-1j * m_s2 * alpha_k * c_q(idx2_best) * ...
                              (1/r_hat_l - u_max));
                W2_k(:, idx2_best) = corr2_k .* W0(:, idx2_best);
            end

            % Accumulate scores over subcarriers
            score_1 = score_1 + real(sum(conj(W1_k) .* (z1_cell{k} * W1_k), 1)).';
            score_2 = score_2 + real(sum(conj(W2_k) .* (z2_cell{k} * W2_k), 1)).';
        end

        % ---- Step 2: select best atoms ---------------------------------
        [~, idx1_best] = max(score_1);
        theta1_hat     = theta_grid(idx1_best);

        [~, idx2_best] = max(score_2);
        theta2_hat     = theta_grid(idx2_best);

        % Avoid degenerate case (both subarrays select same angle)
        if abs(theta2_hat - theta1_hat) < 1e-4
            sc2_tmp = score_2;
            nw = max(1, round(Q_th / 30));
            lo = max(1, idx2_best - nw);
            hi = min(Q_th, idx2_best + nw);
            sc2_tmp(lo:hi) = -Inf;
            [~, idx2_best] = max(sc2_tmp);
            theta2_hat     = theta_grid(idx2_best);
        end

        % ---- Step 3: range by law of sines (verbatim from nf_zhang) ---
        % r1/sin(theta2) = delta/sin(theta2 - theta1)  -- eq. (18) Zhang 2024
        denom_r = sin(theta2_hat) - sin(theta1_hat);
        if abs(denom_r) < 1e-6
            r_hat_l = 0.5 * (r_min + r_max);
        else
            % Primary law-of-sines formula (nf_zhang lines 141-144)
            dsin = sin(theta2_hat - theta1_hat);
            if abs(dsin) > 1e-6
                r_hat_l = delta * sin(theta2_hat) / dsin;
            end
        end
        r_los_raw = r_hat_l;                          % save pre-clamp value
        r_hat_l   = min(r_max, max(r_min, abs(r_hat_l)));

        % ---- Step 4: dictionary update (per-path, carrier-freq base) --
        % Only update the selected columns of W1, W2 (Option B from spec).
        % The correction shifts each atom from r_max to r_hat_l at carrier.
        % Per-subcarrier alpha_k correction is applied in score step above.
        u_hat = 1 / r_hat_l;

        correction1 = exp(-1j * m_s2 * c_q(idx1_best) * (u_hat - u_max));
        W1(:, idx1_best) = correction1 .* W0(:, idx1_best);

        correction2 = exp(-1j * m_s2 * c_q(idx2_best) * (u_hat - u_max));
        W2(:, idx2_best) = correction2 .* W0(:, idx2_best);

    end  % K_iter loop

    % Store results for this path
    active_theta(l) = (theta1_hat + theta2_hat) / 2;
    active_r(l)     = r_hat_l;
    r_los_last      = r_los_raw;

    % ---- Per-subcarrier orthogonal-projection deflation ---------------
    % Deflate z1_cell{k} and z2_cell{k} using alpha_k-correct atoms.
    for k = 1:K_s
        alpha_k = P.alpha_k_vec(k);

        % Final subarray-1 atom at subcarrier k
        c1       = (pi*d_ant^2/lambda_c) * sin(theta1_hat)^2;
        kappa1_k = alpha_k * c1 * u_hat;
        omega1_k = alpha_k * (2*pi*d_ant/lambda_c) * cos(theta1_hat);
        a1_lk    = exp(1j * (m_s * omega1_k - m_s2 * kappa1_k));
        a1_lk    = a1_lk / norm(a1_lk);

        % Final subarray-2 atom at subcarrier k
        c2       = (pi*d_ant^2/lambda_c) * sin(theta2_hat)^2;
        kappa2_k = alpha_k * c2 * u_hat;
        omega2_k = alpha_k * (2*pi*d_ant/lambda_c) * cos(theta2_hat);
        a2_lk    = exp(1j * (m_s * omega2_k - m_s2 * kappa2_k));
        a2_lk    = a2_lk / norm(a2_lk);

        % Deflate residual covariances
        P1_k    = a1_lk * pinv(a1_lk);
        Perp1_k = eye(M_s) - P1_k;
        z1_cell{k} = Perp1_k * R1_cell{k} * Perp1_k';
        z1_cell{k} = (z1_cell{k} + z1_cell{k}') / 2;

        P2_k    = a2_lk * pinv(a2_lk);
        Perp2_k = eye(M_s) - P2_k;
        z2_cell{k} = Perp2_k * R2_cell{k} * Perp2_k';
        z2_cell{k} = (z2_cell{k} + z2_cell{k}') / 2;
    end

end  % path loop

% ---- Outputs -----------------------------------------------------------
theta_hat = active_theta;
r_hat     = active_r;

info.K_iter_used    = K_iter;
info.d_hat          = d;
info.delta_m        = delta;
info.r_los          = r_los_last;    % raw law-of-sines result (pre-clamp)
info.r_per_path     = active_r;
info.theta_per_path = active_theta;

end
