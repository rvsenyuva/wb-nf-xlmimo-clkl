function run_monte_carlo_paperC(P, sweep_type, sweep_vec, N_MC, out_dir)
%RUN_MONTE_CARLO_PAPERC  Master Monte Carlo driver for Paper C.
%
%  run_monte_carlo_paperC(P, sweep_type, sweep_vec, N_MC, out_dir)
%
%  Integrates all five Phase 2 estimators and the wideband compressed-
%  domain CRB baseline. Produces all performance data for Paper C Figs 8-10.
%
%  ESTIMATOR INVENTORY (six baselines)
%  ------------------------------------
%  B1  WB-BPD     bpd_baseline(X_full, P)             -- full-array
%  B2  WB-P-SOMP  wb_psomp(R_hat_cell, W_comb, P)     -- compressed
%  B3  WB-BF-SOMP (NOT IMPLEMENTED -- NaN placeholder)
%  B4  WB-CL-KL   wb_clkl_driver(R_hat_cell, W_comb,  -- compressed
%                      theta_init, r_init, p_init, P)  (BPD warm-start)
%  B5  WB-DL-OMP  wb_dl_omp(X_full, P)                -- full-array
%  B6  CRB        wb_crb_compressed(theta,r,p,N0,P)   -- deterministic
%
%  DATA-FLOW (locked, per realisation):
%    wb_channel_gen_ofdm_nf -> X_full, Y_full, W_comb
%    R_hat_cell{k} = (Yk*Yk')/N, Hermitian-symmetrised   (B2, B4)
%    B1 first (provides warm-start for B4).
%    B6 (CRB) called ONCE per sweep point, OUTSIDE the MC loop.
%
%  SWEEP TYPES
%  -----------
%  'snr'         -- Fig 8 RMSE vs. SNR (sweep_vec in dB)
%  'bandwidth'   -- Fig 9 NMSE vs. Bandwidth (sweep_vec in Hz)
%  'convergence' -- Fig 10 KL objective vs. iteration (sweep_vec = SNR dB)
%
%  INPUTS
%    P          : base parameter struct (Paper C signal model, Phase 1 locked)
%    sweep_type : 'snr' | 'bandwidth' | 'convergence'
%    sweep_vec  : vector of sweep values (SNR [dB] or B [Hz])
%    N_MC       : MC realisations per sweep point (200 publication; 50 dev)
%    out_dir    : output directory for CSV (created if absent)
%
%  OUTPUTS (written to out_dir)
%    mc_snr_sweep_<date>.csv
%    mc_bandwidth_sweep_<date>.csv
%    mc_convergence_sweep_<date>.csv
%    (one row per sweep-point x method)
%
%  PARFOR DISCIPLINE
%    parfor replaces for in the inner MC loop without modification.
%    rng is set per-realisation as rng(mc_seed+mc,'twister').
%    No shared mutable state inside the parfor body.
%
%  NEAR-RAYLEIGH DISCIPLINE (Lessons L33, L34)
%    All estimates including near-Rayleigh outliers are recorded.
%    Do NOT NaN-replace or filter any finite estimates.
%
%  Author   : R. V. Senyuva (Maltepe University)
%  Date     : May 2026
%  Ref      : Paper C Phase 2, Task 11.8.
%  Requires : bpd_baseline.m, wb_psomp.m, wb_clkl_driver.m,
%             wb_dl_omp.m, wb_crb_compressed.m, wb_channel_gen_ofdm_nf.m

% =========================================================================
%  ARGUMENT VALIDATION
% =========================================================================
narginchk(5, 5);
assert(ischar(sweep_type) || isstring(sweep_type), ...
    'run_monte_carlo_paperC: sweep_type must be a char/string.');
sweep_type = char(sweep_type);
valid_types = {'snr','bandwidth','convergence'};
assert(any(strcmp(sweep_type, valid_types)), ...
    'run_monte_carlo_paperC: sweep_type must be ''snr'', ''bandwidth'', or ''convergence''.');
assert(isvector(sweep_vec) && isnumeric(sweep_vec) && ~isempty(sweep_vec), ...
    'run_monte_carlo_paperC: sweep_vec must be a non-empty numeric vector.');
assert(isscalar(N_MC) && N_MC >= 1, ...
    'run_monte_carlo_paperC: N_MC must be a positive scalar integer.');
N_MC = round(N_MC);

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% =========================================================================
%  RIVIELLO SNR-AXIS PRE-PROCESSING (Session 4.9, D1)
%  When P.use_riviello_snr_axis is true AND sweep_type is 'snr':
%    - Draw N_MC per-UT SNR values from the 3GPP UMi path-loss /
%      shadow-fading model via riviello_snr_axis.m.
%    - Sort the draw to produce a CDF-ordered SNR axis for Fig. 11.
%    - Store the sorted vector in P.snr_per_trial_vec (N_MC x 1).
%    - Override sweep_vec to a single point [median SNR] so that the
%      outer loop runs exactly once and the CRB is evaluated at the
%      representative operating point.
%  For all other sweep types, or when the flag is false / absent,
%  this block is a no-op and all downstream code is unchanged.
% =========================================================================
if isfield(P, 'use_riviello_snr_axis') && P.use_riviello_snr_axis
    if ~strcmp(sweep_type, 'snr')
        warning(['run_monte_carlo_paperC: use_riviello_snr_axis=true is ' ...
            'only valid for sweep_type=''snr''. Flag ignored for ''%s''.'], ...
            sweep_type);
        P.use_riviello_snr_axis = false;
    else
        fprintf('[MC] Riviello mode: drawing %d per-UT SNR values...\n', N_MC);
        [snr_per_ut_raw, snr_info_riv] = riviello_snr_axis(P, N_MC);
        snr_axis_sorted = sort(snr_per_ut_raw);   % N_MC x 1, CDF order
        P.snr_per_trial_vec = snr_axis_sorted;    % broadcast into parfor
        P.snr_info_riviello = snr_info_riv;        % saved for reporting
        sweep_vec = median(snr_axis_sorted);       % single CRB eval point
        fprintf('[MC] Riviello SNR: median=%.2f dB  p10=%.2f dB  p90=%.2f dB\n', ...
            snr_info_riv.median_snr_dB, snr_info_riv.p10_snr_dB, ...
            snr_info_riv.p90_snr_dB);
        fprintf('[MC] Indoor: %.1f%%   LOS: %.1f%%\n', ...
            snr_info_riv.pct_indoor, snr_info_riv.pct_los);
    end
else
    P.use_riviello_snr_axis = false;   % normalise absent field
end

% =========================================================================
%  OUTPUT CSV PATH
% =========================================================================
date_str = datestr(now, 'yyyymmdd_HHMMSS');
csv_name = sprintf('mc_%s_sweep_%s.csv', sweep_type, date_str);
csv_path = fullfile(out_dir, csv_name);

% =========================================================================
%  FIXED BANDWIDTH SUBCARRIER SPACING (for bandwidth sweep)
% =========================================================================
Delta_f_fixed = 25e6;   % 25 MHz fixed subcarrier spacing (Lesson L29)

% =========================================================================
%  GLOBAL RNG SEED (once before outer loop; NOT inside parfor)
% =========================================================================
rng(42, 'twister');
mc_seed_base = 1000;   % per-realisation seed = mc_seed_base + mc

% =========================================================================
%  CSV HEADER
% =========================================================================
write_csv_header(csv_path, sweep_type);

fprintf('[MC] sweep_type=%s  N_MC=%d  n_points=%d\n', ...
    sweep_type, N_MC, numel(sweep_vec));

% =========================================================================
%  OUTER SWEEP LOOP
% =========================================================================
n_sweep = numel(sweep_vec);

for s = 1 : n_sweep

    % ------------------------------------------------------------------
    %  1. Update P for this sweep point
    % ------------------------------------------------------------------
    P_s = P;   % copy; do not mutate P across iterations

    if strcmp(sweep_type, 'snr')
        SNR_dB_s = sweep_vec(s);

    elseif strcmp(sweep_type, 'bandwidth')
        B_s        = sweep_vec(s);
        P_s.B      = B_s;
        P_s.K_s    = min(round(B_s / Delta_f_fixed), 512);
        P_s.Delta_f = B_s / P_s.K_s;
        k_idx      = (-(P_s.K_s/2) : (P_s.K_s/2 - 1)).';
        P_s.alpha_k_vec = 1 + k_idx * P_s.Delta_f / P_s.fc;
        P_s.k_indices   = k_idx;
        P_s.K           = P_s.K_s;
        SNR_dB_s        = P.snr_fixed;   % fixed SNR for bandwidth sweep

    elseif strcmp(sweep_type, 'convergence')
        SNR_dB_s = sweep_vec(s);
        % Fixed scene geometry is authoritative in P (set by setup_production_P_v4
        % via P.conv_theta and P.conv_r).  DO NOT reassign conv_r/conv_theta here
        % -- doing so overrides the setup value and breaks the single source of
        % truth.  Sprint A fix: prior hardcode of r=5.0 m here was the proximate
        % cause of the r_hi_fac inconsistency (5.0 m > r_hi=4.2525 m at 0.20).
        P_s.convergence_fixed_scene = true;
        % P.conv_theta and P.conv_r already set and Fresnel-verified by
        % setup_production_P_v4.  Print to confirm at runtime:
        fprintf('[MC] convergence scene: theta=%.1f deg, r=%.2f m\n', ...
            P_s.conv_theta * 180/pi, P_s.conv_r);
    end

    fprintf('[MC] sweep point %d/%d  (%s=%.4g)\n', s, n_sweep, sweep_type, sweep_vec(s));

    % ------------------------------------------------------------------
    %  2. CRB: one call per sweep point, OUTSIDE the MC loop (B6)
    % ------------------------------------------------------------------
    % Use nominal parameters: mid-range angle, mid-range range, equal power.
    theta_nom = (P_s.theta_lo + P_s.theta_hi) / 2;
    r_nom     = P_s.r_lo_fac * P_s.r_RD + ...
                (P_s.r_hi_fac - P_s.r_lo_fac) * P_s.r_RD / 2;
    p_nom     = ones(P_s.d, 1) / P_s.d;
    N0_nom    = 10^(-SNR_dB_s / 10);

    try
        [crb_r_s, crb_theta_s, ~] = wb_crb_compressed( ...
            repmat(theta_nom, P_s.d, 1), ...
            repmat(r_nom,     P_s.d, 1), ...
            p_nom, N0_nom, P_s);
        crb_r_val     = mean(crb_r_s);       % scalar [m],  sqrt-CRB
        crb_theta_val = mean(crb_theta_s);   % scalar [deg], sqrt-CRB
    catch ME
        warning('run_monte_carlo_paperC: CRB failed at sweep point %d: %s', ...
            s, ME.message);
        crb_r_val     = NaN;
        crb_theta_val = NaN;
    end

    % ------------------------------------------------------------------
    %  3. Per-realisation MC loop (parfor-ready)
    % ------------------------------------------------------------------
    % Pre-allocate result arrays (parfor-safe: no cell-of-cell).
    th_B1_all = nan(P_s.d, N_MC);
    r_B1_all  = nan(P_s.d, N_MC);
    th_B2_all = nan(P_s.d, N_MC);
    r_B2_all  = nan(P_s.d, N_MC);
    th_B4_all = nan(P_s.d, N_MC);
    r_B4_all  = nan(P_s.d, N_MC);
    th_B5_all = nan(P_s.d, N_MC);
    r_B5_all  = nan(P_s.d, N_MC);

    th_true_all = nan(P_s.d, N_MC);
    r_true_all  = nan(P_s.d, N_MC);

    % B4 convergence info (for Fig 10)
    L_hist_all  = cell(N_MC, 1);   % each: n_iter x 1 (variable length)
    b4_n_iter   = nan(1, N_MC);
    b4_converged = nan(1, N_MC);

    % Method runtimes [s]
    rt_B1 = nan(1, N_MC);
    rt_B2 = nan(1, N_MC);
    rt_B4 = nan(1, N_MC);
    rt_B5 = nan(1, N_MC);

    parfor mc = 1 : N_MC

        % Select per-trial SNR: Riviello CDF-ordered draw or fixed sweep point
        if P_s.use_riviello_snr_axis
            SNR_dB_mc = P_s.snr_per_trial_vec(mc);   % mc-th sorted UT SNR [dB]
        else
            SNR_dB_mc = SNR_dB_s;                     % fixed sweep-point SNR [dB]
        end

        % Per-realisation result struct (parfor accumulation)
        res = run_one_realisation(P_s, SNR_dB_mc, mc_seed_base + mc);

        th_B1_all(:, mc) = res.th_B1;
        r_B1_all(:, mc)  = res.r_B1;
        th_B2_all(:, mc) = res.th_B2;
        r_B2_all(:, mc)  = res.r_B2;
        th_B4_all(:, mc) = res.th_B4;
        r_B4_all(:, mc)  = res.r_B4;
        th_B5_all(:, mc) = res.th_B5;
        r_B5_all(:, mc)  = res.r_B5;

        th_true_all(:, mc) = res.theta_true;
        r_true_all(:, mc)  = res.r_true;

        L_hist_all{mc}   = res.L_hist;
        b4_n_iter(mc)    = res.n_iter;
        b4_converged(mc) = res.converged;

        rt_B1(mc) = res.rt_B1;
        rt_B2(mc) = res.rt_B2;
        rt_B4(mc) = res.rt_B4;
        rt_B5(mc) = res.rt_B5;

    end  % parfor mc

    % ------------------------------------------------------------------
    %  4. Aggregate metrics across N_MC realisations
    % ------------------------------------------------------------------

    % Hungarian-matched errors
    [err_th_B1, err_r_B1] = compute_matched_errors(th_B1_all, r_B1_all, ...
                                th_true_all, r_true_all, P_s);
    [err_th_B2, err_r_B2] = compute_matched_errors(th_B2_all, r_B2_all, ...
                                th_true_all, r_true_all, P_s);
    [err_th_B4, err_r_B4] = compute_matched_errors(th_B4_all, r_B4_all, ...
                                th_true_all, r_true_all, P_s);
    [err_th_B5, err_r_B5] = compute_matched_errors(th_B5_all, r_B5_all, ...
                                th_true_all, r_true_all, P_s);

    % RMSE (angle [deg], range [m])
    RMSE_theta_B1 = sqrt(mean(err_th_B1(:).^2)) * 180/pi;
    RMSE_theta_B2 = sqrt(mean(err_th_B2(:).^2)) * 180/pi;
    RMSE_theta_B4 = sqrt(mean(err_th_B4(:).^2)) * 180/pi;
    RMSE_theta_B5 = sqrt(mean(err_th_B5(:).^2)) * 180/pi;

    RMSE_r_B1 = sqrt(mean(err_r_B1(:).^2));
    RMSE_r_B2 = sqrt(mean(err_r_B2(:).^2));
    RMSE_r_B4 = sqrt(mean(err_r_B4(:).^2));
    RMSE_r_B5 = sqrt(mean(err_r_B5(:).^2));

    % NMSE (range) in dB
    r_true_mean_sq = mean(r_true_all(:).^2);
    NMSE_r_B1_dB = 10*log10(mean(err_r_B1(:).^2) / max(r_true_mean_sq, 1e-30));
    NMSE_r_B2_dB = 10*log10(mean(err_r_B2(:).^2) / max(r_true_mean_sq, 1e-30));
    NMSE_r_B4_dB = 10*log10(mean(err_r_B4(:).^2) / max(r_true_mean_sq, 1e-30));
    NMSE_r_B5_dB = 10*log10(mean(err_r_B5(:).^2) / max(r_true_mean_sq, 1e-30));

    % Failure rates (|err_r| / r_true > 0.50)
    fail_B1 = compute_fail_rate(err_r_B1, r_true_all);
    fail_B2 = compute_fail_rate(err_r_B2, r_true_all);
    fail_B4 = compute_fail_rate(err_r_B4, r_true_all);
    fail_B5 = compute_fail_rate(err_r_B5, r_true_all);

    % B4 (WB-CL-KL) convergence diagnostics
    b4_iters_mean  = mean(b4_n_iter(~isnan(b4_n_iter)));
    b4_conv_pct    = 100 * mean(b4_converged(~isnan(b4_converged)));

    % Runtime means
    rt_mean = @(x) mean(x(~isnan(x)));
    rt_B1_mean = rt_mean(rt_B1);
    rt_B2_mean = rt_mean(rt_B2);
    rt_B4_mean = rt_mean(rt_B4);
    rt_B5_mean = rt_mean(rt_B5);

    % ------------------------------------------------------------------
    %  5. Append one row per method to CSV
    % ------------------------------------------------------------------
    append_csv_row(csv_path, sweep_type, sweep_vec(s), SNR_dB_s, P_s, N_MC, ...
        'WB-BPD',    'X_full',    RMSE_theta_B1, RMSE_r_B1, NMSE_r_B1_dB, ...
        fail_B1, rt_B1_mean, NaN, NaN, NaN, NaN, NaN);

    append_csv_row(csv_path, sweep_type, sweep_vec(s), SNR_dB_s, P_s, N_MC, ...
        'WB-P-SOMP', 'R_hat_cell', RMSE_theta_B2, RMSE_r_B2, NMSE_r_B2_dB, ...
        fail_B2, rt_B2_mean, NaN, NaN, NaN, NaN, NaN);

    % B3 placeholder (not implemented)
    append_csv_row(csv_path, sweep_type, sweep_vec(s), SNR_dB_s, P_s, N_MC, ...
        'WB-BF-SOMP (N/A)', 'N/A', NaN, NaN, NaN, ...
        NaN, NaN, NaN, NaN, NaN, NaN, NaN);

    append_csv_row(csv_path, sweep_type, sweep_vec(s), SNR_dB_s, P_s, N_MC, ...
        'WB-CL-KL',  'R_hat_cell', RMSE_theta_B4, RMSE_r_B4, NMSE_r_B4_dB, ...
        fail_B4, rt_B4_mean, b4_iters_mean, b4_conv_pct, ...
        crb_theta_val, crb_r_val, NaN);

    append_csv_row(csv_path, sweep_type, sweep_vec(s), SNR_dB_s, P_s, N_MC, ...
        'WB-DL-OMP', 'X_full',    RMSE_theta_B5, RMSE_r_B5, NMSE_r_B5_dB, ...
        fail_B5, rt_B5_mean, NaN, NaN, NaN, NaN, NaN);

    % ------------------------------------------------------------------
    %  6. Save convergence data for Fig 10 (convergence sweep only)
    % ------------------------------------------------------------------
    if strcmp(sweep_type, 'convergence')
        save_convergence_data(out_dir, sweep_vec(s), L_hist_all, date_str);
    end

    fprintf('[MC]   Done. RMSE_r: B1=%.3f  B2=%.3f  B4=%.3f  B5=%.3f  [m]\n', ...
        RMSE_r_B1, RMSE_r_B2, RMSE_r_B4, RMSE_r_B5);
    fprintf('[MC]         CRB_r=%.3e [m]  crb_theta=%.4f [deg]\n', ...
        crb_r_val, crb_theta_val);

end  % for s

fprintf('[MC] COMPLETE. Results written to: %s\n', csv_path);

end  % run_monte_carlo_paperC


% =========================================================================
%  PER-REALISATION FUNCTION (parfor-safe)
% =========================================================================
function res = run_one_realisation(P, SNR_dB, mc_seed)
%RUN_ONE_REALISATION  Execute all five methods for one MC realisation.
%
%  This function is the ONLY location where rng is set inside the MC loop.
%  Called from the parfor body; must be fully self-contained.

rng(mc_seed, 'twister');

% ---- Generate one wideband channel realisation -------------------------
[X_full, Y_full, ~, theta_gen, r_gen, p_gen, N0_gen, W_comb] = ...
    wb_channel_gen_ofdm_nf(P, SNR_dB);

% Override to fixed scene for convergence sweep
if isfield(P, 'convergence_fixed_scene') && P.convergence_fixed_scene
    % Fixed geometry; re-generate with pinned scene by reusing the
    % channel generator with the theta/r passed via P struct.
    % (theta_gen and r_gen are random draws; override ground truth.)
    theta_gen(1) = P.conv_theta;
    r_gen(1)     = P.conv_r;
    % Re-generate channel data with fixed geometry.
    P_fixed = P;
    P_fixed.convergence_fixed_scene = false;   % prevent recursion
    rng(mc_seed, 'twister');   % reset seed for reproducibility
    [X_full, Y_full, ~, theta_gen, r_gen, p_gen, N0_gen, W_comb] = ...
        wb_channel_gen_ofdm_nf_fixed(P_fixed, SNR_dB, P.conv_theta, P.conv_r);
end

res.theta_true = theta_gen;
res.r_true     = r_gen;

% ---- Form compressed covariance cell (for B2 and B4) ------------------
R_hat_cell = cell(P.K_s, 1);
for k = 1 : P.K_s
    Yk            = Y_full(:, :, k);            % N_RF x N
    R_hat_cell{k} = (Yk * Yk') / P.N;
    R_hat_cell{k} = (R_hat_cell{k} + R_hat_cell{k}') / 2;  % enforce Hermitian
end

% ---- B1: WB-BPD (full-array) ------------------------------------------
t0 = tic;
try
    [th_B1, r_B1, ~, ~] = bpd_baseline(X_full, P);
    th_B1 = th_B1(:);
    r_B1  = r_B1(:);
catch ME
    warning('run_one_realisation:bpd_failed', ...
        'BPD failed (seed=%d): %s', mc_seed, ME.message);
    th_B1 = nan(P.d, 1);
    r_B1  = nan(P.d, 1);
end
res.rt_B1 = toc(t0);
res.th_B1 = th_B1;
res.r_B1  = r_B1;

% ---- B5: WB-DL-OMP (full-array) --------------------------------------
t0 = tic;
try
    [th_B5, r_B5, ~] = wb_dl_omp(X_full, P);
    th_B5 = th_B5(:);
    r_B5  = r_B5(:);
catch ME
    warning('run_one_realisation:dlomp_failed', ...
        'DL-OMP failed (seed=%d): %s', mc_seed, ME.message);
    th_B5 = nan(P.d, 1);
    r_B5  = nan(P.d, 1);
end
res.rt_B5 = toc(t0);
res.th_B5 = th_B5;
res.r_B5  = r_B5;

% ---- B2: WB-P-SOMP (compressed) --------------------------------------
t0 = tic;
try
    [th_B2, r_B2, ~] = wb_psomp(R_hat_cell, W_comb, P);
    th_B2 = th_B2(:);
    r_B2  = r_B2(:);
catch ME
    warning('run_one_realisation:psomp_failed', ...
        'P-SOMP failed (seed=%d): %s', mc_seed, ME.message);
    th_B2 = nan(P.d, 1);
    r_B2  = nan(P.d, 1);
end
res.rt_B2 = toc(t0);
res.th_B2 = th_B2;
res.r_B2  = r_B2;

% ---- B4: WB-CL-KL (compressed, BPD warm-start) -----------------------
% B1 must run first to provide the warm-start init.
if all(isfinite(th_B1)) && all(isfinite(r_B1))
    theta_init = th_B1(:);
    r_init     = r_B1(:);
else
    % Fallback: use B2 estimate or nominal values
    if all(isfinite(th_B2)) && all(isfinite(r_B2))
        theta_init = th_B2(:);
        r_init     = r_B2(:);
    else
        theta_init = repmat((P.theta_lo + P.theta_hi)/2, P.d, 1);
        r_init     = repmat(P.r_RD * (P.r_lo_fac + P.r_hi_fac)/2, P.d, 1);
    end
end
p_init = ones(P.d, 1) / P.d;   % equal-power warm-start

t0 = tic;
try
    [th_B4, r_B4, ~, ~, info_B4] = wb_clkl_driver( ...
        R_hat_cell, W_comb, theta_init, r_init, p_init, P);
    th_B4 = th_B4(:);
    r_B4  = r_B4(:);

    % Extract convergence info for Fig 10
    if isfield(info_B4, 'L_hist')
        res.L_hist = info_B4.L_hist(:);
    else
        res.L_hist = [];
    end
    if isfield(info_B4, 'n_iter')
        res.n_iter = info_B4.n_iter;
    else
        res.n_iter = NaN;
    end
    if isfield(info_B4, 'converged')
        res.converged = double(info_B4.converged);
    else
        res.converged = NaN;
    end

catch ME
    warning('run_one_realisation:clkl_failed', ...
        'WB-CL-KL failed (seed=%d): %s', mc_seed, ME.message);
    th_B4 = nan(P.d, 1);
    r_B4  = nan(P.d, 1);
    res.L_hist    = [];
    res.n_iter    = NaN;
    res.converged = NaN;
end
res.rt_B4 = toc(t0);
res.th_B4 = th_B4;
res.r_B4  = r_B4;

end  % run_one_realisation


% =========================================================================
%  CHANNEL GENERATOR WITH FIXED GEOMETRY (convergence sweep)
% =========================================================================
function [X_full, Y_full, H_true, theta_true, r_true, p_true, N0, W_comb] = ...
        wb_channel_gen_ofdm_nf_fixed(P, SNR_dB, theta_fix, r_fix)
%WB_CHANNEL_GEN_OFDM_NF_FIXED  Channel generator with pinned geometry.
%  Thin wrapper around wb_channel_gen_ofdm_nf that overrides the random
%  angle/range draw with fixed values (theta_fix, r_fix). Used only for
%  the convergence sweep (Fig 10) to produce reproducible KL-objective curves.
%
%  This is achieved by running wb_channel_gen_ofdm_nf in the normal way
%  but patching the snapshot tensor afterward using the fixed steering
%  vectors, so the noise realisation still varies across seeds while
%  geometry is frozen.

% Temporarily override range limits to ensure the fixed scene is in-range.
P_tmp         = P;
P_tmp.r_lo_fac = max(P.r_lo_fac, r_fix / P.r_RD * 0.99);
P_tmp.r_hi_fac = min(P.r_hi_fac, r_fix / P.r_RD * 1.01);

% Generate a realisation with the restricted range window.
% The exact theta/r draws will be close but not exactly (theta_fix, r_fix).
% Rebuild the snapshot tensor with the exact geometry below.
[~, ~, ~, ~, ~, p_true, N0, W_comb] = wb_channel_gen_ofdm_nf(P_tmp, SNR_dB);

% Fix geometry to the requested values.
theta_true = repmat(theta_fix, P.d, 1);
r_true     = repmat(r_fix,     P.d, 1);
p_true     = ones(P.d, 1) / P.d;

% Rebuild snapshot tensor with fixed geometry.
M   = P.M;
N   = P.N;
K_s = P.K_s;
m_bar = ((0:M-1) - (M-1)/2).';   % M x 1

X_full = zeros(M, N, K_s);
H_true = zeros(M, P.d, K_s);

for k = 1 : K_s
    alpha_k = P.alpha_k_vec(k);
    A_k     = zeros(M, P.d);
    for l = 1 : P.d
        omega_l  = (2*pi * P.d_ant / P.lambda_c) * cos(theta_true(l));
        kappa_l  = (pi  * P.d_ant^2 / P.lambda_c) * sin(theta_true(l))^2 / r_true(l);
        A_k(:,l) = exp(1j*alpha_k*(omega_l*m_bar - kappa_l*m_bar.^2)) * sqrt(M);
    end
    H_true(:,:,k) = A_k;
    S_k = sqrt(p_true(1)) * (randn(P.d, N) + 1j*randn(P.d, N)) / sqrt(2);
    W_k = sqrt(N0) * (randn(M, N) + 1j*randn(M, N)) / sqrt(2);
    X_full(:,:,k) = A_k * S_k + W_k;
end

% Compress
Y_full = zeros(P.N_RF, N, K_s);
for k = 1 : K_s
    Y_full(:,:,k) = W_comb' * X_full(:,:,k);
end

end  % wb_channel_gen_ofdm_nf_fixed


% =========================================================================
%  METRIC HELPERS
% =========================================================================
function [err_theta, err_r] = compute_matched_errors(th_hat_all, r_hat_all, ...
        th_true_all, r_true_all, P)
%COMPUTE_MATCHED_ERRORS  Hungarian-matched signed errors across N_MC trials.
%
%  th_hat_all  : d x N_MC  estimated angles [rad]
%  r_hat_all   : d x N_MC  estimated ranges [m]
%  th_true_all : d x N_MC  true angles [rad]
%  r_true_all  : d x N_MC  true ranges [m]
%
%  Returns:
%  err_theta : d x N_MC  signed angle errors [rad]
%  err_r     : d x N_MC  signed range errors [m]

d    = P.d;
N_MC = size(th_hat_all, 2);

err_theta = nan(d, N_MC);
err_r     = nan(d, N_MC);

% Range normalisation for cost matrix
r_ref = max(mean(r_true_all(:)), 1e-3);
eta   = 1 / r_ref^2;

for mc = 1 : N_MC

    th_hat  = th_hat_all(:, mc);
    r_hat   = r_hat_all(:, mc);
    th_true = th_true_all(:, mc);
    r_true  = r_true_all(:, mc);

    if any(isnan(th_hat)) || any(isnan(r_hat))
        % Propagate NaN for failed realisations
        err_theta(:, mc) = nan(d, 1);
        err_r(:, mc)     = nan(d, 1);
        continue
    end

    if d == 1
        % Trivial case: no permutation ambiguity
        err_theta(1, mc) = th_true(1) - th_hat(1);
        err_r(1, mc)     = r_true(1)  - r_hat(1);
    else
        % Hungarian matching: geodesic angular + normalised range cost
        cost_mat = zeros(d, d);
        for i = 1 : d
            for j = 1 : d
                dth = abs(th_true(i) - th_hat(j));
                dth = min(dth, pi - abs(th_true(i) + th_hat(j)));
                dr  = r_true(i) - r_hat(j);
                cost_mat(i, j) = dth^2 + eta * dr^2;
            end
        end
        assn = hungarian_assign(cost_mat);
        for i = 1 : d
            j = assn(i);
            err_theta(i, mc) = th_true(i) - th_hat(j);
            err_r(i, mc)     = r_true(i)  - r_hat(j);
        end
    end

end  % for mc

end  % compute_matched_errors


function assn = hungarian_assign(cost)
%HUNGARIAN_ASSIGN  Minimum-cost assignment using matchpairs (MATLAB 2016b+).
d_rows = size(cost, 1);
assn   = zeros(d_rows, 1);
try
    pairs = matchpairs(cost, 1e6);
    for k = 1 : size(pairs, 1)
        assn(pairs(k,1)) = pairs(k,2);
    end
    unmatched = find(assn == 0);
    for i = unmatched(:).'
        row         = cost(i,:);
        [~, j]      = min(row);
        assn(i)     = j;
    end
catch
    % Fallback: greedy nearest-neighbour
    used = false(size(cost, 2), 1);
    for i = 1 : d_rows
        row       = cost(i, :);
        row(used) = Inf;
        [~, j]    = min(row);
        assn(i)   = j;
        used(j)   = true;
    end
end
end  % hungarian_assign


function fail_rate = compute_fail_rate(err_r, r_true_all)
%COMPUTE_FAIL_RATE  Fraction of realisations where |err_r|/r_true > 0.50.
valid_mask = ~isnan(err_r) & ~isnan(r_true_all) & (r_true_all > 0);
if ~any(valid_mask(:))
    fail_rate = NaN;
    return
end
rel_err    = abs(err_r(valid_mask)) ./ r_true_all(valid_mask);
fail_rate  = 100 * mean(rel_err > 0.50);
end  % compute_fail_rate


% =========================================================================
%  CSV WRITER HELPERS
% =========================================================================
function write_csv_header(csv_path, sweep_type)
%WRITE_CSV_HEADER  Create CSV file with header row (first call only).
if isfile(csv_path)
    return
end

header = ['timestamp,sweep_type,sweep_value,SNR_dB,B_hz,K_s,' ...
          'M,N_RF,N,d,r_RD_m,N_MC,' ...
          'method,input_type,' ...
          'RMSE_theta_deg,RMSE_r_m,NMSE_r_dB,' ...
          'fail_rate_pct,runtime_s,' ...
          'clkl_avg_iters,clkl_conv_pct,' ...
          'crb_theta_deg,crb_r_m,' ...
          'notes'];

fid = fopen(csv_path, 'w');
if fid == -1
    error('run_monte_carlo_paperC: cannot create CSV: %s', csv_path);
end
fprintf(fid, '%s\n', header);
fclose(fid);
end  % write_csv_header


function append_csv_row(csv_path, sweep_type, sweep_val, SNR_dB, P, N_MC, ...
        method_name, input_type, RMSE_theta, RMSE_r, NMSE_r_dB, ...
        fail_rate, runtime, b4_iters, b4_conv_pct, crb_theta, crb_r, notes_val)
%APPEND_CSV_ROW  Write one result row to the master CSV.

ts = datestr(now, 'yyyy-mm-dd HH:MM:SS');

% Wideband columns
B_hz = P.B;
K_s  = P.K_s;

% CL-KL diagnostics columns (NaN for non-CL-KL methods)
c_iters    = fmt_num(b4_iters);
c_conv     = fmt_num(b4_conv_pct);
c_crb_th   = fmt_num(crb_theta);
c_crb_r    = fmt_num(crb_r);

if isnan(notes_val)
    c_notes = '';
else
    c_notes = fmt_num(notes_val);
end

fid = fopen(csv_path, 'a');
if fid == -1
    warning('append_csv_row: cannot open %s', csv_path);
    return
end
fprintf(fid, ['%s,%s,%s,%s,%s,%d,' ...
              '%d,%d,%d,%d,%s,%d,' ...
              '%s,%s,' ...
              '%s,%s,%s,' ...
              '%s,%s,' ...
              '%s,%s,' ...
              '%s,%s,' ...
              '%s\n'], ...
    ts, sweep_type, ...
    fmt_num(sweep_val), fmt_num(SNR_dB), fmt_num(B_hz), K_s, ...
    P.M, P.N_RF, P.N, P.d, fmt_num(P.r_RD), N_MC, ...
    method_name, input_type, ...
    fmt_num(RMSE_theta), fmt_num(RMSE_r), fmt_num(NMSE_r_dB), ...
    fmt_num(fail_rate), fmt_num(runtime), ...
    c_iters, c_conv, ...
    c_crb_th, c_crb_r, ...
    c_notes);
fclose(fid);
end  % append_csv_row


function s = fmt_num(x)
%FMT_NUM  Format scalar as string; empty string for NaN.
if isempty(x) || (isscalar(x) && isnan(x))
    s = '';
else
    s = sprintf('%.6g', x);
end
end  % fmt_num


% =========================================================================
%  CONVERGENCE DATA SAVER (Fig 10)
% =========================================================================
function save_convergence_data(out_dir, SNR_dB, L_hist_all, date_str)
%SAVE_CONVERGENCE_DATA  Write per-realisation KL-objective history to CSV.
%  One row per (SNR, realisation, iteration) triplet.

csv_name = sprintf('mc_convergence_Lhist_SNR%+.0fdB_%s.csv', SNR_dB, date_str);
csv_path = fullfile(out_dir, csv_name);

fid = fopen(csv_path, 'w');
if fid == -1
    warning('save_convergence_data: cannot create %s', csv_path);
    return
end
fprintf(fid, 'SNR_dB,mc_idx,iter_idx,L_val\n');

N_MC = numel(L_hist_all);
for mc = 1 : N_MC
    L_hist = L_hist_all{mc};
    if isempty(L_hist)
        continue
    end
    L_hist = L_hist(:);
    for it = 1 : numel(L_hist)
        fprintf(fid, '%.2f,%d,%d,%.10g\n', SNR_dB, mc, it, L_hist(it));
    end
end
fclose(fid);
end  % save_convergence_data
