function t40_e2b_driver(stage, out_dir, grid_mode)
%T40_E2B_DRIVER  T-40 / E2b discriminating experiment (Paper C strengthening).
%
%  t40_e2b_driver(stage)
%  t40_e2b_driver(stage, out_dir)
%  t40_e2b_driver(stage, out_dir, grid_mode)     grid_mode = 'legacy'|'refined'
%
%  REVISION T-54 (A6), 2026-10. Adds the finer-grid self-test specified in
%  PaperC_T53_T54_Specs.md Sec. 4, with its Addendum A and Addendum D.
%  DIAGNOSTIC ONLY: no T-40 or T-54 number is a result and no arm-B number
%  is usable (F-092). No T-40 constant, seed, gate geometry, threshold or
%  tolerance is changed; S6 keeps its definition and its 0.50 tolerance.
%    T54-1  loc_safe_clkl returns info.PhaseD_select and info.L_PhaseD; they
%           are written per trial (arm A) and per trial and subcarrier (arms
%           C_shared, C_percarrier), with a validity flag for every absorbed
%           error and every fallback anchor.
%    T54-2  wb_clkl_driver_pc takes the optional field P.phaseD_anchor_off
%           (default false); stage P2_SELFTEST_NOANCHOR sets it.
%    T54-3  R.degen_Cp, the degeneracy indicator of arm C_percarrier: a
%           control, REPORTED NEVER GATED, never called S6 (Addendum A.2).
%  T-54 stages, in run order. Each needs an explicit grid_mode and exactly
%  one file t54_runtag_<tag>.txt in out_dir (the run sheet creates it):
%    'P0'                    S-0: smoke pass through the S-1..S-4 code path,
%                                 benchmark at both scan grids, stop rule.
%    'T54_GT3'                    gate G-T3 and the T54-2 behaviour check.
%    'P2_SELFTEST'           S-1 (legacy), S-2 (refined): S1-S6 over
%                                 CFG.B_grid, S6 per leg and per bandwidth.
%    'P2_SELFTEST_NOANCHOR'  S-3 (refined), S-4 (legacy): the same with the
%                                 Phase D anchor disabled (T54-2).
%  In the stage list below, 'P0_bench' is the unchanged T-40 benchmark and
%  'P2_selftest' is now the T-54 self-test.
%
%  REVISION v3, 2026-09-03. Supersedes v2. v2's P0/P4 stages stand; v2's
%  gate design is amended after the T-40 P4 consultant review (Sol + Flash,
%  independent). See PaperC_T40_Gate_Amendment_2026-09-03.md, which must be
%  filed BEFORE this file is run. Fixes 6-10 below are the amendment.
%
%    FIX 6  GLOBAL-ERROR INDICATOR. The published rule |err_r|/r_true > 0.50
%           is retained for continuity but is NOT a wrong-mode indicator: a
%           far-grid-node lock at r_true=3.966 -> r_hat=4.819 is 21.5 percent
%           relative error and scores as non-gross. The primary indicator is
%           now (r_hat outside the scene box [r_lo, r_hi]) OR the legacy rule.
%           Both are recorded separately. Measured on the P4 data, the new
%           rule cuts the block-to-block spread of the conditional statistic
%           from 13.95/13.97/15.92 dB to 1.09/2.80/3.84 dB.
%    FIX 7  GATE STATISTIC is conditional RMSE (CRMSE) on each arm's OWN
%           inliers, not median absolute error. Both consultants rejected the
%           median independently: it is nearly blind to rare events, its
%           sampling distribution differs from RMSE so the pre-registered
%           3.0/1.0 dB thresholds would not transfer, and conditioning on
%           jointly-non-gross trials couples the arms. CRMSE is an L2
%           statistic in the same units as RMSE, so the thresholds DO
%           transfer, and per-arm conditioning preserves benchmark
%           independence. The paired common-success estimand is reported as a
%           clearly-labelled secondary, never as the gate.
%    FIX 8  RELIABILITY GATE is a PAIRED McNemar-type comparison on the
%           discordant pairs, with a pre-registered absolute risk-difference
%           margin. Overlapping marginal Wilson intervals waste the pairing
%           and test nothing.
%    FIX 9  FRESH SEED BASE for every gate stage (2000000). P4 was run at
%           700000. The amended estimands are therefore confirmed on data not
%           used to choose them.
%    FIX 10 SECOND PINNED GEOMETRY at r = 3.90 m, in the far quartile where
%           the failure mode lives. REPORTED, NOT GATED: it was added after
%           seeing P4 and is labelled post hoc throughout.
%
%  Superseded v2 notes retained below.
%
%  REVISION v2, 2026-09-02. Supersedes the 2026-09-01 harness, whose MAIN
%  results are VOID. Five defects were found in that run and are fixed here.
%
%    FIX 1  Whitening used a non-conjugate transpose (Tw.') in the arm-B path
%           while the covariance used Tw'. Arm B was fitting a mismatched
%           manifold and every arm-B number from run 20260901 is void. The
%           whitened pair is now built ONCE inside loc_whiten, and the self
%           test exercises that same object rather than a re-derivation.
%    FIX 2  loc_pack never populated n_iter or converged (all NaN in the
%           first run) and there was no boundary flag. Both are now recorded
%           per arm, plus an explicit clamp indicator.
%    FIX 3  New mandatory stage P4_REPRO: arm A alone against the published
%           production CSV. The first ladder had no reproduction gate, so it
%           ran clean while arm A sat 13.4 dB off its own published baseline.
%    FIX 4  The gate statistic moves from RMSE to a two-part gate: median
%           absolute range error on jointly non-gross trials (precision),
%           plus a symmetric gross-rate comparison with Wilson intervals
%           (robustness). RMSE is still reported but is no longer the gate:
%           in run 20260901 arms A and C_shared had IDENTICAL median error at
%           every point and the entire RMSE difference was a 1 to 8 percent
%           clamp-failure tail. Thresholds are carried over unchanged
%           (3.0 / 1.0 dB), so this is a change of statistic, not of bar.
%    FIX 5  Degeneracy diagnostics. Run 20260901 had disp_omega identically
%           zero: the K_s per-subcarrier estimates were bit-identical, so the
%           rival arm received no data diversity at all. The fraction of
%           trials with degenerate per-subcarrier estimates is now measured
%           and reported for arms B and C, and gated in P2 (test S6).
%
%  Also changed: CFG.use_refined now defaults to FALSE. The refined scan grid
%  is under suspicion of having manufactured the far-clamp failure mode; it is
%  now an experimental factor, exercised by running P4_REPRO at both settings.
%
%  STAGES (run in this order; each is a gate)
%  ------------------------------------------
%    'P0_bench'    micro-benchmark; runtime projection. No verdict.
%    'P4_repro'    REPRODUCTION GATE. Arm A alone vs the production CSV.
%                  Run at BOTH grid modes before any other stage.
%    'P1_floor'    grid-quantisation floor probe at SNR = 30 dB.
%    'P2_selftest' harness self-consistency checks (no science).
%    'P3_flat'     frequency-flat control, alpha_k == 1. Yields Delta_arch.
%    'MAIN_SNR'    SNR sweep, pinned geometry (primary) + randomised.
%    'MAIN_BW'     bandwidth sweep at SNR* (mechanism gate).
%    'ESCALATE'    nominal point only, N_MC = 2400 (INTERMEDIATE band only).
%
%  ARMS
%  ----
%    A            wb_clkl_driver_pc over all K_s subcarriers jointly.
%    C_shared     wb_clkl_driver_pc run K_s times, one subcarrier each,
%                 sharing the SINGLE wideband BPD warm start given to arm A.
%    C_percarrier as C_shared, but BPD is re-run on X_full(:,:,k) alone.
%    B_matched    nf_clkl_pc at subcarrier k, lambda_eff = lambda_c/alpha_k.
%    B_naive      nf_clkl_pc with lambda_c for every k, alpha_k removed in
%                 the combination step.
%
%  COMBINATION (pre-registered, unchanged)
%    PRIMARY  : arithmetic mean in (omega, kappa) referred to the carrier.
%    SECONDARY: median; 20 percent trimmed mean; KL-argmin selection.
%
%  REQUIRES on path
%    setup_production_P_v4.m, wb_channel_gen_ofdm_nf.m, bpd_baseline.m,
%    wb_nf_fresnel_steer.m, nf_usw_steer.m,
%    wb_clkl_driver_pc.m   (patched copy -- spec Sec. 1.3, patch W1)
%    nf_clkl_pc.m          (patched copy -- spec Sec. 1.3, patch B1/B2)
%    wb_clkl_estimator.m   (called by wb_clkl_driver_pc)
%  T-54 stages also need on path
%    wb_clkl_driver_pc_prepatch.m  (gate G-T3: the pre-T-54 copy, renamed)
%    run_monte_carlo_paperC.m      (gate G-T1: read as text, never called)
%
%  Author-side execution. MATLAB R2025b. Provenance of every CSV: P-RUN.
%  Encoding: 7-bit ASCII only.

narginchk(1, 3);
if nargin < 2 || isempty(out_dir);   out_dir   = 'results_t40'; end
if nargin < 3 || isempty(grid_mode); grid_mode = '';            end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

stage = upper(char(stage));
rng(20260901, 'twister');

CFG = t40_config();
if ~isempty(grid_mode)
    switch lower(grid_mode)
        case 'legacy';  CFG.use_refined = false;
        case 'refined'; CFG.use_refined = true;
        otherwise; error('grid_mode must be ''legacy'' or ''refined''.');
    end
end
gtag = loc_tern(CFG.use_refined, 'refined', 'legacy');

stamp = datestr(now, 'yyyymmdd_HHMMSS');

% ---- T-54 (A6) stages: run tag and file prefix ----------------------------
% The T-54 stages of one run share one run tag, read from the single file
% t54_runtag_<tag>.txt that the run sheet creates in out_dir before S-0.
% Every file a T-54 stage writes carries the tag in its name, so that one
% capture selects the whole run and nothing else, also across midnight.
t54_names = {'P0', 'T54_GT3', 'P2_SELFTEST', 'P2_SELFTEST_NOANCHOR'};
is_t54    = any(strcmp(stage, t54_names));
run_tag   = '';
fpre      = '';
if is_t54
    if isempty(grid_mode)
        error('t40_e2b_driver:T54Grid', ...
              'T-54 stage %s needs an explicit grid_mode (''legacy'' or ''refined'').', stage);
    end
    run_tag = loc_t54_runtag(out_dir);
    fpre    = fullfile(out_dir, sprintf('t54_%s_%s_%s', run_tag, lower(stage), gtag));
    diary([fpre '_console_' stamp '.txt']);
else
    diary(fullfile(out_dir, sprintf('t40_%s_%s_%s_console.txt', ...
          lower(stage), gtag, stamp)));
end
F54 = struct('pre', fpre, 'stamp', stamp, 'run_tag', run_tag, ...
             'stage', stage, 'gtag', gtag);
diary on;
fprintf('=============================================================\n');
fprintf('  T-40 / E2b  v2   stage = %s   grid = %s   %s\n', stage, gtag, stamp);
fprintf('=============================================================\n');

switch stage
    case 'P0_BENCH';    t40_bench(CFG, out_dir, stamp, gtag);
    case 'P4_REPRO';    t40_repro(CFG, out_dir, stamp, gtag);
    case 'P1_FLOOR';    t40_floor(CFG, out_dir, stamp, gtag);
    case {'P0', 'T54_GT3', 'P2_SELFTEST', 'P2_SELFTEST_NOANCHOR'}
        t54_dispatch(CFG, F54, grid_mode);
    case 'P3_FLAT';     t40_sweep(CFG, out_dir, stamp, gtag, 'snr', ...
                            CFG.SNR_star, CFG.N_MC, true, 'pinned', 'P3_flat');
    case 'MAIN_SNR'
        t40_sweep(CFG, out_dir, stamp, gtag, 'snr', CFG.SNR_grid, CFG.N_MC, ...
                  false, 'pinned',     'MAIN_SNR_pinned');
        t40_sweep(CFG, out_dir, stamp, gtag, 'snr', CFG.SNR_star, CFG.N_MC, ...
                  false, 'randomised', 'MAIN_SNR_random');
        t40_sweep(CFG, out_dir, stamp, gtag, 'snr', CFG.SNR_grid, CFG.N_MC, ...
                  false, 'pinned_far', 'MAIN_SNR_pinnedfar');
    case 'MAIN_BW'
        t40_sweep(CFG, out_dir, stamp, gtag, 'bandwidth', CFG.B_grid, ...
                  CFG.N_MC, false, 'pinned', 'MAIN_BW');
    case 'ESCALATE'
        t40_sweep(CFG, out_dir, stamp, gtag, 'snr', CFG.SNR_star, ...
                  CFG.N_MC_esc, false, 'pinned', 'ESCALATE');
    otherwise
        error('t40_e2b_driver: unknown stage %s', stage);
end

fprintf('\n[T-40] stage %s (%s grid) complete.\n', stage, gtag);
diary off;
end


% =========================================================================
%  CONFIGURATION -- every pre-registered constant lives here and nowhere else
% =========================================================================
function CFG = t40_config()

CFG.r_lo_fac   = 0.05;        % LOCKED. README's [0.2126,...] implies 0.01;
                              % setup_production_P_v4 and F-062 both say 0.05.
CFG.r_hi_fac   = 0.20;        % LOCKED box [1.0631, 4.2525] m. See spec A-01.
CFG.SNR_star   = 5;           % dB, nominal operating point (gate is here)
CFG.SNR_grid   = [-5 0 5 10 15];
CFG.B_grid     = [100 200 400 600] * 1e6;   % 800 MHz excluded ex ante (L54)
CFG.N_MC       = 600;
CFG.N_MC_esc   = 2400;
CFG.N_MC_bench = 8;
CFG.N_MC_floor = 100;
CFG.SNR_floor  = 30;
CFG.pin_theta  = 40 * pi/180;
CFG.pin_r      = 2.13;        % gate geometry (pre-registered, unchanged)
CFG.pin_r_far  = 3.90;        % FIX 10: reported, never gated (post hoc)
CFG.seed_base  = 2000000;     % FIX 9: fresh. P4 ran at 700000.
CFG.n_boot     = 2000;
CFG.boot_seed  = 20260901;

% ---- FIX 3: reproduction gate anchors -----------------------------------
% Source: results/regime_probe_rhi020/mc_snr_sweep_20260528_204613.csv,
% method WB-CL-KL, N_MC = 600, randomised box r_hi_fac = 0.20, B = 400 MHz.
CFG.repro_SNR      = [0 5 10];
CFG.repro_RMSE_ref = [0.0444729 0.0279913 0.0198205];   % m
CFG.repro_tol_dB   = 1.5;
CFG.N_MC_repro     = 600;
CFG.repro_seed_blocks = [0 100000 200000 300000];   % offsets on seed_base

% ---- Scan grid ----------------------------------------------------------
% Legacy IS the production grid and is now the default. The refined grid is
% an experimental factor, not a fixed choice: run P4_REPRO at both.
CFG.Q_scan_th_ref = 768;      % refined  (production/legacy: 192)
CFG.Q_scan_u_ref  = 2048;     % refined  (production/legacy: 256)
CFG.use_refined   = false;

CFG.do_C_percarrier = true;
CFG.floor_gate_frac = 0.30;
CFG.gross_frac      = 0.50;   % legacy rule, retained for continuity only
CFG.use_box_rule    = true;   % FIX 6: primary global-error indicator
CFG.rd_margin       = 0.01;   % FIX 8: pre-registered absolute risk-difference
                              % margin on the global-error probability

% ---- T-54 (A6): finer-grid self-test, DIAGNOSTIC ONLY ---------------------
% Spec = PaperC_T53_T54_Specs.md Sec. 4 with Addendum A. A line marked
% [D-gN] holds the option of gap gN that Addendum D registers (author
% decision before any A6 data exist). No constant above this block is
% changed by T-54.
CFG.t54_nd           = 12;        % S6 trials per cell (Spec 4.4; unchanged)
CFG.t54_s6_offset    = 500;       % S6 seeds = seed_base + 500 + (1:nd); unchanged
CFG.t54_s6_tol       = 0.50;      % S6 tolerance (F-091); NOT amended
CFG.t54_n_large      = 200;       % large-sample estimate (Spec 4.4); never S6
CFG.t54_large_offset = 400000;    % [D-g4] its seeds = seed_base + 400000 + (1:n)
CFG.t54_gt3_n        = 50;        % G-T3 trials (Spec 4.7)
CFG.t54_gt3_offset   = 450000;    % [D-g5] G-T3 seeds = seed_base + 450000 + (1:50)
CFG.t54_legs         = {'pinned', 'pinned_far'};            % [D-g9] legs run
CFG.t54_reg_leg      = 'pinned';  % [D-g6] registered cell: leg
CFG.t54_reg_B        = 400e6;     % [D-g6] registered cell: bandwidth [Hz]
CFG.t54_anchor_arms  = {'A', 'C_shared', 'C_percarrier'};   % [D-g8] arms of T54-2
CFG.t54_gt2_mode     = 'threshold';   % [D-g3] 'threshold' | 'exact' | 'report'
CFG.t54_stop_ratio   = 30;        % stop rule (Spec 4.7); NOT amended
CFG.t54_invalid_halt = true;      % [D-g18] invalid S6 trial at the registered cell halts
CFG.t54_dflt_halt    = true;      % [D-g13] a rejected setup call or a defaulted field halts
CFG.t54_du_ref       = 8.615869224890557e-04;   % [D-g24] G-T4 reference, recomputed from the locked box
CFG.t54_du_tol       = 1e-12;     % G-T4 relative tolerance (Spec 4.7); NOT amended
CFG.t54_du_box       = (1/(0.05*21.2625/2) - 1/(0.20*21.2625*2)) / 2047;
                                  % reported beside G-T4: the step the locked box
                                  % (r_lo_fac, r_hi_fac, r_RD, u_margin = 2) implies
CFG.t54_smoke_nd     = 2;         % S-0 smoke pass: trials per set (not a statistic)
CFG.t54_smoke_B      = [100 400] * 1e6;   % S-0 smoke pass: bandwidths [Hz]
% The four flags below are set by the T-54 stage functions only.
CFG.t54_anchor_off   = false;     % true in stage P2_SELFTEST_NOANCHOR (T54-2)
CFG.t54_prepatch     = false;     % true in T54_GT3 only: call the pre-patch copy
CFG.t54_force_false  = false;     % true in T54_GT3 only: P.phaseD_anchor_off = false
CFG.t54_smoke        = false;     % true in the S-0 smoke pass only
end


% =========================================================================
%  P STRUCT
% =========================================================================
function [P, pinfo] = t40_P(CFG, B_hz)
%T40_P  Self-sufficient parameter struct. Any field the live
%       setup_production_P[_v4] does not supply is filled from the locked
%       set; every derived quantity is recomputed; the box is pinned here.
%       pinfo (T-54) records whether the live setup call was rejected and
%       which fields were filled from the locked set, so that a T-54 stage
%       can refuse a silently defaulted model (L-34).

pinfo = struct('setup_rejected', false, 'n_filled', 0, 'filled', '');
P = struct();
try
    if exist('setup_production_P_v4', 'file') == 2
        P = setup_production_P_v4('snr', 'full');
    elseif exist('setup_production_P', 'file') == 2
        P = setup_production_P('snr', 'full');
    end
catch ME
    fprintf('  [P] live setup rejected the call (%s)\n', ME.message);
    fprintf('  [P] proceeding on locked defaults.\n');
    P = struct();
    pinfo.setup_rejected = true;
end
if ~isstruct(P); P = struct(); end

[P, miss] = loc_fill_P(P, CFG);
pinfo.n_filled = numel(miss);
if ~isempty(miss); pinfo.filled = strjoin(miss, ' '); end

if nargin >= 2 && ~isempty(B_hz)
    Delta_f_fixed = 25e6;                       % Lesson L29, locked
    P.B       = B_hz;
    P.K_s     = min(round(B_hz / Delta_f_fixed), 512);
    P.Delta_f = B_hz / P.K_s;
    k_idx     = (-(P.K_s/2) : (P.K_s/2 - 1)).';
    P.alpha_k_vec = 1 + k_idx * P.Delta_f / P.fc;
    P.k_indices   = k_idx;
    P.K           = P.K_s;
end

if CFG.use_refined
    P.Q_scan_th = CFG.Q_scan_th_ref;
    P.Q_scan_u  = CFG.Q_scan_u_ref;
else
    P.Q_scan_th = 192;
    P.Q_scan_u  = 256;
end

fprintf('  [P] M=%d N_RF=%d N=%d d=%d K_s=%d B=%.0f MHz\n', ...
        P.M, P.N_RF, P.N, P.d, P.K_s, P.B/1e6);
fprintf('  [P] r_RD=%.4f m  scene box=[%.4f, %.4f] m  (r_lo_fac=%.2f r_hi_fac=%.2f)\n', ...
        P.r_RD, P.r_lo_fac*P.r_RD, P.r_hi_fac*P.r_RD, P.r_lo_fac, P.r_hi_fac);
fprintf('  [P] u=[%.6f, %.6f] 1/m -> admissible r=[%.4f, %.4f] m\n', ...
        P.u_min, P.u_max, 1/P.u_max, 1/P.u_min);
fprintf('  [P] Q_scan=[%d, %d]\n', P.Q_scan_th, P.Q_scan_u);
end


function [P, miss] = loc_fill_P(P, CFG)
%LOC_FILL_P  Locked parameter set. Fills only ABSENT fields; recomputes every
%            derived quantity unconditionally; pins the box.
D.c = 3e8;          D.fc = 28e9;
D.M = 64;           D.N_RF = 8;      D.N = 64;    D.d = 1;
D.theta_lo = 20*pi/180;              D.theta_hi = 60*pi/180;
D.u_margin = 2.0;
D.B = 400e6;        D.K_s = 16;
D.Q_theta = 256;    D.beta_delta = 1.2;
D.G_theta = 128;    D.G_r = 128;
D.lambda_reg = 1e-4; D.max_iter = 200; D.tol_clkl = 1e-5;
D.alpha_p = 0.5;    D.ls_beta = 0.5;  D.ls_sigma = 1e-4;
D.eps_reg = 1e-3;
D.N_seed = 50;      D.rng_seed_W = 0;
D.wb_gen_write_csv = false;
D.bpd_write_csv = false;
D.use_riviello_snr_axis = false;

f = fieldnames(D); miss = {};
for i = 1:numel(f)
    if ~isfield(P, f{i}); P.(f{i}) = D.(f{i}); miss{end+1} = f{i}; end %#ok<AGROW>
end

P.c0       = P.c;
P.lambda   = P.c / P.fc;
P.lambda_c = P.lambda;
P.d_ant    = P.lambda_c / 2;
P.r_RD     = 2 * ((P.M - 1) * P.d_ant)^2 / P.lambda_c;
P.r_lo_fac = CFG.r_lo_fac;
P.r_hi_fac = CFG.r_hi_fac;
P.u_min    = 1 / (P.r_hi_fac * P.r_RD * P.u_margin);
P.u_max    = 1 / (P.r_lo_fac * P.r_RD / P.u_margin);
P.Delta_f  = P.B / P.K_s;
k_idx      = (-(P.K_s/2) : (P.K_s/2 - 1)).';
P.k_indices   = k_idx;
P.alpha_k_vec = 1 + k_idx * P.Delta_f / P.fc;
P.K        = P.K_s;
P.d_max    = P.d;

if ~isempty(miss)
    fprintf('  [P] filled from locked defaults (%d): %s\n', ...
            numel(miss), strjoin(miss, ', '));
end

assert(P.M == 64 && P.N_RF == 8 && P.N == 64 && P.d == 1, ...
    't40_P: live setup returned a non-production array configuration.');
assert(abs(P.fc - 28e9) < 1, 't40_P: carrier is not 28 GHz.');
assert(abs(P.r_RD - 21.2625) < 1e-3, 't40_P: r_RD is not 21.2625 m.');
end


% =========================================================================
%  P0 -- MICRO-BENCHMARK
% =========================================================================
function t40_bench(CFG, out_dir, stamp, gtag)
P = t40_P(CFG, 400e6);
fprintf('\n### P0 micro-benchmark (N_MC = %d) ###\n', CFG.N_MC_bench);
t = zeros(CFG.N_MC_bench, 5);
for i = 1:CFG.N_MC_bench
    R = t40_trial(P, CFG, CFG.SNR_star, CFG.seed_base + i, 'pinned', false);
    t(i,:) = [R.rt_gen R.rt_bpd R.rt_A R.rt_C R.rt_B];
end
m = mean(t, 1);
lab = {'gen','bpd_wb','armA','armsC','armsB'};
for i = 1:5; fprintf('  %-8s %.3f s/trial\n', lab{i}, m(i)); end
tot = sum(m);
np  = max(1, feature('numcores') - 2);
fprintf('  total %.3f s/trial | %d workers assumed\n', tot, np);
fprintf('  P4_REPRO  (%d pts x %d MC, arm A only) : %.1f min\n', ...
        numel(CFG.repro_SNR), CFG.N_MC_repro, ...
        (m(1)+m(2)+m(3))*numel(CFG.repro_SNR)*CFG.N_MC_repro/60/np);
fprintf('  MAIN_SNR  (%d pts x %d MC) : %.1f min\n', ...
        numel(CFG.SNR_grid), CFG.N_MC, tot*numel(CFG.SNR_grid)*CFG.N_MC/60/np);
fprintf('  MAIN_BW   (%d pts x %d MC) : %.1f min\n', ...
        numel(CFG.B_grid), CFG.N_MC, tot*numel(CFG.B_grid)*CFG.N_MC/60/np);
fprintf('  GATE P0: if MAIN_SNR projection > 600 min, set\n');
fprintf('           CFG.do_C_percarrier = false and re-run P0.\n');
writetable(array2table([m tot], 'VariableNames', [lab {'total'}]), ...
    fullfile(out_dir, sprintf('t40_P0_bench_%s_%s.csv', gtag, stamp)));
end


% =========================================================================
%  P4 -- REPRODUCTION GATE (FIX 3). Arm A alone, randomised geometry.
% =========================================================================
function t40_repro(CFG, out_dir, stamp, gtag)
fprintf('\n### P4 reproduction gate (arm A only, randomised box) ###\n');
fprintf('  Reference: mc_snr_sweep_20260528_204613.csv, WB-CL-KL, N_MC=600.\n');
fprintf('  Seed blocks: %s (offsets on seed_base = %d)\n', ...
        mat2str(CFG.repro_seed_blocks), CFG.seed_base);

TT = table(); TR = table(); first = true;

for s = 1:numel(CFG.repro_SNR)
    P      = t40_P(CFG, 400e6);
    SNR_dB = CFG.repro_SNR(s);
    ref    = CFG.repro_RMSE_ref(s);
    N_MC   = CFG.N_MC_repro;

    for bkt = 1:numel(CFG.repro_seed_blocks)
        base = CFG.seed_base + CFG.repro_seed_blocks(bkt);
        er = nan(N_MC,1); rt = nan(N_MC,1); tt = nan(N_MC,1);
        rh = nan(N_MC,1); ni = nan(N_MC,1); cv = nan(N_MC,1);
        bd = nan(N_MC,1); gr = nan(N_MC,1);

        parfor mc = 1:N_MC
            R = t40_trial(P, CFG, SNR_dB, base + mc, 'randomised', true); %#ok<PFBNS>
            er(mc) = R.err_r.A;    rt(mc) = R.r_true;  tt(mc) = R.theta_true;
            rh(mc) = R.r_hat.A;    ni(mc) = R.n_iter.A;
            cv(mc) = R.converged.A; bd(mc) = R.boundary.A; gr(mc) = R.gross.A;
        end

        rmse  = sqrt(mean(er.^2, 'omitnan'));
        clean = gr == 0;
        medae = median(abs(er(clean)), 'omitnan');
        rmse_c = sqrt(mean(er(clean).^2, 'omitnan'));
        [mx, imx] = max(abs(er));
        dev   = 20*log10(rmse/ref);

        fprintf(['  SNR=%+5.1f blk%d  RMSE=%.5f (dev %+6.2f dB) | ' ...
                 'RMSE_nongross=%.5f | med|e|=%.5f\n'], ...
                SNR_dB, bkt, rmse, dev, rmse_c, medae);
        fprintf(['                 max|e|=%.4f m at trial %d ' ...
                 '(r_true=%.3f, r_hat=%.3f, th_true=%.2f deg) | ' ...
                 'gross=%.4f conv=%.3f\n'], ...
                mx, imx, rt(imx), rh(imx), tt(imx)*180/pi, ...
                mean(gr,'omitnan'), mean(cv,'omitnan'));

        TT = [TT; table(string(gtag), SNR_dB, bkt, base, N_MC, rmse, rmse_c, ...
              ref, dev, medae, mx, imx, rt(imx), rh(imx), tt(imx)*180/pi, ...
              mean(gr,'omitnan'), mean(bd,'omitnan'), mean(cv,'omitnan'), ...
              'VariableNames', {'grid','SNR_dB','block','seed_base','N_MC', ...
              'RMSE_r_m','RMSE_nongross_m','RMSE_ref_m','dev_dB', ...
              'med_abs_err_m','max_abs_err_m','argmax_trial','r_true_at_max', ...
              'r_hat_at_max','theta_true_at_max_deg','gross_rate', ...
              'clamp_rate','conv_rate'})]; %#ok<AGROW>

        TRb = table(repmat(string(gtag),N_MC,1), repmat(SNR_dB,N_MC,1), ...
              repmat(bkt,N_MC,1), (base+(1:N_MC)).', tt*180/pi, rt, rh, er, ...
              ni, cv, gr, bd, ...
              'VariableNames', {'grid','SNR_dB','block','seed', ...
              'theta_true_deg','r_true_m','r_hat_m','err_r_m','n_iter', ...
              'converged','gross','boundary'});
        writetable(TRb, fullfile(out_dir, ...
            sprintf('t40_P4_repro_%s_trials_%s.csv', gtag, stamp)), ...
            'WriteMode', loc_tern(first,'overwrite','append'));
        first = false;
    end

    ix = TT.SNR_dB == SNR_dB;
    fprintf('  --> SNR=%+5.1f: RMSE spread across blocks = %.2f dB\n\n', ...
            SNR_dB, 20*log10(max(TT.RMSE_r_m(ix))/min(TT.RMSE_r_m(ix))));
end

writetable(TT, fullfile(out_dir, ...
    sprintf('t40_P4_repro_%s_%s.csv', gtag, stamp)));
fprintf('  Verdict is adjudicated outside the run. See T-45 addendum.\n');
end

% =========================================================================
%  P1 -- GRID-QUANTISATION FLOOR PROBE
% =========================================================================
function t40_floor(CFG, out_dir, stamp, gtag)
P = t40_P(CFG, 400e6);
du    = (P.u_max - P.u_min) / (P.Q_scan_u - 1);
u_pin = 1 / CFG.pin_r;
dr    = du / u_pin^2;
pred  = dr / sqrt(12);
fprintf('\n### P1 floor probe ###\n');
fprintf('  Q_scan_u = %d | du = %.6e 1/m | dr(r*) = %.6e m\n', P.Q_scan_u, du, dr);
fprintf('  predicted uniform-quantisation RMSE floor = %.6e m  [P-CALC]\n', pred);

n = CFG.N_MC_floor;
eA = nan(n,1); eB = nan(n,1); gA = nan(n,1); gB = nan(n,1);
parfor i = 1:n
    R = t40_trial(P, CFG, CFG.SNR_floor, CFG.seed_base + 900000 + i, ...
                  'pinned', false); %#ok<PFBNS>
    eA(i) = R.err_r.A;  eB(i) = R.err_r.B_matched_mean;
    gA(i) = R.gross.A;  gB(i) = R.gross.B_matched_mean;
end
% The floor must be measured on NON-GROSS trials, else it is a failure-rate
% statistic wearing a quantisation label. Run 20260901 failed exactly here.
flA_all   = sqrt(mean(eA.^2, 'omitnan'));
flA_clean = sqrt(mean(eA(gA == 0).^2, 'omitnan'));
flB_clean = sqrt(mean(eB(gB == 0).^2, 'omitnan'));
fprintf('  arm A at SNR=%d dB: RMSE(all)=%.6e | RMSE(non-gross)=%.6e | gross=%.4f\n', ...
        CFG.SNR_floor, flA_all, flA_clean, mean(gA,'omitnan'));
fprintf('  arm B at SNR=%d dB: RMSE(non-gross)=%.6e | gross=%.4f\n', ...
        CFG.SNR_floor, flB_clean, mean(gB,'omitnan'));
fprintf('  GATE P1: the NON-GROSS floor must be < %.2f x CRMSE at SNR*.\n', ...
        CFG.floor_gate_frac);
fprintf('  Compare against MAIN CRMSE_own_m, NOT against RMSE_r_m.\n');
fprintf('  NOTE: gross here is the FIX 6 rule (outside scene box OR legacy).\n');
writetable(table(string(gtag), P.Q_scan_u, du, dr, pred, flA_all, flA_clean, ...
    flB_clean, mean(gA,'omitnan'), mean(gB,'omitnan'), ...
    'VariableNames', {'grid','Q_scan_u','du_per_m','dr_step_m', ...
    'pred_floor_m','rmse_A_all_m','floor_A_nongross_m', ...
    'floor_B_nongross_m','gross_A','gross_B'}), ...
    fullfile(out_dir, sprintf('t40_P1_floor_%s_%s.csv', gtag, stamp)));
end


% =========================================================================
%  P2 -- HARNESS SELF-TESTS
% =========================================================================
function [GT, halt] = t40_selftest(CFG, F, GT, anchor_off, smoke)
%T40_SELFTEST  P2 self-tests as extended by T-54 (Spec Sec. 4.4; Add. A.2).
%  S1-S5 run at every bandwidth of the stage on the pinned leg, with the
%  Phase D anchor on: they check the harness, not the T-54 arm. S6, with its
%  definition and its 0.50 tolerance unchanged, and the T-54 reporting
%  statistics are produced per leg and per bandwidth by t54_cell.
%  anchor_off = true is stage P2_SELFTEST_NOANCHOR (change T54-2).
%  smoke = true is the S-0 smoke pass: CFG.t54_smoke_B and CFG.t54_smoke_nd
%  trials per set; nothing it prints is a statistic.
CFG.t54_anchor_off = logical(anchor_off);
CFG.t54_smoke      = logical(smoke);
slab = loc_t54_slab(CFG.use_refined, CFG.t54_anchor_off, CFG.t54_smoke);
if CFG.t54_smoke
    B_list = CFG.t54_smoke_B;  nd = CFG.t54_smoke_nd;  nl = CFG.t54_smoke_nd;
else
    B_list = CFG.B_grid;       nd = CFG.t54_nd;        nl = CFG.t54_n_large;
end
C5 = CFG;  C5.t54_anchor_off = false;     % S1-S5: harness checks, anchor on
halt = {};
Gall = cell(0, 9);
REG  = struct('s6', [], 'large', []);
ok_harness = true;
n_cells = 0;  n_s6 = 0;  n_sel = 0;  n_inv = 0;  n_tr = 0;
fprintf('\n### P2 self-tests -- T-54 stage %s (grid %s, Phase D anchor %s) ###\n', ...
        slab, loc_tern(CFG.use_refined, 'refined', 'legacy'), ...
        loc_tern(CFG.t54_anchor_off, 'OFF', 'ON'));
loc_t54_scope();
fprintf(['  SEEDS (L-26): S6 set = seed_base + %d + (1:%d); large-sample set = ' ...
         'seed_base + %d + (1:%d). The same seeds serve every stage, leg and ' ...
         'bandwidth: stages and cells are paired, not independent replicates.\n'], ...
        CFG.t54_s6_offset, nd, CFG.t54_large_offset, nl);

for b = 1:numel(B_list)
    [P, pinfo] = t40_P(CFG, B_list(b));
    K = struct('slab', slab, 'nd', nd, 'nl', nl, 'first', b == 1, ...
               'smoke', CFG.t54_smoke);
    GT = loc_t54_Pgates(P, pinfo, CFG, F, K, GT);
    G = {};
    fprintf('\n--- %s: S1-S5 at B = %.0f MHz, K_s = %d (pinned leg, anchor on) ---\n', ...
            slab, P.B/1e6, P.K_s);

    % S1 (FIX 1): the whitened pair actually used by the arm-B path.
    rng(1, 'twister');
    Wtest = (1/sqrt(P.M)) * exp(1j*2*pi*rand(P.M, P.N_RF));
    [~, Tw, Wt] = loc_whiten(Wtest);
    G(end+1,:) = {'S1a_whiten_identity', norm(Wt'*Wt - eye(P.N_RF), 'fro'), 1e-10}; %#ok<AGROW>

    % S1b: the dictionary map is the SAME congruence as the covariance map.
    a_test = exp(1j*2*pi*rand(P.M,1));
    G(end+1,:) = {'S1b_atom_congruence', ...
                  norm(Wt'*a_test - Tw*(Wtest'*a_test)), 1e-10}; %#ok<AGROW>

    % S2: per-subcarrier path at K_s = 1 reproduces the joint driver at K_s = 1.
    Pk = P; Pk.K_s = 1; Pk.K = 1; Pk.alpha_k_vec = P.alpha_k_vec(1);
    Pk.k_indices = P.k_indices(1); Pk.B = P.Delta_f;
    R1 = t40_trial(Pk, C5, CFG.SNR_star, CFG.seed_base + 1, 'pinned', false);
    G(end+1,:) = {'S2_Ks1_A_eq_Cshared', ...
                  abs(R1.r_hat.A - R1.r_hat.C_shared_mean), 1e-9}; %#ok<AGROW>

    % S3: combination rule round trip.
    om = 1.234 * ones(8,1); ka = 5.6e-3 * ones(8,1);
    [th_c, r_c] = loc_combine(om, ka, P, 'mean');
    G(end+1,:) = {'S3_combine_roundtrip', ...
                  abs(r_c - (pi*P.d_ant^2/P.lambda_c)*sin(th_c)^2/ka(1)), 1e-9}; %#ok<AGROW>

    % S4: alpha de-scaling round trip.
    G(end+1,:) = {'S4_alpha_descale', abs(mean(P.alpha_k_vec) - 1), 5e-3}; %#ok<AGROW>

    % S5: shared scene integrity.
    R2 = t40_trial(P, C5, CFG.SNR_star, CFG.seed_base + 2, 'pinned', false);
    G(end+1,:) = {'S5_shared_scene', abs(R2.r_true - CFG.pin_r), 1e-12}; %#ok<AGROW>

    for i = 1:size(G,1)
        p = G{i,2} < G{i,3};
        ok_harness = ok_harness && p;
        fprintf('  %-24s val=%.3e tol=%.1e  %s\n', G{i,1}, G{i,2}, G{i,3}, ...
                loc_tern(p, 'PASS', 'FAIL'));
        Gall(end+1,:) = {slab, P.B, P.K_s, 'pinned', 'gate geometry', ...
                         G{i,1}, G{i,2}, G{i,3}, loc_tern(p, 'PASS', 'FAIL')}; %#ok<AGROW>
    end

    % S6 (FIX 5) and the T-54 reporting statistics, per leg (Spec Sec. 4.4).
    % S6 is a capability check on the EXPERIMENT, not on the code: if all K_s
    % estimates are bit-identical, the rival arm receives no data diversity.
    % In T-54 an S6 FAIL is the diagnostic observation, not a halt, and its
    % cause is not asserted (F-091 is contested by F-127).
    for g = 1:numel(CFG.t54_legs)
        C  = loc_t54_ctx(CFG, F, slab, P, CFG.t54_legs{g});
        st = t54_cell(P, CFG, F, C, nd, nl);
        n_cells = n_cells + 1;
        n_s6    = n_s6  + st.n_s6_lines;
        n_sel   = n_sel + st.n_sel_lines;
        n_inv   = n_inv + st.s6.n_invalid + st.large.n_invalid;
        n_tr    = n_tr  + st.s6.n + st.large.n;
        p6 = st.s6.degen_C_mean < CFG.t54_s6_tol;
        Gall(end+1,:) = {slab, P.B, P.K_s, C.geom, C.leg_label_csv, ...
                         'S6_C_not_degenerate', st.s6.degen_C_mean, ...
                         CFG.t54_s6_tol, loc_tern(p6, 'PASS', 'FAIL')}; %#ok<AGROW>
        if C.is_reg
            REG.s6 = st.s6;  REG.large = st.large;
        end
    end
end

% ---- stage verdict lines: each prints its own verdict (L-36) and names its
% ---- comparison domain (L-43) ---------------------------------------------
fprintf('\n--- %s: stage verdict lines ---\n', slab);
fprintf(['  GATE P2 [%s; S1-S5, harness self-consistency, %d bandwidth(s); a test ' ...
         'passes when its value is below its tolerance]: %s\n'], ...
        slab, numel(B_list), loc_tern(ok_harness, 'PASS', 'FAIL -- do not proceed'));
GT = loc_t54_gate(GT, 'GATE-P2-S1toS5', slab, 'all bandwidths of the stage', ...
                  'value below tolerance', double(~ok_harness), 0, ok_harness, 'yes', '');
if ~ok_harness; halt{end+1} = [slab ': GATE P2 (S1-S5) FAIL']; end

n_exp = numel(B_list) * numel(CFG.t54_legs);
ok5   = (n_cells == n_exp) && (n_s6 == n_exp) && (n_sel == 2 * n_exp);
fprintf(['  G-T5 [%s; integer domain; non-halting]: cells run %d of %d; S6 PASS/FAIL ' ...
         'lines printed %d; Phase D selection-rate blocks printed %d (two trial sets ' ...
         'per cell): %s\n'], ...
        slab, n_cells, n_exp, n_s6, n_sel, loc_tern(ok5, 'PASS', 'FAIL'));
GT = loc_t54_gate(GT, 'G-T5', slab, 'all cells of the stage', 'integer', ...
                  n_s6, n_exp, ok5, 'no', '');
fprintf(['  S6 NOTE (T-54): an S6 FAIL is the diagnostic observation of T-54, not a ' ...
         'halt; its cause is not asserted here (F-091 is contested by F-127).\n']);

fprintf(['  VALIDITY-ALL [%s; integer domain]: trials with a validity flag set = %d of ' ...
         '%d, over all cells and both trial sets (%s): %s\n'], ...
        slab, n_inv, n_tr, ...
        loc_tern(CFG.t54_smoke, 'halting in the smoke pass', 'reported; not halting'), ...
        loc_tern(n_inv == 0, 'PASS', 'FAIL'));
GT = loc_t54_gate(GT, 'VALIDITY-ALL', slab, 'all cells; both trial sets', 'integer', ...
                  n_inv, 0, n_inv == 0, loc_tern(CFG.t54_smoke, 'yes', 'no'), '');
if CFG.t54_smoke && n_inv > 0
    halt{end+1} = [slab ': invalid trial in the smoke pass'];
end

if isempty(REG.s6)
    fprintf('  REGISTERED CELL [%s]: NOT RUN (leg %s, B = %.0f MHz): FAIL\n', ...
            slab, CFG.t54_reg_leg, CFG.t54_reg_B/1e6);
    GT = loc_t54_gate(GT, 'REGISTERED-CELL', slab, 'registered cell', 'integer', ...
                      0, 1, false, 'yes', 'cell not run');
    halt{end+1} = [slab ': registered cell not run'];
else
    s6 = REG.s6;  sl = REG.large;
    okv = (s6.n_invalid == 0);
    fprintf(['  VALIDITY-REG [%s | registered cell; integer domain]: S6 trials with a ' ...
             'validity flag set = %d of %d (large-sample set: %d of %d): %s\n'], ...
            slab, s6.n_invalid, s6.n, sl.n_invalid, sl.n, loc_tern(okv, 'PASS', 'FAIL'));
    GT = loc_t54_gate(GT, 'VALIDITY-REG', slab, 'registered cell; S6 set', 'integer', ...
                      s6.n_invalid, 0, okv, loc_tern(CFG.t54_invalid_halt, 'yes', 'no'), '');
    if ~okv && CFG.t54_invalid_halt
        halt{end+1} = [slab ': VALIDITY-REG FAIL'];
    end

    if ~CFG.use_refined && ~CFG.t54_anchor_off      % stage S-1 (and its smoke pass)
        d1 = abs(s6.degen_C_mean - 1.000);
        switch CFG.t54_gt2_mode
            case 'exact'
                ok2 = (d1 <= 1e-12);
                h2 = 'yes';  rule = 'PASS iff abs(S6 - 1.000) <= 1e-12';
            case 'threshold'
                ok2 = (s6.degen_C_mean >= CFG.t54_s6_tol);
                h2 = 'yes';  rule = 'PASS iff S6 >= 0.50';
            case 'report'
                ok2 = true;
                h2 = 'no';   rule = 'reported only';
            otherwise
                error('t40_e2b_driver:T54Config', ...
                      'CFG.t54_gt2_mode = %s is not a registered mode.', CFG.t54_gt2_mode);
        end
        if CFG.t54_smoke
            v2 = 'NOT EVALUATED (smoke pass)';  h2 = 'no';
        elseif strcmp(CFG.t54_gt2_mode, 'report')
            v2 = 'REPORTED';
        else
            v2 = loc_tern(ok2, 'PASS', 'FAIL');
        end
        fprintf(['  G-T2 [%s | registered cell; integer domain, k of nd]: legacy S6 = ' ...
                 '%d/%d = %.6f | mode %s: %s | F-091 reference 1.000 (run of 2026-09-03, ' ...
                 'measured before the F-088 patch reached the pinned generator; ' ...
                 'historical comparison), abs difference %.3e: %s\n'], ...
                slab, s6.k_degen_C, s6.n, s6.degen_C_mean, CFG.t54_gt2_mode, rule, d1, v2);
        GT = loc_t54_gate(GT, 'G-T2', slab, 'registered cell; S6 set', 'integer', ...
                          s6.degen_C_mean, 1, v2, h2, rule);
        if strcmp(v2, 'FAIL'); halt{end+1} = [slab ': G-T2 FAIL']; end
    end

    fprintf(['  Q-INPUT [%s | registered cell: leg %s, B = %.0f MHz, SNR = %g dB]: ' ...
             'S6 = %d/%d = %.4f (tolerance %.2f: %s) | large-sample degen_C = %d/%d = ' ...
             '%.4f, Wilson 95%% [%.4f, %.4f] | Phase D B+C rate, arm C_shared, per ' ...
             'call: large-sample set %.4f, S6 set %.4f | arm A, per trial: ' ...
             'large-sample set %.4f | degen_Cp: S6 set %d/%d, large-sample set %d/%d\n'], ...
            slab, CFG.t54_reg_leg, CFG.t54_reg_B/1e6, CFG.SNR_star, ...
            s6.k_degen_C, s6.n, s6.degen_C_mean, CFG.t54_s6_tol, ...
            loc_tern(s6.degen_C_mean < CFG.t54_s6_tol, 'PASS', 'FAIL'), ...
            sl.k_degen_C, sl.n, sl.degen_C_mean, sl.degen_C_wilson_lo, sl.degen_C_wilson_hi, ...
            sl.selC_BC_rate_calls, s6.selC_BC_rate_calls, sl.selA_BC_rate_trials, ...
            s6.k_degen_Cp, s6.n, sl.k_degen_Cp, sl.n);
    fprintf(['  Q-INPUT note: the harness prints the inputs only. The outcome Q1 / Q2 / ' ...
             'Q3 / UNRESOLVED is classified by the scoring session, with the rule ' ...
             'Addendum D registers.\n']);
end

Tt = cell2table(Gall, 'VariableNames', {'stage_label', 'B_hz', 'K_s', 'leg', ...
     'leg_label', 'test', 'value', 'tol', 'verdict'});
Tt.run_tag    = repmat(string(F.run_tag), size(Gall, 1), 1);
Tt.scope_note = repmat(string(loc_t54_scope_text()), size(Gall, 1), 1);
loc_t54_write(Tt, loc_t54_file(F, 'tests'));
end


% =========================================================================
%  SWEEP ENGINE
% =========================================================================
function t40_sweep(CFG, out_dir, stamp, gtag, sweep_type, sweep_vec, N_MC, ...
                   flat_alpha, geom, tag)

trial_csv = fullfile(out_dir, ...
    sprintf('t40_%s_%s_trials_%s.csv',  tag, gtag, stamp));
summ_csv  = fullfile(out_dir, ...
    sprintf('t40_%s_%s_summary_%s.csv', tag, gtag, stamp));
arms  = t40_armlist(CFG);
nA    = numel(arms);
first = true;

for s = 1:numel(sweep_vec)

    if strcmp(sweep_type, 'snr')
        P = t40_P(CFG, 400e6); SNR_dB = sweep_vec(s);
    else
        P = t40_P(CFG, sweep_vec(s)); SNR_dB = CFG.SNR_star;
    end
    if flat_alpha
        P.alpha_k_vec = ones(P.K_s, 1);
    end

    fprintf('\n[T-40] %s point %d/%d  (%s = %.4g, SNR = %g dB, K_s = %d)\n', ...
            tag, s, numel(sweep_vec), sweep_type, sweep_vec(s), SNR_dB, P.K_s);

    ERR_R  = nan(N_MC,nA); ERR_TH = nan(N_MC,nA); RH  = nan(N_MC,nA);
    NIT    = nan(N_MC,nA); CONV   = nan(N_MC,nA);
    GROSS  = nan(N_MC,nA); BND    = nan(N_MC,nA);
    GLEG   = nan(N_MC,nA); OUTB   = nan(N_MC,nA);
    DISP_O = nan(N_MC,1);  DISP_K = nan(N_MC,1);  DISP_OB = nan(N_MC,1);
    DEGC   = nan(N_MC,1);  DEGB   = nan(N_MC,1);
    DEGCP  = nan(N_MC,1);                          % T54-3 (Addendum A.2)
    KARG   = nan(N_MC,1);  CONDB  = nan(N_MC,1);
    RT     = nan(N_MC,1);  TT     = nan(N_MC,1);
    RTIME  = nan(N_MC,5);

    parfor mc = 1:N_MC
        R = t40_trial(P, CFG, SNR_dB, CFG.seed_base + mc, geom, false); %#ok<PFBNS>
        [er, et, rh, ni, cv, gr, bd, gl, ob] = loc_pack(R, arms);
        ERR_R(mc,:) = er; ERR_TH(mc,:) = et; RH(mc,:) = rh;
        NIT(mc,:)   = ni; CONV(mc,:)   = cv;
        GROSS(mc,:) = gr; BND(mc,:)    = bd;
        GLEG(mc,:)  = gl; OUTB(mc,:)   = ob;
        DISP_O(mc)  = R.disp_omega;   DISP_K(mc) = R.disp_kappa;
        DISP_OB(mc) = R.disp_omega_B; DEGC(mc)   = R.degen_C;
        DEGB(mc)    = R.degen_B;      KARG(mc)   = R.k_argmin_L;
        CONDB(mc)   = R.cond_B_max;
        DEGCP(mc)   = R.degen_Cp;
        RTIME(mc,:) = [R.rt_gen R.rt_bpd R.rt_A R.rt_C R.rt_B];
        RT(mc) = R.r_true; TT(mc) = R.theta_true;
    end

    % ---- per-trial CSV --------------------------------------------------
    Ttr = table();
    for a = 1:nA
        Ttr = [Ttr; table( ...
            repmat(string(tag),N_MC,1), repmat(string(gtag),N_MC,1), ...
            repmat(sweep_vec(s),N_MC,1), repmat(SNR_dB,N_MC,1), ...
            repmat(P.B,N_MC,1), repmat(P.K_s,N_MC,1), ...
            repmat(string(geom),N_MC,1), repmat(double(flat_alpha),N_MC,1), ...
            (1:N_MC).', (CFG.seed_base+(1:N_MC)).', ...
            repmat(string(arms{a}),N_MC,1), ...
            TT*180/pi, RT, RH(:,a), ERR_R(:,a), ERR_TH(:,a)*180/pi, ...
            NIT(:,a), CONV(:,a), GROSS(:,a), GLEG(:,a), OUTB(:,a), BND(:,a), ...
            DISP_O, DISP_K, DISP_OB, DEGC, DEGB, KARG, CONDB, DEGCP, ...
            'VariableNames', {'tag','grid','sweep_value','SNR_dB','B_hz', ...
            'K_s','geom','flat_alpha','trial','seed','arm','theta_true_deg', ...
            'r_true_m','r_hat_m','err_r_m','err_theta_deg','n_iter', ...
            'converged','gross','gross_legacy','outbox','boundary', ...
            'disp_omega','disp_kappa','disp_omega_B','degen_C','degen_B', ...
            'k_argmin_L','cond_B_max','degen_Cp'})]; %#ok<AGROW>
    end
    writetable(Ttr, trial_csv, 'WriteMode', loc_tern(first,'overwrite','append'));

    % ---- summary: FIX 7 / FIX 8 two-component gate statistics -----------
    iA = find(strcmp(arms,'A'), 1);
    Ts = table();
    for a = 1:nA
        % PRIMARY precision statistic: CRMSE on each arm's OWN inliers.
        [dloc, llo, lhi, nA_in, nR_in] = loc_boot_crmse( ...
            ERR_R(:,iA), ERR_R(:,a), GROSS(:,iA), GROSS(:,a), ...
            CFG.n_boot, CFG.boot_seed);
        % SECONDARY, clearly labelled: paired common-success estimand.
        [dpair, plo, phi, npair] = loc_boot_crmse_paired( ...
            ERR_R(:,iA), ERR_R(:,a), GROSS(:,iA), GROSS(:,a), ...
            CFG.n_boot, CFG.boot_seed);
        % SECONDARY: median, and all-trial RMSE. Reported, never gated.
        [dmed, mlo, mhi, ~] = loc_paired_boot_med( ...
            ERR_R(:,iA), ERR_R(:,a), GROSS(:,iA), GROSS(:,a), ...
            CFG.n_boot, CFG.boot_seed);
        [drms, rlo, rhi] = loc_paired_boot_rmse(ERR_R(:,iA), ERR_R(:,a), ...
            CFG.n_boot, CFG.boot_seed);
        % RELIABILITY: paired McNemar risk difference (arm minus A).
        [rd, rdlo, rdhi, n10, n01] = loc_mcnemar(GROSS(:,iA), GROSS(:,a));
        [glo, ghi] = loc_wilson(sum(GROSS(:,a)==1), sum(isfinite(GROSS(:,a))));
        clean  = GROSS(:,a) == 0;
        crmse  = sqrt(mean(ERR_R(clean,a).^2, 'omitnan'));
        Ts = [Ts; table(string(tag), string(gtag), sweep_vec(s), SNR_dB, ...
            P.B, P.K_s, string(geom), double(flat_alpha), N_MC, ...
            string(arms{a}), ...
            sqrt(mean(ERR_R(:,a).^2,'omitnan')), ...
            sqrt(mean(ERR_TH(:,a).^2,'omitnan'))*180/pi, ...
            crmse, sum(clean), median(abs(ERR_R(clean,a)),'omitnan'), ...
            dloc, llo, lhi, nA_in, nR_in, ...
            dpair, plo, phi, npair, ...
            dmed, mlo, mhi, drms, rlo, rhi, ...
            mean(GROSS(:,a),'omitnan'), glo, ghi, ...
            mean(GLEG(:,a),'omitnan'), mean(OUTB(:,a),'omitnan'), ...
            rd, rdlo, rdhi, n10, n01, ...
            mean(BND(:,a),'omitnan'), mean(CONV(:,a),'omitnan'), ...
            median(NIT(:,a),'omitnan'), ...
            mean(DEGC,'omitnan'), mean(DEGB,'omitnan'), ...
            mean(RTIME(:),'omitnan'), mean(DEGCP,'omitnan'), ...
            'VariableNames', {'tag','grid','sweep_value','SNR_dB','B_hz', ...
            'K_s','geom','flat_alpha','N_MC','arm','RMSE_r_m', ...
            'RMSE_theta_deg','CRMSE_own_m','n_own_inlier','med_abs_err_m', ...
            'delta_local_dB','delta_local_lo','delta_local_hi', ...
            'n_inlier_A','n_inlier_R', ...
            'delta_paired_dB','delta_paired_lo','delta_paired_hi', ...
            'n_pair_clean','delta_med_dB','delta_med_lo','delta_med_hi', ...
            'delta_rmse_dB','delta_rmse_lo','delta_rmse_hi', ...
            'gross_rate','gross_lo','gross_hi','gross_legacy_rate', ...
            'outbox_rate','risk_diff','rd_lo','rd_hi','n10_AonlyGross', ...
            'n01_RonlyGross','clamp_rate','conv_rate','n_iter_med', ...
            'degen_C_frac','degen_B_frac','rt_mean_s','degen_Cp_frac'})]; %#ok<AGROW>
    end
    writetable(Ts, summ_csv, 'WriteMode', loc_tern(first,'overwrite','append'));
    first = false;

    fprintf('  %-18s %9s %9s %21s %22s\n', 'arm', 'RMSE_all', 'CRMSE', ...
            'Delta_local [CI] dB', 'riskdiff vs A [CI]');
    for a = 1:nA
        fprintf('  %-18s %9.5f %9.5f  %+6.2f [%+6.2f,%+6.2f]  %+.4f [%+.4f,%+.4f]\n', ...
            arms{a}, Ts.RMSE_r_m(a), Ts.CRMSE_own_m(a), ...
            Ts.delta_local_dB(a), Ts.delta_local_lo(a), Ts.delta_local_hi(a), ...
            Ts.risk_diff(a), Ts.rd_lo(a), Ts.rd_hi(a));
    end
    fprintf('  gross rates (box OR legacy): ');
    for a = 1:nA; fprintf('%s=%.4f ', arms{a}, Ts.gross_rate(a)); end
    fprintf('\n');
    fprintf('  degeneracy (bit-identical per-subcarrier estimates): C %.3f | B %.3f\n', ...
            mean(DEGC,'omitnan'), mean(DEGB,'omitnan'));
    fprintf('  degeneracy, arm C_percarrier (control; reported, never gated; not S6): %.3f\n', ...
            mean(DEGCP,'omitnan'));
end
fprintf('\n  -> %s\n  -> %s\n', trial_csv, summ_csv);
end


function arms = t40_armlist(CFG)
arms = {'A', 'C_shared_mean', 'C_shared_med', 'C_shared_trim', ...
        'C_shared_klsel', 'B_matched_mean', 'B_matched_med', ...
        'B_matched_trim', 'B_naive_mean'};
if CFG.do_C_percarrier
    arms = [arms, {'C_percarrier_mean'}];
end
end


% =========================================================================
%  ONE TRIAL -- all arms on identical data
% =========================================================================
function R = t40_trial(P, CFG, SNR_dB, seed, geom, only_A)

if nargin < 6 || isempty(only_A); only_A = false; end
rng(seed, 'twister');

t0 = tic;
if strcmp(geom, 'pinned')
    [X, Y, theta_t, r_t, N0, W] = loc_gen_fixed(P, SNR_dB, ...
                                    CFG.pin_theta, CFG.pin_r, seed);
elseif strcmp(geom, 'pinned_far')
    [X, Y, theta_t, r_t, N0, W] = loc_gen_fixed(P, SNR_dB, ...
                                    CFG.pin_theta, CFG.pin_r_far, seed);
else
    [X, Y, ~, theta_t, r_t, ~, N0, W] = wb_channel_gen_ofdm_nf(P, SNR_dB);
    theta_t = theta_t(1); r_t = r_t(1);
end
R.rt_gen = toc(t0);
R.theta_true = theta_t; R.r_true = r_t; R.N0 = N0;

K_s = P.K_s;
Rh  = cell(K_s,1);
for k = 1:K_s
    Yk = Y(:,:,k);  Rk = (Yk*Yk')/P.N;  Rh{k} = (Rk + Rk')/2;
end

% ---- shared wideband BPD warm start --------------------------------------
% T-54 (validity): a failed BPD still falls back to the central anchor, as
% before, but the trial is flagged (bpd_fail) and the message is kept.
t0 = tic;
bpd_fail = 0;  msg1 = '';
try
    [th_bpd, r_bpd] = bpd_baseline(X, P);
    th_bpd = th_bpd(1); r_bpd = r_bpd(1);
catch ME_w
    th_bpd = (P.theta_lo + P.theta_hi)/2;
    r_bpd  = 1/(0.5*(P.u_min + P.u_max));
    bpd_fail = 1;
    msg1 = loc_t54_first(msg1, 'bpd_baseline (wideband)', ME_w.message);
end
R.rt_bpd = toc(t0);
p0 = ones(P.d,1)/P.d;

% ---- ARM A ---------------------------------------------------------------
t0 = tic;
PA = loc_t54_P(P, CFG, 'A');
[thA, rA, ~, niA, cvA, selA, LdA, exA, emA] = ...
    loc_safe_clkl(Rh, W, th_bpd, r_bpd, p0, PA, CFG.t54_prepatch);
msg1 = loc_t54_first(msg1, 'wb_clkl_driver_pc (arm A)', emA);
R.rt_A = toc(t0);

r_hat = struct(); th_hat = struct(); n_it = struct(); conv = struct();
r_hat.A = rA; th_hat.A = thA; n_it.A = niA; conv.A = cvA;

if only_A
    R.rt_C = 0; R.rt_B = 0;
    R.disp_omega = NaN; R.disp_kappa = NaN; R.disp_omega_B = NaN;
    R.degen_C = NaN; R.degen_B = NaN; R.k_argmin_L = NaN; R.cond_B_max = NaN;
    R.degen_Cp = NaN;
    R = loc_finish(R, r_hat, th_hat, n_it, conv, P, CFG);
    return
end

% ---- ARMS C -------------------------------------------------------------
t0 = tic;
omC  = nan(K_s,1); kaC  = nan(K_s,1); LkC = nan(K_s,1);
niC  = nan(K_s,1); cvC  = nan(K_s,1);
omCp = nan(K_s,1); kaCp = nan(K_s,1); niCp = nan(K_s,1); cvCp = nan(K_s,1);
% T-54 change T54-1: per-subcarrier Phase D log (selection code, L_PhaseD)
% and validity flags (absorbed error; failed per-carrier BPD).
thC  = nan(K_s,1); rC  = nan(K_s,1); selC  = zeros(K_s,1); LdC  = nan(K_s,3);
thCp = nan(K_s,1); rCp = nan(K_s,1); selCp = zeros(K_s,1); LdCp = nan(K_s,3);
exC  = zeros(K_s,1); exCp = zeros(K_s,1);
tbv  = nan(K_s,1); rbv = nan(K_s,1); bpdk_fail = zeros(K_s,1);
for k = 1:K_s
    Pk  = loc_Pk(P, k);
    PkC = loc_t54_P(Pk, CFG, 'C_shared');
    [thk, rk, Lk, nik, cvk, sk, Ldk, exk, emk] = ...
        loc_safe_clkl(Rh(k), W, th_bpd, r_bpd, p0, PkC, CFG.t54_prepatch);
    [omC(k), kaC(k)] = loc_to_natural(thk, rk, P);
    LkC(k) = Lk; niC(k) = nik; cvC(k) = cvk;
    thC(k) = thk; rC(k) = rk; selC(k) = loc_t54_selcode(sk);
    LdC(k,:) = Ldk; exC(k) = exk;
    msg1 = loc_t54_first(msg1, sprintf('wb_clkl_driver_pc (C_shared, k = %d)', k), emk);
    if CFG.do_C_percarrier
        try
            [tb, rb] = bpd_baseline(X(:,:,k), Pk);
            tb = tb(1); rb = rb(1);
        catch ME_k
            % The fallback to the SHARED anchor is kept, but it defeats the
            % C_percarrier control (Addendum A.2): the trial is flagged.
            tb = th_bpd; rb = r_bpd;
            bpdk_fail(k) = 1;
            msg1 = loc_t54_first(msg1, sprintf('bpd_baseline (subcarrier %d)', k), ...
                                 ME_k.message);
        end
        tbv(k) = tb; rbv(k) = rb;
        PkCp = loc_t54_P(Pk, CFG, 'C_percarrier');
        [thk2, rk2, ~, ni2, cv2, sk2, Ldk2, exk2, emk2] = ...
            loc_safe_clkl(Rh(k), W, tb, rb, p0, PkCp, CFG.t54_prepatch);
        [omCp(k), kaCp(k)] = loc_to_natural(thk2, rk2, P);
        niCp(k) = ni2; cvCp(k) = cv2;
        thCp(k) = thk2; rCp(k) = rk2; selCp(k) = loc_t54_selcode(sk2);
        LdCp(k,:) = Ldk2; exCp(k) = exk2;
        msg1 = loc_t54_first(msg1, sprintf('wb_clkl_driver_pc (C_percarrier, k = %d)', k), emk2);
    end
end
R.rt_C = toc(t0);

% ---- ARMS B (FIX 1: one whitened pair, conjugate transpose throughout) ---
t0 = tic;
[~, Tw, Wt] = loc_whiten(W);
omBm = nan(K_s,1); kaBm = nan(K_s,1); niB = nan(K_s,1); cvB = nan(K_s,1);
omBn = nan(K_s,1); kaBn = nan(K_s,1);
LkB  = nan(K_s,1); condB = nan(K_s,1);
for k = 1:K_s
    Rt = Tw * Rh{k} * Tw';  Rt = (Rt + Rt')/2;
    condB(k) = cond(Rt);
    Pb = loc_Pb(P, P.lambda_c / P.alpha_k_vec(k));
    [tb, rb, Lb, nib, cvb] = loc_safe_nf(Rt, Wt, Pb);
    [omBm(k), kaBm(k)] = loc_to_natural(tb, rb, P);
    LkB(k) = Lb; niB(k) = nib; cvB(k) = cvb;

    Pb2 = loc_Pb(P, P.lambda_c);
    [tb2, rb2] = loc_safe_nf(Rt, Wt, Pb2);
    [o2, k2] = loc_to_natural(tb2, rb2, P);
    omBn(k) = o2 / P.alpha_k_vec(k);
    kaBn(k) = k2 / P.alpha_k_vec(k);
end
R.rt_B = toc(t0);

% ---- combination ---------------------------------------------------------
sets = {'C_shared',  omC,  kaC,  LkC, niC, cvC; ...
        'B_matched', omBm, kaBm, LkB, niB, cvB};
for i = 1:size(sets,1)
    nm = sets{i,1}; om = sets{i,2}; ka = sets{i,3};
    Lv = sets{i,4}; nv = sets{i,5}; cvv = sets{i,6};
    for rule = {'mean','med','trim'}
        [t_, r_] = loc_combine(om, ka, P, rule{1});
        th_hat.([nm '_' rule{1}]) = t_;  r_hat.([nm '_' rule{1}]) = r_;
        n_it.([nm '_' rule{1}])   = median(nv,'omitnan');
        conv.([nm '_' rule{1}])   = mean(cvv,'omitnan');
    end
    [~, kbest] = min(Lv);
    if isempty(kbest) || isnan(kbest); kbest = 1; end
    [t_, r_] = loc_combine(om(kbest), ka(kbest), P, 'mean');
    th_hat.([nm '_klsel']) = t_; r_hat.([nm '_klsel']) = r_;
    n_it.([nm '_klsel'])   = nv(kbest);
    conv.([nm '_klsel'])   = cvv(kbest);
end
[t_, r_] = loc_combine(omBn, kaBn, P, 'mean');
th_hat.B_naive_mean = t_; r_hat.B_naive_mean = r_;
n_it.B_naive_mean = median(niB,'omitnan');
conv.B_naive_mean = mean(cvB,'omitnan');

if CFG.do_C_percarrier
    [t_, r_] = loc_combine(omCp, kaCp, P, 'mean');
    th_hat.C_percarrier_mean = t_; r_hat.C_percarrier_mean = r_;
    n_it.C_percarrier_mean = median(niCp,'omitnan');
    conv.C_percarrier_mean = mean(cvCp,'omitnan');
end

% ---- FIX 5: degeneracy diagnostics --------------------------------------
R.disp_omega   = std(omC,  'omitnan');
R.disp_kappa   = std(kaC,  'omitnan');
R.disp_omega_B = std(omBm, 'omitnan');
R.degen_C = double(loc_alleq(omC)  && loc_alleq(kaC));
R.degen_B = double(loc_alleq(omBm) && loc_alleq(kaBm));
% T-54 change T54-3 (Addendum A.2): the same indicator for arm C_percarrier.
% A control statistic: REPORTED, NEVER GATED. It is not S6, does not enter
% S6, and amends neither S6's nd = 12 definition nor its 0.50 tolerance.
if CFG.do_C_percarrier
    R.degen_Cp = double(loc_alleq(omCp) && loc_alleq(kaCp));
else
    R.degen_Cp = NaN;
end
[~, R.k_argmin_L] = min(LkC);
R.cond_B_max = max(condB);

% ---- T-54 change T54-1: Phase D log and validity flags of this trial ------
% Selections are stored as codes (loc_t54_selcode): 0 '?', 1 'A', 2 'B',
% 3 'C', 4 'a'. Read by t54_one only; no T-40 stage uses R.t54.
T54 = struct();
T54.th_bpd = th_bpd;  T54.r_bpd = r_bpd;  T54.bpd_fail = bpd_fail;
T54.thA = thA;  T54.rA = rA;  T54.selA = loc_t54_selcode(selA);
T54.LdA = LdA;  T54.exA = exA;
T54.thC = thC;  T54.rC = rC;  T54.omC = omC;  T54.kaC = kaC;
T54.selC = selC;  T54.LdC = LdC;  T54.exC = exC;  T54.niC = niC;  T54.cvC = cvC;
T54.tbv = tbv;  T54.rbv = rbv;  T54.bpdk_fail = bpdk_fail;
T54.thCp = thCp;  T54.rCp = rCp;  T54.omCp = omCp;  T54.kaCp = kaCp;
T54.selCp = selCp;  T54.LdCp = LdCp;  T54.exCp = exCp;
T54.niCp = niCp;  T54.cvCp = cvCp;
T54.msg = msg1;
R.t54 = T54;

R = loc_finish(R, r_hat, th_hat, n_it, conv, P, CFG);
end


function R = loc_finish(R, r_hat, th_hat, n_it, conv, P, CFG)
%LOC_FINISH  Errors and error-class flags for every arm (FIX 2, FIX 6).
%
%  Three indicators are recorded per arm and per trial:
%    gross_legacy  |err_r|/r_true > CFG.gross_frac. The published definition.
%                  Retained for continuity. NOT a wrong-mode indicator.
%    outbox        r_hat outside the SCENE box [r_lo_fac*r_RD, r_hi_fac*r_RD].
%                  The true target is always inside that box by construction,
%                  so an estimate outside it is a mode error by definition.
%    gross         outbox OR gross_legacy. This is the PRIMARY indicator and
%                  the one the gate conditions on.
%    boundary      r_hat sits on the admissible-set bound (u_min or u_max).
R.r_hat = r_hat; R.th_hat = th_hat; R.n_iter = n_it; R.converged = conv;
r_adm_lo = 1/P.u_max;            r_adm_hi = 1/P.u_min;
r_box_lo = P.r_lo_fac * P.r_RD;  r_box_hi = P.r_hi_fac * P.r_RD;
fn = fieldnames(r_hat);
R.err_r = struct(); R.err_theta = struct();
R.gross = struct(); R.gross_legacy = struct();
R.outbox = struct(); R.boundary = struct();
for i = 1:numel(fn)
    f  = fn{i};
    rh = r_hat.(f);
    er = rh - R.r_true;
    R.err_r.(f)       = er;
    R.err_theta.(f)   = th_hat.(f) - R.theta_true;
    gl = double(abs(er)/max(R.r_true,eps) > CFG.gross_frac);
    ob = double(isfinite(rh) && (rh > r_box_hi + 1e-9 || rh < r_box_lo - 1e-9));
    R.gross_legacy.(f) = gl;
    R.outbox.(f)       = ob;
    if CFG.use_box_rule
        R.gross.(f) = double(gl == 1 || ob == 1);
    else
        R.gross.(f) = gl;
    end
    R.boundary.(f) = double(abs(rh - r_adm_lo) < 1e-3 || ...
                            abs(rh - r_adm_hi) < 1e-3);
end
end


% =========================================================================
%  LOCAL HELPERS
% =========================================================================
function Pk = loc_Pk(P, k)
Pk = P;
Pk.K_s = 1; Pk.K = 1;
Pk.alpha_k_vec = P.alpha_k_vec(k);
Pk.k_indices   = P.k_indices(k);
Pk.B = P.Delta_f;
end

function Pb = loc_Pb(P, lam)
Pb = P;
Pb.lambda   = lam;
Pb.lambda_c = lam;
Pb.nf_eps_matched = true;
end

function [th, r, L, ni, cv, sel, Ld, ex, emsg] = ...
        loc_safe_clkl(Rc, W, th0, r0, p0, P, use_prepatch)
%LOC_SAFE_CLKL  One wb_clkl_driver_pc call. T-54 change T54-1 extends the
%  return list: sel is info.PhaseD_select (one character: 'A', 'B' or 'C';
%  'a' when P.phaseD_anchor_off is true; '?' if absent), Ld is info.L_PhaseD
%  (1x3 double; NaN if absent), ex is 1 if the call threw -- the error is
%  absorbed as before, but no longer silently -- and emsg is its message.
%  use_prepatch = true (gate G-T3 only) calls the pre-patch copy instead.
th = NaN; r = NaN; L = Inf; ni = NaN; cv = 0;
sel = '?'; Ld = nan(1,3); ex = 0; emsg = '';
if nargin < 7 || isempty(use_prepatch); use_prepatch = false; end
try
    if use_prepatch
        [th, r, ~, ~, info] = wb_clkl_driver_pc_prepatch(Rc, W, th0, r0, p0, P);
    else
        [th, r, ~, ~, info] = wb_clkl_driver_pc(Rc, W, th0, r0, p0, P);
    end
    th = th(1); r = r(1);
    if isfield(info,'L_hist') && ~isempty(info.L_hist); L = info.L_hist(end); end
    if isfield(info,'n_iter');    ni = info.n_iter;              end
    if isfield(info,'converged'); cv = double(info.converged);   end
    if isfield(info,'PhaseD_select') && ischar(info.PhaseD_select) ...
            && numel(info.PhaseD_select) == 1
        sel = info.PhaseD_select;
    end
    if isfield(info,'L_PhaseD') && isnumeric(info.L_PhaseD) ...
            && numel(info.L_PhaseD) == 3
        Ld = double(reshape(info.L_PhaseD, 1, 3));
    end
catch ME
    ex = 1; emsg = ME.message;
end
end

function [th, r, L, ni, cv] = loc_safe_nf(Rt, Wt, Pb)
th = NaN; r = NaN; L = Inf; ni = NaN; cv = 0;
try
    [th, r, ~, ~, info] = nf_clkl_pc(Rt, Wt, Pb);
    th = th(1); r = r(1);
    if isfield(info,'L_hist') && ~isempty(info.L_hist); L = info.L_hist(end); end
    if isfield(info,'n_iter');    ni = info.n_iter;              end
    if isfield(info,'converged'); cv = double(info.converged);   end
catch
end
end

function [om, ka] = loc_to_natural(th, r, P)
c_lin  = 2*pi*P.d_ant   / P.lambda_c;
c_quad =   pi*P.d_ant^2 / P.lambda_c;
om = c_lin  * cos(th);
ka = c_quad * sin(th)^2 / r;
end

function [th, r] = loc_combine(om, ka, P, rule)
om = om(isfinite(om)); ka = ka(isfinite(ka));
if isempty(om) || isempty(ka); th = NaN; r = NaN; return; end
switch rule
    case 'mean';  ob = mean(om);   kb = mean(ka);
    case 'med';   ob = median(om); kb = median(ka);
    case 'trim';  ob = loc_trimmean(om,20); kb = loc_trimmean(ka,20);
    otherwise;    error('loc_combine: unknown rule %s', rule);
end
c_lin  = 2*pi*P.d_ant   / P.lambda_c;
c_quad =   pi*P.d_ant^2 / P.lambda_c;
arg = min(1-1e-9, max(1e-9, ob/c_lin));
th  = acos(arg);
u   = kb / (c_quad * sin(th)^2);
u   = min(P.u_max, max(P.u_min, u));
r   = 1/u;
end

function m = loc_trimmean(v, pct)
v = sort(v(:)); n = numel(v);
c = floor(n*pct/200);
if n - 2*c < 1; m = median(v); else; m = mean(v(1+c:n-c)); end
end

function tf = loc_alleq(v)
v = v(isfinite(v));
tf = ~isempty(v) && all(v == v(1));
end

function [Gw, Tw, Wt] = loc_whiten(W)
%LOC_WHITEN  Single definition of the whitened pair (FIX 1).
%   W'W = L L^H,  Tw = L^{-1},  Wt = W Tw'  ==>  Wt' Wt = I  and
%   Wt' a = Tw (W' a), the SAME congruence applied to the covariance by
%   Tw R Tw'. Note Tw' (conjugate transpose), never Tw.': Tw is complex.
Gw = W'*W; Gw = (Gw + Gw')/2;
L  = chol(Gw, 'lower');
Tw = L \ eye(size(L,1));
Wt = W * Tw';
end

function [d_dB, lo, hi, nA_in, nR_in] = loc_boot_crmse(eA, eR, gA, gR, nboot, seed)
%LOC_BOOT_CRMSE  PRIMARY precision statistic (FIX 7).
%   Conditional RMSE on each arm's OWN inlier set, so neither arm's score
%   depends on the other arm's failures (benchmark independence). The
%   bootstrap resamples the COMMON trial index, retaining the pairing in the
%   data while keeping the conditioning per-arm.
%   d_dB = 20 log10( CRMSE(rival) / CRMSE(A) ).  Positive => A is better.
ok = isfinite(eA) & isfinite(eR) & isfinite(gA) & isfinite(gR);
eA = eA(ok); eR = eR(ok); gA = gA(ok); gR = gR(ok);
n  = numel(eA);
iA = gA == 0; iR = gR == 0;
nA_in = sum(iA); nR_in = sum(iR);
if nA_in < 20 || nR_in < 20
    d_dB = NaN; lo = NaN; hi = NaN; return
end
cA = sqrt(mean(eA(iA).^2));  cR = sqrt(mean(eR(iR).^2));
d_dB = 20*log10(cR / max(cA, realmin));
rs = RandStream('twister','Seed',seed);
b  = nan(nboot,1);
for i = 1:nboot
    idx = randi(rs, n, n, 1);
    sA  = eA(idx); sR = eR(idx); hA = gA(idx); hR = gR(idx);
    jA  = hA == 0; jR = hR == 0;
    if sum(jA) < 5 || sum(jR) < 5; continue; end
    b(i) = 20*log10(sqrt(mean(sR(jR).^2)) / max(sqrt(mean(sA(jA).^2)), realmin));
end
b = sort(b(isfinite(b)));
if numel(b) < 50; lo = NaN; hi = NaN; return; end
lo = b(max(1, round(0.025*numel(b))));
hi = b(min(numel(b), round(0.975*numel(b))));
end

function [d_dB, lo, hi, npair] = loc_boot_crmse_paired(eA, eR, gA, gR, nboot, seed)
%LOC_BOOT_CRMSE_PAIRED  SECONDARY, clearly labelled (FIX 7).
%   Common-success conditional estimand: the CRMSE ratio on trials where BOTH
%   arms are non-gross. This is NOT the marginal performance difference. It is
%   reported alongside loc_boot_crmse and is never the gate.
ok = isfinite(eA) & isfinite(eR) & (gA == 0) & (gR == 0);
eA = eA(ok); eR = eR(ok); npair = numel(eA);
if npair < 20; d_dB = NaN; lo = NaN; hi = NaN; return; end
d_dB = 20*log10(sqrt(mean(eR.^2)) / max(sqrt(mean(eA.^2)), realmin));
rs = RandStream('twister','Seed',seed);
b  = zeros(nboot,1);
for i = 1:nboot
    idx  = randi(rs, npair, npair, 1);
    b(i) = 20*log10(sqrt(mean(eR(idx).^2)) / max(sqrt(mean(eA(idx).^2)), realmin));
end
bs = sort(b);
lo = bs(max(1, round(0.025*nboot)));
hi = bs(min(nboot, round(0.975*nboot)));
end

function [rd, lo, hi, n10, n01] = loc_mcnemar(gA, gR)
%LOC_MCNEMAR  PAIRED reliability comparison (FIX 8).
%   The arms see identical data, so global-error probabilities are compared on
%   the discordant pairs, not by two marginal Wilson intervals. Returns the
%   paired risk difference rd = p_R - p_A with a 95 percent Wald interval on
%   the discordant counts, plus the raw counts.
%   rd > 0 means the RIVAL fails more often than arm A.
ok  = isfinite(gA) & isfinite(gR);
gA  = gA(ok) == 1; gR = gR(ok) == 1;
n   = numel(gA);
n10 = sum(gA & ~gR);      % A gross, rival not
n01 = sum(~gA & gR);      % rival gross, A not
if n == 0; rd = NaN; lo = NaN; hi = NaN; return; end
rd = (n01 - n10) / n;
z  = 1.959963984540054;
se = sqrt(max(n01 + n10 - (n01 - n10)^2 / n, 0)) / n;
if (n01 + n10) == 0
    lo = 0; hi = 0;       % no discordant pairs: no evidence either way
else
    lo = max(-1, rd - z*se);
    hi = min( 1, rd + z*se);
end
end

function [d_dB, lo, hi, npair] = loc_paired_boot_med(eA, eR, gA, gR, nboot, seed)
%LOC_PAIRED_BOOT_MED  SECONDARY summary only; superseded as the gate by
%   loc_boot_crmse (FIX 7). Median absolute range error on trials where BOTH
%   arms are non-gross. Retained because it is nearly insensitive to the rare
%   mode and is therefore a useful cross-check on the bulk.
ok = isfinite(eA) & isfinite(eR) & (gA == 0) & (gR == 0);
eA = abs(eA(ok)); eR = abs(eR(ok)); npair = numel(eA);
if npair < 20; d_dB = NaN; lo = NaN; hi = NaN; return; end
d_dB = 20*log10(median(eR)/max(median(eA), realmin));
rs = RandStream('twister','Seed',seed);
b  = zeros(nboot,1);
for i = 1:nboot
    idx  = randi(rs, npair, npair, 1);
    b(i) = 20*log10(median(eR(idx))/max(median(eA(idx)), realmin));
end
bs = sort(b);
lo = bs(max(1, round(0.025*nboot)));
hi = bs(min(nboot, round(0.975*nboot)));
end

function [d_dB, lo, hi] = loc_paired_boot_rmse(eA, eR, nboot, seed)
%LOC_PAIRED_BOOT_RMSE  Reported, not gated (FIX 4).
ok = isfinite(eA) & isfinite(eR);
eA = eA(ok); eR = eR(ok); n = numel(eA);
d_dB = 20*log10(sqrt(mean(eR.^2))/max(sqrt(mean(eA.^2)), realmin));
if n < 20; lo = NaN; hi = NaN; return; end
rs = RandStream('twister','Seed',seed);
b  = zeros(nboot,1);
for i = 1:nboot
    idx  = randi(rs, n, n, 1);
    b(i) = 20*log10(sqrt(mean(eR(idx).^2))/max(sqrt(mean(eA(idx).^2)), realmin));
end
bs = sort(b);
lo = bs(max(1, round(0.025*nboot)));
hi = bs(min(nboot, round(0.975*nboot)));
end

function [lo, hi] = loc_wilson(k, n)
%LOC_WILSON  95 percent Wilson score interval for a binomial proportion.
if n == 0; lo = NaN; hi = NaN; return; end
z = 1.959963984540054;
p = k/n; den = 1 + z^2/n;
c = (p + z^2/(2*n)) / den;
h = z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / den;
lo = max(0, c-h); hi = min(1, c+h);
end

function [er, et, rh, ni, cv, gr, bd, gl, ob] = loc_pack(R, arms)
%LOC_PACK  FIX 2 / FIX 6: every diagnostic column populated, three error
%          classes recorded separately.
n = numel(arms);
er = nan(1,n); et = nan(1,n); rh = nan(1,n);
ni = nan(1,n); cv = nan(1,n); gr = nan(1,n); bd = nan(1,n);
gl = nan(1,n); ob = nan(1,n);
for a = 1:n
    f = arms{a};
    if isfield(R.err_r, f)
        er(a) = R.err_r.(f);      et(a) = R.err_theta.(f);
        rh(a) = R.r_hat.(f);      gr(a) = R.gross.(f);
        bd(a) = R.boundary.(f);
        ni(a) = R.n_iter.(f);     cv(a) = R.converged.(f);
        gl(a) = R.gross_legacy.(f); ob(a) = R.outbox.(f);
    end
end
end

function s = loc_tern(c, a, b)
if c; s = a; else; s = b; end
end

% -------------------------------------------------------------------------
%  loc_gen_fixed -- pinned-geometry generator.
%  TRANSCRIBED from wb_channel_gen_ofdm_nf_fixed inside run_monte_carlo_paperC.m.
%  Author MUST diff this body against that source and record the diff.
% -------------------------------------------------------------------------
function [X_full, Y_full, theta_true, r_true, N0, W_comb] = ...
        loc_gen_fixed(P, SNR_dB, theta_fix, r_fix, seed)
rng(seed, 'twister');
P_tmp = P;
P_tmp.r_lo_fac = max(P.r_lo_fac, r_fix / P.r_RD * 0.99);
P_tmp.r_hi_fac = min(P.r_hi_fac, r_fix / P.r_RD * 1.01);
[~, ~, ~, ~, ~, ~, N0, W_comb] = wb_channel_gen_ofdm_nf(P_tmp, SNR_dB);

theta_true = theta_fix; r_true = r_fix;
p_true = ones(P.d,1)/P.d;
M = P.M; N = P.N; K_s = P.K_s;
m_bar = ((0:M-1) - (M-1)/2).';
X_full = zeros(M, N, K_s);
for k = 1:K_s
    ak  = P.alpha_k_vec(k);
    A_k = zeros(M, P.d);
    for l = 1:P.d
        A_k(:,l) = wb_nf_fresnel_steer(theta_true, 1/r_true, ak, P) * sqrt(M);
    end
    S_k = sqrt(p_true(1)) * (randn(P.d,N) + 1j*randn(P.d,N)) / sqrt(2);
    W_k = sqrt(N0)        * (randn(M,N)   + 1j*randn(M,N))   / sqrt(2);
    X_full(:,:,k) = A_k*S_k + W_k;
end
Y_full = zeros(P.N_RF, N, K_s);
for k = 1:K_s
    Y_full(:,:,k) = W_comb' * X_full(:,:,k);
end
end


% =========================================================================
%  T-54 (A6) -- finer-grid self-test of the T-40 harness. DIAGNOSTIC ONLY.
%  Spec: PaperC_T53_T54_Specs.md Sec. 4, Addendum A, Addendum D.
%  Everything below is used by the T-54 stages only.
% =========================================================================
function t54_dispatch(CFG, F, grid_mode)
%T54_DISPATCH  Runs one T-54 stage: preamble gates, the stage, the gates
%  CSV, the halting conditions. On an error it prints the message into the
%  console file, closes the diary and rethrows.
t_stage = tic;
try
    [GT, halt] = loc_t54_preamble(F, grid_mode);
    if ~isempty(halt)
        loc_t54_finish(GT, halt, F);
    end
    switch F.stage
        case 'P0'
            [GT, halt] = t54_bench(CFG, F, GT);
        case 'T54_GT3'
            [GT, halt] = t54_gt3(CFG, F, GT);
        case 'P2_SELFTEST'
            [GT, halt] = t40_selftest(CFG, F, GT, false, false);
        case 'P2_SELFTEST_NOANCHOR'
            [GT, halt] = t40_selftest(CFG, F, GT, true, false);
        otherwise
            error('t40_e2b_driver:T54Stage', '%s is not a T-54 stage.', F.stage);
    end
    fprintf('\n[T-54] stage %s (%s grid) elapsed %.1f s (wall clock; reporting only).\n', ...
            F.stage, F.gtag, toc(t_stage));
    loc_t54_scope();
    loc_t54_finish(GT, halt, F);
catch ME
    fprintf('\n[T-54] stage %s STOPPED after %.1f s: %s\n', F.stage, toc(t_stage), ME.message);
    if ~isempty(ME.stack)
        fprintf('[T-54] raised in %s, line %d\n', ME.stack(1).name, ME.stack(1).line);
    end
    diary off;
    rethrow(ME);
end
end


function [GT, halt] = t54_bench(CFG, F, GT)
%T54_BENCH  T-54 stage S-0 ('P0'; Spec Sec. 4.4, 4.7). Part 1, smoke pass:
%  the S-1..S-4 code path with a few trials, so that a coding error stops
%  the run in its first minutes. Part 2: per-trial wall-clock time at the
%  legacy and at the refined scan grid on the same trials, and the stop
%  rule. No statistic of this stage is read.
halt = {};
fprintf('\n### T-54 S-0, part 1: smoke pass (legacy grid, %d trials per set; NOT a statistic) ###\n', ...
        CFG.t54_smoke_nd);
C = CFG;  C.use_refined = false;
[GT, h1] = t40_selftest(C, F, GT, false, true);
[GT, h2] = t40_selftest(C, F, GT, true,  true);
halt = [halt, h1, h2];
ok_s = isempty(halt);
fprintf('\n  S0-SMOKE [integer domain]: halting conditions fired in the two smoke passes = %d: %s\n', ...
        numel(halt), loc_tern(ok_s, 'PASS', 'FAIL'));
GT = loc_t54_gate(GT, 'S0-SMOKE', 'S0', 'smoke pass; anchor on and off', 'integer', ...
                  numel(halt), 0, ok_s, 'yes', '');
if ~ok_s
    fprintf('  The benchmark is not run.\n');
    return
end

fprintf('\n### T-54 S-0, part 2: per-trial time at both scan grids; stop rule (Spec Sec. 4.7) ###\n');
nb    = CFG.N_MC_bench;
seeds = CFG.seed_base + (1:nb);
gm    = {'legacy'; 'refined'};
m     = nan(2, 5);
qs    = nan(2, 2);
ninv  = 0;
for g = 1:2
    C = CFG;  C.use_refined = (g == 2);
    [P, pinfo] = t40_P(C, CFG.t54_reg_B);
    K  = struct('slab', 'S0', 'nd', CFG.t54_nd, 'nl', CFG.t54_n_large, ...
                'first', g == 1, 'smoke', false);
    GT = loc_t54_Pgates(P, pinfo, C, F, K, GT);
    if g == 1
        t54_run_set(P, C, CFG.t54_reg_leg, seeds(1), false);   % warm-up, untimed
    end
    TR = t54_run_set(P, C, CFG.t54_reg_leg, seeds, false);
    m(g,:)  = [mean(loc_t54_col(TR, 'trial', 'rt_gen')), ...
               mean(loc_t54_col(TR, 'trial', 'rt_bpd')), ...
               mean(loc_t54_col(TR, 'trial', 'rt_A')), ...
               mean(loc_t54_col(TR, 'trial', 'rt_C')), ...
               mean(loc_t54_col(TR, 'trial', 'rt_B'))];
    qs(g,:) = [P.Q_scan_th, P.Q_scan_u];
    ninv    = ninv + sum(loc_t54_col(TR, 'trial', 'valid') == 0);
    fprintf(['  BENCH grid=%-7s Q_scan=[%d, %d]: gen %.3f | bpd_wb %.3f | armA %.3f | ' ...
             'armsC %.3f | armsB %.3f | total %.3f s/trial (n = %d trials, serial, ' ...
             'seeds %d..%d)\n'], ...
            gm{g}, qs(g,1), qs(g,2), m(g,1), m(g,2), m(g,3), m(g,4), m(g,5), ...
            sum(m(g,:)), nb, seeds(1), seeds(end));
end
tot   = sum(m, 2);
ratio = tot(2) / tot(1);
go    = isfinite(ratio) && (ratio <= CFG.t54_stop_ratio);
fprintf(['  STOP RULE [Spec Sec. 4.7; domain: ratio of the mean wall-clock seconds per ' ...
         'trial, all arms; leg %s, B = %.0f MHz, SNR = %g dB, n = %d trials per grid]: ' ...
         'legacy %.3f s | refined %.3f s | ratio %.3f | limit %g: %s\n'], ...
        CFG.t54_reg_leg, CFG.t54_reg_B/1e6, CFG.SNR_star, nb, tot(1), tot(2), ratio, ...
        CFG.t54_stop_ratio, loc_tern(go, 'CONTINUE', 'STOP -- do not run S-1 to S-4; report'));
GT = loc_t54_gate(GT, 'STOP-RULE', 'S0', 'benchmark; both grids', 'linear ratio', ...
                  ratio, CFG.t54_stop_ratio, loc_tern(go, 'CONTINUE', 'STOP'), 'yes', '');
if ~go; halt{end+1} = 'S0: STOP RULE (cost)'; end

fprintf('  BENCH VALIDITY [integer domain]: benchmark trials with a validity flag set = %d of %d: %s\n', ...
        ninv, 2 * nb, loc_tern(ninv == 0, 'PASS', 'FAIL'));
GT = loc_t54_gate(GT, 'BENCH-VALIDITY', 'S0', 'benchmark; both grids', 'integer', ...
                  ninv, 0, ninv == 0, 'yes', '');
if ninv > 0; halt{end+1} = 'S0: invalid trial in the benchmark'; end

n_stage = numel(CFG.t54_legs) * numel(CFG.B_grid) * (CFG.t54_nd + CFG.t54_n_large);
fprintf(['  PROJECTION (reporting only; from the %.0f MHz per-trial totals): %d trials ' ...
         'per stage; if run serially, a legacy stage takes about %.0f min and a ' ...
         'refined stage about %.0f min; the large-sample sets run in parfor.\n'], ...
        CFG.t54_reg_B/1e6, n_stage, n_stage * tot(1) / 60, n_stage * tot(2) / 60);

Tb = table(repmat(string(F.run_tag), 2, 1), string(gm), qs(:,1), qs(:,2), ...
           repmat(nb, 2, 1), m(:,1), m(:,2), m(:,3), m(:,4), m(:,5), tot, ...
           repmat(ratio, 2, 1), repmat(CFG.t54_stop_ratio, 2, 1), ...
           repmat(string(loc_tern(go, 'CONTINUE', 'STOP')), 2, 1), ...
           repmat(string(loc_t54_scope_text()), 2, 1), ...
           'VariableNames', {'run_tag', 'grid', 'Q_scan_th', 'Q_scan_u', 'n_trials', ...
           'gen_s', 'bpd_wb_s', 'armA_s', 'armsC_s', 'armsB_s', 'total_s', ...
           'ratio_refined_over_legacy', 'stop_limit', 'stop_rule', 'scope_note'});
writetable(Tb, loc_t54_file(F, 'bench'));
end


function [GT, halt] = t54_gt3(CFG, F, GT)
%T54_GT3  Gate G-T3 (Spec Sec. 4.7) and the T54-2 behaviour check.
%  G-T3: with P.phaseD_anchor_off absent, and with it false, the T-54 copy
%  wb_clkl_driver_pc gives, on every field the harness reports for arms A,
%  C_shared and C_percarrier, the values of the pre-patch copy
%  wb_clkl_driver_pc_prepatch (isequaln; the wall-clock columns excluded).
%  T54-2 check: with the field true, every call of an arm in
%  CFG.t54_anchor_arms is labelled 'a', L_PhaseD is unchanged, and a call
%  that selected candidate A with the anchor on returns the same estimate.
halt  = {};
n     = CFG.t54_gt3_n;
seeds = CFG.seed_base + CFG.t54_gt3_offset + (1:n);
geom  = CFG.t54_reg_leg;
fprintf('\n### T-54 gate G-T3 and T54-2 check (%d trials, seeds %d..%d, leg %s, B = %.0f MHz, serial) ###\n', ...
        n, seeds(1), seeds(end), geom, CFG.t54_reg_B/1e6);
[P, pinfo] = t40_P(CFG, CFG.t54_reg_B);
K  = struct('slab', 'GT3', 'nd', CFG.t54_nd, 'nl', CFG.t54_n_large, ...
            'first', true, 'smoke', false);
GT = loc_t54_Pgates(P, pinfo, CFG, F, K, GT);

C0 = CFG;  C0.t54_anchor_off  = false;    % field absent
C1 = C0;   C1.t54_prepatch    = true;     % pre-patch copy
C2 = C0;   C2.t54_force_false = true;     % field present and false
C3 = C0;   C3.t54_anchor_off  = true;     % field true (T54-2)
[TR0, SC0] = t54_run_set(P, C0, geom, seeds, false);
[TR1, SC1] = t54_run_set(P, C1, geom, seeds, false);
[TR2, SC2] = t54_run_set(P, C2, geom, seeds, false);
[TR3, SC3] = t54_run_set(P, C3, geom, seeds, false);

tn  = t54_cols('trial');
sn  = t54_cols('sub');
cmp = ~ismember(tn, {'rt_gen', 'rt_bpd', 'rt_A', 'rt_C', 'rt_B'});
e01 = false(n, 1);
e21 = false(n, 1);
for i = 1:n
    e01(i) = isequaln(TR0(i, cmp), TR1(i, cmp)) && isequaln(SC0(:, :, i), SC1(:, :, i));
    e21(i) = isequaln(TR2(i, cmp), TR1(i, cmp)) && isequaln(SC2(:, :, i), SC1(:, :, i));
end
v0  = (loc_t54_col(TR0, 'trial', 'valid') == 1);
v1  = (loc_t54_col(TR1, 'trial', 'valid') == 1);
v2  = (loc_t54_col(TR2, 'trial', 'valid') == 1);
v3  = (loc_t54_col(TR3, 'trial', 'valid') == 1);
nv  = sum(v0 & v1 & v2);
ok3 = all(e01) && all(e21) && (nv == n);
fprintf(['  G-T3 [bitwise domain: isequaln on %d per-trial and %d per-subcarrier columns ' ...
         'of arms A, C_shared, C_percarrier; the five wall-clock columns excluded; ' ...
         'grid %s]: field absent equals the pre-patch copy in %d of %d trials; field ' ...
         'false equals the pre-patch copy in %d of %d trials; trials valid in all ' ...
         'three variants %d of %d: %s\n'], ...
        sum(cmp), numel(sn), F.gtag, sum(e01), n, sum(e21), n, nv, n, ...
        loc_tern(ok3, 'PASS', 'FAIL'));
GT = loc_t54_gate(GT, 'G-T3', 'GT3', 'registered cell; 50-trial seed set', ...
                  'bitwise (isequaln)', sum(e01) + sum(e21), 2 * n, ok3, 'yes', '');
if ~ok3; halt{end+1} = 'G-T3 FAIL'; end

la = ismember(tn, {'A_L_A', 'A_L_B', 'A_L_C'});
lc = ismember(sn, {'C_L_A', 'C_L_B', 'C_L_C', 'Cp_L_A', 'Cp_L_B', 'Cp_L_C'});
ja = ismember(tn, {'A_theta_deg', 'A_r_m'});
a0 = (loc_t54_col(TR0, 'trial', 'A_sel') == 1);
c0 = (loc_t54_sub(SC0, 'C_sel') == 1);
r0 = loc_t54_sub(SC0, 'C_r_m');
r3 = loc_t54_sub(SC3, 'C_r_m');
h0 = loc_t54_sub(SC0, 'C_theta_deg');
h3 = loc_t54_sub(SC3, 'C_theta_deg');
sel_ok = (sum(loc_t54_col(TR3, 'trial', 'v_sel_unexpected')) == 0);
L_ok   = isequaln(TR3(:, la), TR0(:, la)) && isequaln(SC3(:, lc, :), SC0(:, lc, :));
keep_A = isequaln(TR3(a0, ja), TR0(a0, ja));
keep_C = isequaln(r3(c0), r0(c0)) && isequaln(h3(c0), h0(c0));
ok2    = sel_ok && L_ok && keep_A && keep_C && all(v3);
fprintf(['  T54-2 CHECK [bitwise and integer domain; field true against field absent]: ' ...
         'labels as the stage expects in every call = %d | L_PhaseD unchanged = %d | ' ...
         'arm A estimate unchanged in the %d of %d trials that selected A = %d | ' ...
         'C_shared estimate unchanged in the %d of %d calls that selected A = %d | ' ...
         'valid trials %d of %d: %s\n'], ...
        sel_ok, L_ok, sum(a0), n, keep_A, sum(c0(:)), numel(c0), keep_C, sum(v3), n, ...
        loc_tern(ok2, 'PASS', 'FAIL'));
GT = loc_t54_gate(GT, 'T54-2-CHECK', 'GT3', 'registered cell; 50-trial seed set', ...
                  'bitwise and integer', double(~ok2), 0, ok2, 'yes', '');
if ~ok2; halt{end+1} = 'T54-2 CHECK FAIL'; end

lab = {'GT3-field-absent', 'GT3-prepatch', 'GT3-field-false', 'GT3-field-true'};
TRs = {TR0, TR1, TR2, TR3};
SCs = {SC0, SC1, SC2, SC3};
Cs  = {C0, C1, C2, C3};
for q = 1:4
    Cq = loc_t54_ctx(Cs{q}, F, 'GT3', P, geom);
    loc_t54_write(loc_t54_trial_table(TRs{q}, Cq, lab{q}), loc_t54_file(F, 'trials'));
    loc_t54_write(loc_t54_sub_table(SCs{q}, seeds, Cq, lab{q}), ...
                  loc_t54_file(F, 'subcarriers'));
end
end


function st = t54_cell(P, CFG, F, C, nd, nl)
%T54_CELL  One (stage, leg, bandwidth) cell of T-54 (Spec Sec. 4.4).
%  S6 set: nd trials, serial, seeds CFG.seed_base + 500 + (1:nd). S6 is
%  mean(degen_C) over this set and passes iff it is below 0.50: definition
%  and tolerance unchanged (FIX 5; F-091). Large-sample set: nl trials, in
%  parfor, seeds CFG.seed_base + CFG.t54_large_offset + (1:nl): a reporting
%  statistic with a Wilson interval, never called S6, never gated.
fprintf('\n--- %s: cell [%s] ---\n', C.slab, C.cell_label);

seeds6 = CFG.seed_base + CFG.t54_s6_offset + (1:nd);
t0 = tic;
[TR6, SC6] = t54_run_set(P, CFG, C.geom, seeds6, false);
s6 = loc_t54_stats(TR6, SC6, C, 'S6', seeds6, toc(t0));
p6 = s6.degen_C_mean < CFG.t54_s6_tol;
fprintf('  %-24s val=%.3e tol=%.1e  %s   [%s | nd=%d seeds %d..%d | k=%d]\n', ...
        'S6_C_not_degenerate', s6.degen_C_mean, CFG.t54_s6_tol, ...
        loc_tern(p6, 'PASS', 'FAIL'), C.cell_label, nd, seeds6(1), seeds6(end), ...
        s6.k_degen_C);
loc_t54_print(s6, C);

seedsL = CFG.seed_base + CFG.t54_large_offset + (1:nl);
t0 = tic;
[TRL, SCL] = t54_run_set(P, CFG, C.geom, seedsL, true);
sL = loc_t54_stats(TRL, SCL, C, 'LARGE', seedsL, toc(t0));
loc_t54_print(sL, C);

loc_t54_write(loc_t54_trial_table(TR6, C, 'S6'),          loc_t54_file(F, 'trials'));
loc_t54_write(loc_t54_trial_table(TRL, C, 'LARGE'),       loc_t54_file(F, 'trials'));
loc_t54_write(loc_t54_sub_table(SC6, seeds6, C, 'S6'),    loc_t54_file(F, 'subcarriers'));
loc_t54_write(loc_t54_sub_table(SCL, seedsL, C, 'LARGE'), loc_t54_file(F, 'subcarriers'));
loc_t54_write(struct2table(s6, 'AsArray', true),          loc_t54_file(F, 'cells'));
loc_t54_write(struct2table(sL, 'AsArray', true),          loc_t54_file(F, 'cells'));

st = struct('n_s6_lines', 1, 'n_sel_lines', 2);
st.s6    = s6;
st.large = sL;
end


function [TR, SC] = t54_run_set(P, CFG, geom, seeds, use_par)
%T54_RUN_SET  Runs t54_one for every seed. TR is n x NT, SC is K_s x NS x n.
%  Every trial seeds its own generator (t40_trial), so the result does not
%  depend on the order of the iterations.
n   = numel(seeds);
NT  = numel(t54_cols('trial'));
NS  = numel(t54_cols('sub'));
K_s = P.K_s;
TR  = nan(n, NT);
SC  = nan(K_s, NS, n);
if use_par
    parfor i = 1:n
        [a, s] = t54_one(P, CFG, geom, seeds(i)); %#ok<PFBNS>
        TR(i, :)    = a;
        SC(:, :, i) = s;
    end
else
    for i = 1:n
        [a, s] = t54_one(P, CFG, geom, seeds(i));
        TR(i, :)    = a;
        SC(:, :, i) = s;
    end
end
end


function [tr, sc] = t54_one(P, CFG, geom, seed)
%T54_ONE  One T-54 trial as numeric rows: tr is 1 x numel(t54_cols('trial')),
%  sc is K_s x numel(t54_cols('sub')). A trial is VALID when every v_* flag
%  is zero: no absorbed error, no fallback anchor, no non-finite estimate or
%  objective, and every Phase D label is one the stage expects.
R    = t40_trial(P, CFG, CFG.SNR_star, seed, geom, false);
T    = R.t54;
K_s  = P.K_s;
doCp = double(CFG.do_C_percarrier);
if ~isempty(T.msg)
    fprintf('  T54 ABSORBED ERROR, seed %d (first of the trial): %s\n', seed, T.msg);
end
offA  = CFG.t54_anchor_off && any(strcmp('A',            CFG.t54_anchor_arms));
offC  = CFG.t54_anchor_off && any(strcmp('C_shared',     CFG.t54_anchor_arms));
offCp = CFG.t54_anchor_off && any(strcmp('C_percarrier', CFG.t54_anchor_arms));

v = struct();
v.seed             = seed;
v.theta_true_deg   = R.theta_true * 180/pi;
v.r_true_m         = R.r_true;
v.bpd_theta_deg    = T.th_bpd * 180/pi;
v.bpd_r_m          = T.r_bpd;
v.A_theta_deg      = T.thA * 180/pi;
v.A_r_m            = T.rA;
v.A_sel            = T.selA;
v.A_L_A            = T.LdA(1);
v.A_L_B            = T.LdA(2);
v.A_L_C            = T.LdA(3);
v.A_n_iter         = R.n_iter.A;
v.A_converged      = R.converged.A;
v.degen_C          = R.degen_C;
v.degen_Cp         = R.degen_Cp;
v.degen_B          = R.degen_B;
v.disp_omega       = R.disp_omega;
v.disp_kappa       = R.disp_kappa;
v.disp_omega_Cp    = std(T.omCp, 'omitnan');
v.disp_kappa_Cp    = std(T.kaCp, 'omitnan');
v.n_u_C            = numel(unique(T.rC(isfinite(T.rC))));
v.n_u_Cp           = numel(unique(T.rCp(isfinite(T.rCp))));
v.v_bpd_wb_fail    = T.bpd_fail;
v.v_A_exc          = T.exA;
v.v_A_nonfinite    = double(~(isfinite(T.thA) && isfinite(T.rA)));
v.v_C_exc          = sum(T.exC);
v.v_Cp_bpd_fail    = sum(T.bpdk_fail);
v.v_Cp_exc         = sum(T.exCp);
v.v_C_nonfinite    = sum(~isfinite(T.omC) | ~isfinite(T.kaC));
v.v_Cp_nonfinite   = doCp * sum(~isfinite(T.omCp) | ~isfinite(T.kaCp));
v.v_L_nonfinite    = double(any(~isfinite(T.LdA))) + sum(any(~isfinite(T.LdC), 2)) ...
                     + doCp * sum(any(~isfinite(T.LdCp), 2));
v.v_sel_unexpected = double(~loc_t54_selok(T.selA, offA)) ...
                     + sum(~loc_t54_selok(T.selC, offC)) ...
                     + doCp * sum(~loc_t54_selok(T.selCp, offCp));
v.valid            = double((v.v_bpd_wb_fail + v.v_A_exc + v.v_A_nonfinite ...
                     + v.v_C_exc + v.v_Cp_bpd_fail + v.v_Cp_exc + v.v_C_nonfinite ...
                     + v.v_Cp_nonfinite + v.v_L_nonfinite + v.v_sel_unexpected) == 0);
v.rt_gen           = R.rt_gen;
v.rt_bpd           = R.rt_bpd;
v.rt_A             = R.rt_A;
v.rt_C             = R.rt_C;
v.rt_B             = R.rt_B;

w = struct();
w.k                = (1:K_s).';
w.alpha_k          = P.alpha_k_vec(:);
w.C_theta_deg      = T.thC * 180/pi;
w.C_r_m            = T.rC;
w.C_omega          = T.omC;
w.C_kappa          = T.kaC;
w.C_sel            = T.selC;
w.C_L_A            = T.LdC(:, 1);
w.C_L_B            = T.LdC(:, 2);
w.C_L_C            = T.LdC(:, 3);
w.C_n_iter         = T.niC;
w.C_converged      = T.cvC;
w.C_exc            = T.exC;
w.Cp_bpd_theta_deg = T.tbv * 180/pi;
w.Cp_bpd_r_m       = T.rbv;
w.Cp_bpd_fail      = T.bpdk_fail;
w.Cp_theta_deg     = T.thCp * 180/pi;
w.Cp_r_m           = T.rCp;
w.Cp_omega         = T.omCp;
w.Cp_kappa         = T.kaCp;
w.Cp_sel           = T.selCp;
w.Cp_L_A           = T.LdCp(:, 1);
w.Cp_L_B           = T.LdCp(:, 2);
w.Cp_L_C           = T.LdCp(:, 3);
w.Cp_n_iter        = T.niCp;
w.Cp_converged     = T.cvCp;
w.Cp_exc           = T.exCp;

tn = t54_cols('trial');
sn = t54_cols('sub');
assert(numel(fieldnames(v)) == numel(tn) && numel(fieldnames(w)) == numel(sn), ...
       't54_one: the row structs and t54_cols disagree.');
tr = nan(1, numel(tn));
for j = 1:numel(tn)
    tr(j) = v.(tn{j});
end
sc = nan(K_s, numel(sn));
for j = 1:numel(sn)
    sc(:, j) = w.(sn{j});
end
end


function names = t54_cols(kind)
%T54_COLS  Column names of the T-54 numeric rows. 'trial': one row per
%  trial. 'sub': one row per trial and subcarrier. A *_sel column holds a
%  Phase D selection as a code: 0 '?', 1 'A', 2 'B', 3 'C', 4 'a'; the CSV
%  writers replace it by a *_PhaseD_select column of characters.
switch kind
    case 'trial'
        names = {'seed', 'theta_true_deg', 'r_true_m', 'bpd_theta_deg', 'bpd_r_m', ...
                 'A_theta_deg', 'A_r_m', 'A_sel', 'A_L_A', 'A_L_B', 'A_L_C', ...
                 'A_n_iter', 'A_converged', 'degen_C', 'degen_Cp', 'degen_B', ...
                 'disp_omega', 'disp_kappa', 'disp_omega_Cp', 'disp_kappa_Cp', ...
                 'n_u_C', 'n_u_Cp', 'v_bpd_wb_fail', 'v_A_exc', 'v_A_nonfinite', ...
                 'v_C_exc', 'v_Cp_bpd_fail', 'v_Cp_exc', 'v_C_nonfinite', ...
                 'v_Cp_nonfinite', 'v_L_nonfinite', 'v_sel_unexpected', 'valid', ...
                 'rt_gen', 'rt_bpd', 'rt_A', 'rt_C', 'rt_B'};
    case 'sub'
        names = {'k', 'alpha_k', 'C_theta_deg', 'C_r_m', 'C_omega', 'C_kappa', ...
                 'C_sel', 'C_L_A', 'C_L_B', 'C_L_C', 'C_n_iter', 'C_converged', ...
                 'C_exc', 'Cp_bpd_theta_deg', 'Cp_bpd_r_m', 'Cp_bpd_fail', ...
                 'Cp_theta_deg', 'Cp_r_m', 'Cp_omega', 'Cp_kappa', 'Cp_sel', ...
                 'Cp_L_A', 'Cp_L_B', 'Cp_L_C', 'Cp_n_iter', 'Cp_converged', 'Cp_exc'};
    otherwise
        error('t40_e2b_driver:T54Cols', 't54_cols: unknown kind %s.', kind);
end
end


% -------------------------------------------------------------------------
%  T-54 local helpers
% -------------------------------------------------------------------------
function tag = loc_t54_runtag(out_dir)
%LOC_T54_RUNTAG  The tag of a T-54 run: exactly one file named
%  t54_runtag_<tag>.txt must exist in out_dir. The run sheet creates it
%  before stage P0; the stages never create it.
d = dir(fullfile(out_dir, 't54_runtag_*.txt'));
if numel(d) ~= 1
    error('t40_e2b_driver:T54RunTag', ...
          ['T-54 stages need exactly one file t54_runtag_<tag>.txt in %s; found %d. ' ...
           'The run sheet creates it before stage P0.'], out_dir, numel(d));
end
nm  = d(1).name;
tag = nm(numel('t54_runtag_') + 1 : end - numel('.txt'));
if isempty(tag) || ~all(isstrprop(tag, 'alphanum') | tag == '_')
    error('t40_e2b_driver:T54RunTag', ...
          'The T-54 run tag ''%s'' is empty or holds a character other than a letter, a digit or an underscore.', tag);
end
end

function [GT, halt] = loc_t54_preamble(F, grid_mode)
%LOC_T54_PREAMBLE  Top of every T-54 console: the scope sentence, the run
%  tag, and the gates that need no trial -- the path check, PC-T3 and G-T1.
%  They read the files the path resolves NOW, i.e. the code this stage runs.
GT   = cell(0, 9);
halt = {};
loc_t54_scope();
fprintf('  Governing record: PaperC_T40_Gate_Amendment_2026-09-03.md; its disclosure record travels with any write-up of T-54.\n');
fprintf('  [T-54] run tag %s | stage %s | grid argument ''%s'' | console stamp %s\n', ...
        F.run_tag, F.stage, grid_mode, F.stamp);
fprintf('  [T-54] MATLAB version %s\n', version);

fn = {'t40_e2b_driver', 'wb_clkl_driver_pc', 'wb_clkl_driver_pc_prepatch', ...
      'nf_clkl_pc', 'wb_clkl_estimator', 'setup_production_P_v4', 'bpd_baseline', ...
      'wb_channel_gen_ofdm_nf', 'wb_nf_fresnel_steer', 'nf_usw_steer', ...
      'run_monte_carlo_paperC'};
nbad = 0;
for i = 1:numel(fn)
    wh = which(fn{i}, '-all');
    if ischar(wh)
        if isempty(wh); wh = {}; else; wh = {wh}; end
    end
    first = '(not found)';
    if ~isempty(wh); first = wh{1}; end
    fprintf('  PATH %-28s hits=%d  %s\n', fn{i}, numel(wh), first);
    nbad = nbad + double(numel(wh) ~= 1);
end
ok = (nbad == 0);
fprintf('  PATH-CHECK [integer domain; exactly one file on the path per function, %d functions]: functions with a hit count other than one = %d: %s\n', ...
        numel(fn), nbad, loc_tern(ok, 'PASS', 'FAIL'));
GT = loc_t54_gate(GT, 'PATH-CHECK', F.stage, 'every function the harness reaches', ...
                  'integer', nbad, 0, ok, 'yes', '');
if ~ok; halt{end+1} = 'PATH-CHECK FAIL'; end

n1 = loc_t54_count('wb_clkl_driver_pc', ['info.L_PhaseD      = ' '[L_A, L_B, L_C];']);
n2 = loc_t54_count('wb_clkl_driver_pc', ['info.PhaseD_select = ' 'PhaseD_labels(sel_idx);']);
ok = (n1 == 1) && (n2 == 1);
fprintf(['  PC-T3 [integer domain; text of the wb_clkl_driver_pc file on the path]: lines ' ...
         'assigning info.L_PhaseD = [L_A, L_B, L_C]: %d; lines assigning ' ...
         'info.PhaseD_select = PhaseD_labels(sel_idx): %d (each must be 1): %s\n'], ...
        n1, n2, loc_tern(ok, 'PASS', 'FAIL'));
GT = loc_t54_gate(GT, 'PC-T3', F.stage, 'wb_clkl_driver_pc on the path', 'integer', ...
                  n1 + n2, 2, ok, 'yes', '');
if ~ok; halt{end+1} = 'PC-T3 FAIL'; end

pat_h = ['A_k(:,l) = wb_nf_fresnel_steer(theta_true, 1/r_true, ' 'ak, P) * sqrt(M);'];
pat_p = ['A_k(:,l) = wb_nf_fresnel_steer(theta_true(l), 1/r_true(l), ' 'alpha_k, P) * sqrt(M);'];
n1 = loc_t54_count('t40_e2b_driver', pat_h);
n2 = loc_t54_count('run_monte_carlo_paperC', pat_p);
ok = (n1 == 1) && (n2 == 1);
fprintf(['  G-T1 [integer domain; text of the files on the path]: F-088 patched steering ' ...
         'line in loc_gen_fixed of t40_e2b_driver: %d; in ' ...
         'wb_channel_gen_ofdm_nf_fixed of run_monte_carlo_paperC: %d (each must be ' ...
         '1): %s\n'], n1, n2, loc_tern(ok, 'PASS', 'FAIL'));
GT = loc_t54_gate(GT, 'G-T1', F.stage, 't40_e2b_driver and run_monte_carlo_paperC on the path', ...
                  'integer', n1 + n2, 2, ok, 'yes', '');
if ~ok; halt{end+1} = 'G-T1 FAIL'; end
end

function GT = loc_t54_Pgates(P, pinfo, CFG, F, K, GT)
%LOC_T54_PGATES  Gates that need the parameter struct. PC-T4 (halting) at
%  every t40_P call of a T-54 stage; G-T4 (non-halting) and the record of
%  every defaulted argument (L-34) once per stage (K.first).
if CFG.use_refined; qreg = [768 2048]; else; qreg = [192 256]; end
ok = (CFG.seed_base == 2000000) && isequal([P.Q_scan_th, P.Q_scan_u], qreg) ...
     && (CFG.pin_r == 2.13) && (CFG.pin_r_far == 3.90) ...
     && (CFG.pin_theta == 40 * pi/180) && (CFG.SNR_star == 5) ...
     && isequal(CFG.B_grid, [100 200 400 600] * 1e6) && CFG.do_C_percarrier ...
     && (~CFG.t54_dflt_halt || (~pinfo.setup_rejected && pinfo.n_filled == 0)) ...
     && (K.smoke || (K.nd == 12 && K.nl == 200));
fprintf(['  PC-T4 [%s | B=%.0f MHz K_s=%d; exact-value domain]: seed_base=%d | ' ...
         'Q_scan=[%d, %d] (grid %s) | nd=%d | n_large=%d | S6 seeds %d..%d | ' ...
         'large-sample seeds %d..%d | SNR=%g dB | legs: pinned r=%.2f m theta=%.2f deg ' ...
         '(gate geometry); pinned_far r=%.2f m theta=%.2f deg (added post hoc, ' ...
         'reported never gated) | legs run: %s | Phase D anchor off=%d, arms {%s} | ' ...
         'fields filled from locked defaults=%d | setup call rejected=%d: %s\n'], ...
        K.slab, P.B/1e6, P.K_s, CFG.seed_base, P.Q_scan_th, P.Q_scan_u, ...
        loc_tern(CFG.use_refined, 'refined', 'legacy'), K.nd, K.nl, ...
        CFG.seed_base + CFG.t54_s6_offset + 1, CFG.seed_base + CFG.t54_s6_offset + K.nd, ...
        CFG.seed_base + CFG.t54_large_offset + 1, CFG.seed_base + CFG.t54_large_offset + K.nl, ...
        CFG.SNR_star, CFG.pin_r, CFG.pin_theta * 180/pi, CFG.pin_r_far, ...
        CFG.pin_theta * 180/pi, strjoin(CFG.t54_legs, ', '), ...
        double(CFG.t54_anchor_off), strjoin(CFG.t54_anchor_arms, ', '), ...
        pinfo.n_filled, double(pinfo.setup_rejected), loc_tern(ok, 'PASS', 'FAIL'));
GT = loc_t54_gate(GT, 'PC-T4', K.slab, sprintf('B=%.0f MHz', P.B/1e6), 'exact value', ...
                  double(~ok), 0, ok, 'yes', pinfo.filled);
if ~ok
    loc_t54_finish(GT, {[K.slab ': PC-T4 FAIL']}, F);
end
if ~K.first; return; end

du  = (P.u_max - P.u_min) / (CFG.Q_scan_u_ref - 1);
rel = (du - CFG.t54_du_ref) / CFG.t54_du_ref;
ok4 = abs(rel) <= CFG.t54_du_tol;
dr  = du * CFG.pin_r^2;
fprintf(['  G-T4 [relative domain; non-halting]: du_refined = (u_max - u_min) / ' ...
         '(Q_scan_u_ref - 1) = %.12e 1/m | reference %.12e 1/m | relative deviation ' ...
         '%+.3e | tolerance %.1e: %s\n'], ...
        du, CFG.t54_du_ref, rel, CFG.t54_du_tol, loc_tern(ok4, 'PASS', 'FAIL'));
fprintf(['  G-T4 arithmetic (reporting; it sizes the Phase C u-scan step only and is not ' ...
         'a floor for the reported range: F-110, F-124): Q_scan_u_ref = %d | ' ...
         'dr_refined at r = %.2f m = du * r^2 = %.9e m | dr / sqrt(12) = %.9e m | ' ...
         'step ratio legacy to refined = %.4f\n'], ...
        CFG.Q_scan_u_ref, CFG.pin_r, dr, dr / sqrt(12), (CFG.Q_scan_u_ref - 1) / 255);
fprintf(['  G-T4 cross-check (reported; not a registered gate): the locked box and ' ...
         'Q_scan_u = 2048 imply du = %.15e 1/m; relative deviation of du_refined ' ...
         'from it %+.3e\n'], CFG.t54_du_box, (du - CFG.t54_du_box) / CFG.t54_du_box);
GT = loc_t54_gate(GT, 'G-T4', K.slab, 'refined u-scan step', 'relative', ...
                  du, CFG.t54_du_ref, ok4, 'no', sprintf('relative deviation %+.3e', rel));

fprintf(['  EFFECTIVE DEFAULTS (L-34): wb_clkl_driver_pc n_ms_starts=%s, ' ...
         'short_iter_ms=%s, use_multi_start=%s, ablation_skip_scan=%s, ' ...
         'ablation_joint_update=%s | bpd_baseline bpd_do_ls_polish=%s, ' ...
         'bpd_write_csv=%s | nf_clkl_pc do_post_loop_scan=true (argument not ' ...
         'passed)\n'], ...
        loc_t54_eff(P, 'n_ms_starts', '4'), loc_t54_eff(P, 'short_iter_ms', '20'), ...
        loc_t54_eff(P, 'use_multi_start', 'true'), ...
        loc_t54_eff(P, 'ablation_skip_scan', 'false'), ...
        loc_t54_eff(P, 'ablation_joint_update', 'true'), ...
        loc_t54_eff(P, 'bpd_do_ls_polish', 'true'), loc_t54_eff(P, 'bpd_write_csv', 'true'));
fprintf(['  DATA MODEL (L-34): pinned legs are generated by loc_gen_fixed with Fresnel ' ...
         'steering (wb_nf_fresnel_steer), the estimator model; it calls ' ...
         'wb_channel_gen_ofdm_nf for N0 and the combiner only (its use_exact ' ...
         'argument defaults to true and its snapshots are discarded).\n']);
end

function loc_t54_finish(GT, halt, F)
%LOC_T54_FINISH  Writes the gates CSV of the stage, then stops the stage if
%  a halting condition fired.
if ~isempty(GT)
    Tg = cell2table(GT, 'VariableNames', {'gate', 'stage_label', 'scope', 'domain', ...
         'value', 'reference', 'verdict', 'halting', 'note'});
    Tg.run_tag    = repmat(string(F.run_tag), size(GT, 1), 1);
    Tg.scope_note = repmat(string(loc_t54_scope_text()), size(GT, 1), 1);
    writetable(Tg, loc_t54_file(F, 'gates'));
end
if ~isempty(halt)
    fprintf('\n  T-54 HALT: %s. Run no further stage; capture what exists as a partial run.\n', ...
            strjoin(halt, ' | '));
    error('t40_e2b_driver:T54Halt', 'T-54 HALT in stage %s: %s', F.stage, strjoin(halt, ' | '));
end
end

function GT = loc_t54_gate(GT, gate, slab, scope, domain, value, ref, ok, halting, note)
%LOC_T54_GATE  One row of the gates table. ok is a logical (PASS or FAIL) or
%  the verdict text itself.
if ischar(ok)
    verdict = ok;
else
    verdict = loc_tern(ok, 'PASS', 'FAIL');
end
GT(end+1, :) = {gate, slab, scope, domain, double(value), double(ref), verdict, halting, note};
end

function n = loc_t54_count(fname, pat)
%LOC_T54_COUNT  Occurrences of the text pat in the file fname resolves to on
%  the path; -1 if fname is not found.
n  = -1;
wh = which(fname);
if isempty(wh); return; end
n = numel(strfind(fileread(wh), pat));
end

function t = loc_t54_eff(P, f, dflt)
%LOC_T54_EFF  Effective value of an optional P field, for the L-34 record.
if isfield(P, f)
    t = [mat2str(P.(f)) ' (set in P)'];
else
    t = [dflt ' (default; field absent from P)'];
end
end

function t = loc_t54_scope_text()
t = 'T-54 is diagnostic only; no T-40 or T-54 number is a result; no arm-B number is usable (F-092).';
end

function loc_t54_scope()
fprintf('  %s\n', loc_t54_scope_text());
end

function slab = loc_t54_slab(use_refined, anchor_off, smoke)
%LOC_T54_SLAB  Stage label of a self-test call (Spec Sec. 4.4): S1 legacy and
%  S2 refined with the anchor on; S3 refined and S4 legacy with it off.
if smoke
    slab = loc_tern(anchor_off, 'S0-SMOKE-NOANCHOR', 'S0-SMOKE-ANCHOR');
elseif anchor_off
    slab = loc_tern(use_refined, 'S3', 'S4');
else
    slab = loc_tern(use_refined, 'S2', 'S1');
end
end

function C = loc_t54_ctx(CFG, F, slab, P, geom)
%LOC_T54_CTX  Context of one cell: what its console lines and CSV rows carry.
C = struct();
C.run_tag    = F.run_tag;
C.slab       = slab;
C.gtag       = loc_tern(CFG.use_refined, 'refined', 'legacy');
C.Q_scan_th  = P.Q_scan_th;
C.Q_scan_u   = P.Q_scan_u;
C.anchor_off = CFG.t54_anchor_off;
C.geom       = geom;
switch geom
    case 'pinned'
        C.leg_label     = sprintf('pinned r=%.2f m theta=%.1f deg (gate geometry)', ...
                                  CFG.pin_r, CFG.pin_theta * 180/pi);
        C.leg_label_csv = 'gate geometry';
    case 'pinned_far'
        C.leg_label     = sprintf('pinned_far r=%.2f m theta=%.1f deg (added post hoc, reported never gated)', ...
                                  CFG.pin_r_far, CFG.pin_theta * 180/pi);
        C.leg_label_csv = 'added post hoc - reported never gated';
    otherwise
        C.leg_label     = 'randomised scene box (reported never gated; exact-USW data, use_exact defaulted true)';
        C.leg_label_csv = 'randomised - reported never gated';
end
C.is_reg     = strcmp(geom, CFG.t54_reg_leg) && (P.B == CFG.t54_reg_B);
C.B_hz       = P.B;
C.K_s        = P.K_s;
C.SNR_dB     = CFG.SNR_star;
C.cell_label = sprintf('%s | grid %s | anchor %s | %s | B=%.0f MHz K_s=%d | SNR=%g dB%s', ...
                       slab, C.gtag, loc_tern(C.anchor_off, 'OFF', 'ON'), C.leg_label, ...
                       P.B/1e6, P.K_s, CFG.SNR_star, ...
                       loc_tern(C.is_reg, ' | REGISTERED CELL', ''));
end

function Pa = loc_t54_P(P, CFG, arm)
%LOC_T54_P  The parameter struct one wb_clkl_driver_pc call of arm 'A',
%  'C_shared' or 'C_percarrier' receives (change T54-2). Outside stage
%  P2_SELFTEST_NOANCHOR and gate G-T3 it returns P unchanged.
Pa = P;
if CFG.t54_anchor_off && any(strcmp(arm, CFG.t54_anchor_arms))
    Pa.phaseD_anchor_off = true;
elseif CFG.t54_force_false
    Pa.phaseD_anchor_off = false;
end
end

function c = loc_t54_selcode(ch)
%LOC_T54_SELCODE  Phase D selection as a code: 'A' 1, 'B' 2, 'C' 3, 'a' 4;
%  anything else, '?' included, 0.
c = 0;
if ischar(ch) && numel(ch) == 1
    k = find('ABCa' == ch, 1);
    if ~isempty(k); c = k; end
end
end

function ok = loc_t54_selok(code, off)
%LOC_T54_SELOK  True where a Phase D code is one the stage expects: 'a' (4)
%  for an arm whose anchor is disabled, otherwise 'A', 'B' or 'C' (1 to 3).
if off
    ok = (code == 4);
else
    ok = (code >= 1) & (code <= 3);
end
end

function m = loc_t54_first(m, where, emsg)
%LOC_T54_FIRST  Keeps the first absorbed error message of a trial.
if isempty(m) && ~isempty(emsg)
    m = [where ': ' emsg];
end
end

function s = loc_t54_stats(TR, SC, C, setname, seeds, elapsed_s)
%LOC_T54_STATS  Statistics of one trial set of one cell (Spec Sec. 4.4).
%  Every count is over ALL trials of the set: no trial is dropped. The
%  *_valid fields repeat the degeneracy counts on the trials whose validity
%  flags are all zero. For the S6 set, degen_C_mean is S6. Means and SDs of
%  L_PhaseD are over the finite values (SD with the N-1 normalisation).
n    = size(TR, 1);
dC   = loc_t54_col(TR, 'trial', 'degen_C');
dCp  = loc_t54_col(TR, 'trial', 'degen_Cp');
val  = (loc_t54_col(TR, 'trial', 'valid') == 1);
cA   = loc_t54_col(TR, 'trial', 'A_sel');
cC   = loc_t54_sub(SC, 'C_sel');
cCp  = loc_t54_sub(SC, 'Cp_sel');
bcC  = (cC == 2) | (cC == 3);
bcCp = (cCp == 2) | (cCp == 3);
LA   = [loc_t54_col(TR, 'trial', 'A_L_A'), loc_t54_col(TR, 'trial', 'A_L_B'), ...
        loc_t54_col(TR, 'trial', 'A_L_C')];
LC   = [reshape(loc_t54_sub(SC, 'C_L_A'), [], 1), reshape(loc_t54_sub(SC, 'C_L_B'), [], 1), ...
        reshape(loc_t54_sub(SC, 'C_L_C'), [], 1)];
LP   = [reshape(loc_t54_sub(SC, 'Cp_L_A'), [], 1), reshape(loc_t54_sub(SC, 'Cp_L_B'), [], 1), ...
        reshape(loc_t54_sub(SC, 'Cp_L_C'), [], 1)];
nu   = loc_t54_col(TR, 'trial', 'n_u_C');
nup  = loc_t54_col(TR, 'trial', 'n_u_Cp');
dO   = loc_t54_col(TR, 'trial', 'disp_omega');
dK   = loc_t54_col(TR, 'trial', 'disp_kappa');
rt   = loc_t54_col(TR, 'trial', 'rt_gen') + loc_t54_col(TR, 'trial', 'rt_bpd') ...
       + loc_t54_col(TR, 'trial', 'rt_A') + loc_t54_col(TR, 'trial', 'rt_C') ...
       + loc_t54_col(TR, 'trial', 'rt_B');

s = struct();
s.run_tag         = string(C.run_tag);
s.stage_label     = string(C.slab);
s.grid            = string(C.gtag);
s.Q_scan_th       = C.Q_scan_th;
s.Q_scan_u        = C.Q_scan_u;
s.anchor_off      = double(C.anchor_off);
s.leg             = string(C.geom);
s.leg_label       = string(C.leg_label_csv);
s.registered_cell = double(C.is_reg);
s.B_hz            = C.B_hz;
s.K_s             = C.K_s;
s.SNR_dB          = C.SNR_dB;
s.trial_set       = string(setname);
s.n               = n;
s.seed_first      = seeds(1);
s.seed_last       = seeds(end);

s.k_degen_C       = sum(dC == 1);
s.degen_C_mean    = mean(dC);
[lo, hi] = loc_wilson(s.k_degen_C, n);
s.degen_C_wilson_lo = lo;
s.degen_C_wilson_hi = hi;
s.k_degen_Cp      = sum(dCp == 1);
s.degen_Cp_mean   = mean(dCp);
[lo, hi] = loc_wilson(s.k_degen_Cp, sum(isfinite(dCp)));
s.degen_Cp_wilson_lo = lo;
s.degen_Cp_wilson_hi = hi;
s.n_valid         = sum(val);
s.n_invalid       = n - sum(val);
s.k_degen_C_valid  = sum(dC(val) == 1);
s.k_degen_Cp_valid = sum(dCp(val) == 1);

s.selA_A   = mean(cA == 1);
s.selA_B   = mean(cA == 2);
s.selA_C   = mean(cA == 3);
s.selA_a   = mean(cA == 4);
s.selA_unk = mean(cA == 0);
s.selA_BC_rate_trials = mean((cA == 2) | (cA == 3));
s.n_calls  = numel(cC);
s.selC_A   = mean(cC(:) == 1);
s.selC_B   = mean(cC(:) == 2);
s.selC_C   = mean(cC(:) == 3);
s.selC_a   = mean(cC(:) == 4);
s.selC_unk = mean(cC(:) == 0);
s.selC_BC_rate_calls     = mean(bcC(:));
s.selC_allBC_rate_trials = mean(all(bcC, 1));
s.selC_allC_rate_trials  = mean(all(cC == 3, 1));
s.selCp_A   = mean(cCp(:) == 1);
s.selCp_B   = mean(cCp(:) == 2);
s.selCp_C   = mean(cCp(:) == 3);
s.selCp_a   = mean(cCp(:) == 4);
s.selCp_unk = mean(cCp(:) == 0);
s.selCp_BC_rate_calls     = mean(bcCp(:));
s.selCp_allBC_rate_trials = mean(all(bcCp, 1));

[mu, sd] = loc_t54_msd(LA(:, 1));  s.LD_A_mean_A = mu;   s.LD_A_sd_A = sd;
[mu, sd] = loc_t54_msd(LA(:, 2));  s.LD_A_mean_B = mu;   s.LD_A_sd_B = sd;
[mu, sd] = loc_t54_msd(LA(:, 3));  s.LD_A_mean_C = mu;   s.LD_A_sd_C = sd;
[mu, sd] = loc_t54_msd(LC(:, 1));  s.LD_C_mean_A = mu;   s.LD_C_sd_A = sd;
[mu, sd] = loc_t54_msd(LC(:, 2));  s.LD_C_mean_B = mu;   s.LD_C_sd_B = sd;
[mu, sd] = loc_t54_msd(LC(:, 3));  s.LD_C_mean_C = mu;   s.LD_C_sd_C = sd;
[mu, sd] = loc_t54_msd(LP(:, 1));  s.LD_Cp_mean_A = mu;  s.LD_Cp_sd_A = sd;
[mu, sd] = loc_t54_msd(LP(:, 2));  s.LD_Cp_mean_B = mu;  s.LD_Cp_sd_B = sd;
[mu, sd] = loc_t54_msd(LP(:, 3));  s.LD_Cp_mean_C = mu;  s.LD_Cp_sd_C = sd;
s.fracLC_A_trials = loc_t54_fracC(LA);
s.fracLC_C_calls  = loc_t54_fracC(LC);
s.fracLC_Cp_calls = loc_t54_fracC(LP);

s.n_u_C_mean  = mean(nu);
s.n_u_C_min   = min(nu);
s.n_u_C_max   = max(nu);
s.n_u_Cp_mean = mean(nup);
s.n_u_Cp_min  = min(nup);
s.n_u_Cp_max  = max(nup);
s.disp_omega_mean   = mean(dO, 'omitnan');
s.disp_omega_median = median(dO, 'omitnan');
s.disp_omega_max    = max(dO);
s.disp_kappa_mean   = mean(dK, 'omitnan');
s.disp_kappa_median = median(dK, 'omitnan');
s.disp_kappa_max    = max(dK);
s.disp_omega_Cp_mean = mean(loc_t54_col(TR, 'trial', 'disp_omega_Cp'), 'omitnan');
s.disp_kappa_Cp_mean = mean(loc_t54_col(TR, 'trial', 'disp_kappa_Cp'), 'omitnan');

s.sum_v_bpd_wb_fail    = sum(loc_t54_col(TR, 'trial', 'v_bpd_wb_fail'));
s.sum_v_A_exc          = sum(loc_t54_col(TR, 'trial', 'v_A_exc'));
s.sum_v_A_nonfinite    = sum(loc_t54_col(TR, 'trial', 'v_A_nonfinite'));
s.sum_v_C_exc          = sum(loc_t54_col(TR, 'trial', 'v_C_exc'));
s.sum_v_Cp_bpd_fail    = sum(loc_t54_col(TR, 'trial', 'v_Cp_bpd_fail'));
s.sum_v_Cp_exc         = sum(loc_t54_col(TR, 'trial', 'v_Cp_exc'));
s.sum_v_C_nonfinite    = sum(loc_t54_col(TR, 'trial', 'v_C_nonfinite'));
s.sum_v_Cp_nonfinite   = sum(loc_t54_col(TR, 'trial', 'v_Cp_nonfinite'));
s.sum_v_L_nonfinite    = sum(loc_t54_col(TR, 'trial', 'v_L_nonfinite'));
s.sum_v_sel_unexpected = sum(loc_t54_col(TR, 'trial', 'v_sel_unexpected'));

s.rt_trial_mean_s = mean(rt);
s.elapsed_s       = elapsed_s;
s.scope_note      = string(loc_t54_scope_text());
end

function loc_t54_print(s, C)
%LOC_T54_PRINT  Console block of one trial set of one cell (Spec Sec. 4.4).
if strcmp(s.trial_set, 'LARGE')
    fprintf('  LARGE-SAMPLE ESTIMATE [%s | n=%d seeds %d..%d] -- reporting statistic, never called S6, never gated:\n', ...
            C.cell_label, s.n, s.seed_first, s.seed_last);
    fprintf('    degen_C  = %d/%d = %.4f, Wilson 95%% [%.4f, %.4f]; among the %d valid trials: %d\n', ...
            s.k_degen_C, s.n, s.degen_C_mean, s.degen_C_wilson_lo, s.degen_C_wilson_hi, ...
            s.n_valid, s.k_degen_C_valid);
else
    fprintf('  TRIAL SET %s [%s | n=%d seeds %d..%d]:\n', ...
            s.trial_set, C.cell_label, s.n, s.seed_first, s.seed_last);
    fprintf('    degen_C  = %d/%d = %.4f; among the %d valid trials: %d\n', ...
            s.k_degen_C, s.n, s.degen_C_mean, s.n_valid, s.k_degen_C_valid);
end
fprintf('    degen_Cp = %d/%d = %.4f, Wilson 95%% [%.4f, %.4f] (C_percarrier control, Add. A.2: REPORTED NEVER GATED; not S6)\n', ...
        s.k_degen_Cp, s.n, s.degen_Cp_mean, s.degen_Cp_wilson_lo, s.degen_Cp_wilson_hi);
fprintf('    PHASE-D SELECTION arm A (per trial, n=%d): A=%.4f B=%.4f C=%.4f a=%.4f ?=%.4f | B+C=%.4f\n', ...
        s.n, s.selA_A, s.selA_B, s.selA_C, s.selA_a, s.selA_unk, s.selA_BC_rate_trials);
fprintf('    PHASE-D SELECTION arm C_shared (per call, n=%d): A=%.4f B=%.4f C=%.4f a=%.4f ?=%.4f | B+C=%.4f | trials with all K_s in {B,C}: %.4f | trials with all K_s = C: %.4f\n', ...
        s.n_calls, s.selC_A, s.selC_B, s.selC_C, s.selC_a, s.selC_unk, ...
        s.selC_BC_rate_calls, s.selC_allBC_rate_trials, s.selC_allC_rate_trials);
fprintf('    PHASE-D SELECTION arm C_percarrier (per call, n=%d): A=%.4f B=%.4f C=%.4f a=%.4f ?=%.4f | B+C=%.4f | trials with all K_s in {B,C}: %.4f\n', ...
        s.n_calls, s.selCp_A, s.selCp_B, s.selCp_C, s.selCp_a, s.selCp_unk, ...
        s.selCp_BC_rate_calls, s.selCp_allBC_rate_trials);
fprintf('    L_PhaseD mean (SD) of candidates A / B / C, arm A: %.6g (%.3g) / %.6g (%.3g) / %.6g (%.3g)\n', ...
        s.LD_A_mean_A, s.LD_A_sd_A, s.LD_A_mean_B, s.LD_A_sd_B, s.LD_A_mean_C, s.LD_A_sd_C);
fprintf('      arm C_shared: %.6g (%.3g) / %.6g (%.3g) / %.6g (%.3g) | arm C_percarrier: %.6g (%.3g) / %.6g (%.3g) / %.6g (%.3g)\n', ...
        s.LD_C_mean_A, s.LD_C_sd_A, s.LD_C_mean_B, s.LD_C_sd_B, s.LD_C_mean_C, s.LD_C_sd_C, ...
        s.LD_Cp_mean_A, s.LD_Cp_sd_A, s.LD_Cp_mean_B, s.LD_Cp_sd_B, s.LD_Cp_mean_C, s.LD_Cp_sd_C);
fprintf('    fraction with L_C <= min(L_A, L_B): arm A (trials) %.4f | C_shared (calls) %.4f | C_percarrier (calls) %.4f\n', ...
        s.fracLC_A_trials, s.fracLC_C_calls, s.fracLC_Cp_calls);
fprintf('    distinct u nodes among the K_s per-subcarrier estimates: C_shared mean %.2f, min %d, max %d | C_percarrier mean %.2f, min %d, max %d\n', ...
        s.n_u_C_mean, s.n_u_C_min, s.n_u_C_max, s.n_u_Cp_mean, s.n_u_Cp_min, s.n_u_Cp_max);
fprintf('    dispersion over the K_s estimates (C_shared): disp_omega mean %.3e, median %.3e, max %.3e | disp_kappa mean %.3e, median %.3e, max %.3e\n', ...
        s.disp_omega_mean, s.disp_omega_median, s.disp_omega_max, ...
        s.disp_kappa_mean, s.disp_kappa_median, s.disp_kappa_max);
fprintf('    VALIDITY: valid %d of %d | bpd_wb_fail %d | A_exc %d | A_nonfinite %d | C_exc %d | Cp_bpd_fail %d | Cp_exc %d | C_nonfinite %d | Cp_nonfinite %d | L_nonfinite %d | sel_unexpected %d\n', ...
        s.n_valid, s.n, s.sum_v_bpd_wb_fail, s.sum_v_A_exc, s.sum_v_A_nonfinite, ...
        s.sum_v_C_exc, s.sum_v_Cp_bpd_fail, s.sum_v_Cp_exc, s.sum_v_C_nonfinite, ...
        s.sum_v_Cp_nonfinite, s.sum_v_L_nonfinite, s.sum_v_sel_unexpected);
fprintf('    time: %.1f s for the set; %.3f s per trial (mean of the per-trial totals)\n', ...
        s.elapsed_s, s.rt_trial_mean_s);
end

function [m, sd] = loc_t54_msd(x)
%LOC_T54_MSD  Mean and SD of the finite entries; NaN if there is none.
x = x(isfinite(x));
if isempty(x)
    m = NaN; sd = NaN;
else
    m = mean(x); sd = std(x);
end
end

function f = loc_t54_fracC(L)
%LOC_T54_FRACC  Fraction of rows with L_C <= min(L_A, L_B), over the rows
%  whose three objective values are finite; NaN if there is none.
ok = all(isfinite(L), 2);
if ~any(ok)
    f = NaN;
else
    f = mean(L(ok, 3) <= min(L(ok, 1), L(ok, 2)));
end
end

function x = loc_t54_col(M, kind, nm)
%LOC_T54_COL  One named column of a T-54 numeric row matrix.
names = t54_cols(kind);
j = find(strcmp(names, nm));
assert(numel(j) == 1, 't40_e2b_driver: T-54 column %s not found exactly once.', nm);
x = M(:, j);
end

function X = loc_t54_sub(SC, nm)
%LOC_T54_SUB  One named per-subcarrier column as a K_s x n matrix.
names = t54_cols('sub');
j = find(strcmp(names, nm));
assert(numel(j) == 1, 't40_e2b_driver: T-54 per-subcarrier column %s not found exactly once.', nm);
X = reshape(SC(:, j, :), size(SC, 1), size(SC, 3));
end

function T = loc_t54_trial_table(TR, C, setname)
%LOC_T54_TRIAL_TABLE  Per-trial rows of one trial set (change T54-1, arm A).
T = loc_t54_decode(array2table(TR, 'VariableNames', t54_cols('trial')));
T = [loc_t54_ctx_table(C, setname, size(TR, 1)), T];
end

function T = loc_t54_sub_table(SC, seeds, C, setname)
%LOC_T54_SUB_TABLE  Per-trial and per-subcarrier rows of one trial set
%  (change T54-1, arms C_shared and C_percarrier).
K_s = size(SC, 1);
NS  = size(SC, 2);
n   = size(SC, 3);
M   = reshape(permute(SC, [1 3 2]), K_s * n, NS);
T   = loc_t54_decode(array2table(M, 'VariableNames', t54_cols('sub')));
X   = loc_t54_ctx_table(C, setname, K_s * n);
X.seed = kron(seeds(:), ones(K_s, 1));
T   = [X, T];
end

function X = loc_t54_ctx_table(C, setname, n)
%LOC_T54_CTX_TABLE  The context columns of n CSV rows.
X = table(repmat(string(C.run_tag), n, 1), repmat(string(C.slab), n, 1), ...
          repmat(string(C.gtag), n, 1), repmat(double(C.anchor_off), n, 1), ...
          repmat(string(C.geom), n, 1), repmat(string(C.leg_label_csv), n, 1), ...
          repmat(double(C.is_reg), n, 1), repmat(C.B_hz, n, 1), repmat(C.K_s, n, 1), ...
          repmat(C.SNR_dB, n, 1), repmat(string(setname), n, 1), ...
          'VariableNames', {'run_tag', 'stage_label', 'grid', 'anchor_off', 'leg', ...
          'leg_label', 'registered_cell', 'B_hz', 'K_s', 'SNR_dB', 'trial_set'});
end

function T = loc_t54_decode(T)
%LOC_T54_DECODE  Replaces every *_sel code column by a *_PhaseD_select
%  column of single characters: 0 '?', 1 'A', 2 'B', 3 'C', 4 'a'.
lab = '?ABCa';
vn  = T.Properties.VariableNames;
for j = 1:numel(vn)
    nm = vn{j};
    if numel(nm) > 4 && strcmp(nm(end-3:end), '_sel')
        c = T.(nm);
        c(~isfinite(c) | c < 0 | c > 4) = 0;
        T.(nm) = cellstr(reshape(lab(c + 1), [], 1));
        T.Properties.VariableNames{j} = [nm(1:end-4) '_PhaseD_select'];
    end
end
end

function loc_t54_write(T, fpath)
%LOC_T54_WRITE  Writes a table to a CSV: creates it, or appends rows to it.
if exist(fpath, 'file') == 2
    writetable(T, fpath, 'WriteMode', 'append');
else
    writetable(T, fpath);
end
end

function fpath = loc_t54_file(F, kind)
%LOC_T54_FILE  Path of one output file of a T-54 stage call:
%  t54_<run tag>_<stage>_<grid>_<kind>_<stamp>.csv in out_dir.
fpath = [F.pre '_' kind '_' F.stamp '.csv'];
end
