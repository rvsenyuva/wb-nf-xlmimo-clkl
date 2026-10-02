function t40_e2b_driver(stage, out_dir, grid_mode)
%T40_E2B_DRIVER  T-40 / E2b discriminating experiment (Paper C strengthening).
%
%  t40_e2b_driver(stage)
%  t40_e2b_driver(stage, out_dir)
%  t40_e2b_driver(stage, out_dir, grid_mode)     grid_mode = 'legacy'|'refined'
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
diary(fullfile(out_dir, sprintf('t40_%s_%s_%s_console.txt', ...
      lower(stage), gtag, stamp)));
diary on;
fprintf('=============================================================\n');
fprintf('  T-40 / E2b  v2   stage = %s   grid = %s   %s\n', stage, gtag, stamp);
fprintf('=============================================================\n');

switch stage
    case 'P0_BENCH';    t40_bench(CFG, out_dir, stamp, gtag);
    case 'P4_REPRO';    t40_repro(CFG, out_dir, stamp, gtag);
    case 'P1_FLOOR';    t40_floor(CFG, out_dir, stamp, gtag);
    case 'P2_SELFTEST'; t40_selftest(CFG, out_dir, stamp, gtag);
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
end


% =========================================================================
%  P STRUCT
% =========================================================================
function P = t40_P(CFG, B_hz)
%T40_P  Self-sufficient parameter struct. Any field the live
%       setup_production_P[_v4] does not supply is filled from the locked
%       set; every derived quantity is recomputed; the box is pinned here.

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
end
if ~isstruct(P); P = struct(); end

P = loc_fill_P(P, CFG);

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


function P = loc_fill_P(P, CFG)
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
function t40_selftest(CFG, out_dir, stamp, gtag)
P = t40_P(CFG, 400e6);
G = {};
fprintf('\n### P2 self-tests ###\n');

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
R1 = t40_trial(Pk, CFG, CFG.SNR_star, CFG.seed_base + 1, 'pinned', false);
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
R2 = t40_trial(P, CFG, CFG.SNR_star, CFG.seed_base + 2, 'pinned', false);
G(end+1,:) = {'S5_shared_scene', abs(R2.r_true - CFG.pin_r), 1e-12}; %#ok<AGROW>

% S6 (FIX 5): the per-subcarrier estimates must NOT be degenerate. This is a
% capability check on the EXPERIMENT, not on the code: if all K_s estimates
% are bit-identical, the rival arm receives no data diversity and the A-vs-C
% comparison cannot answer the question at this grid.
nd = 12; deg = zeros(nd,1);
for i = 1:nd
    Rd = t40_trial(P, CFG, CFG.SNR_star, CFG.seed_base + 500 + i, ...
                   'pinned', false);
    deg(i) = Rd.degen_C;
end
G(end+1,:) = {'S6_C_not_degenerate', mean(deg), 0.50}; %#ok<AGROW>

ok = true;
for i = 1:size(G,1)
    p = G{i,2} < G{i,3};
    ok = ok && p;
    fprintf('  %-24s val=%.3e tol=%.1e  %s\n', G{i,1}, G{i,2}, G{i,3}, ...
            loc_tern(p, 'PASS', 'FAIL'));
end
if ~(mean(deg) < 0.50)
    fprintf('  S6 NOTE: degenerate per-subcarrier estimates. The grid scan\n');
    fprintf('  cannot resolve subcarrier-to-subcarrier variation, so the\n');
    fprintf('  rival arm gets no data diversity. A local off-grid polish is\n');
    fprintf('  required before the gate means anything. STOP and report.\n');
end
fprintf('  GATE P2: %s\n', loc_tern(ok, 'PASS', 'FAIL -- do not proceed'));
writetable(cell2table(G, 'VariableNames', {'test','value','tol'}), ...
    fullfile(out_dir, sprintf('t40_P2_selftest_%s_%s.csv', gtag, stamp)));
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
            DISP_O, DISP_K, DISP_OB, DEGC, DEGB, KARG, CONDB, ...
            'VariableNames', {'tag','grid','sweep_value','SNR_dB','B_hz', ...
            'K_s','geom','flat_alpha','trial','seed','arm','theta_true_deg', ...
            'r_true_m','r_hat_m','err_r_m','err_theta_deg','n_iter', ...
            'converged','gross','gross_legacy','outbox','boundary', ...
            'disp_omega','disp_kappa','disp_omega_B','degen_C','degen_B', ...
            'k_argmin_L','cond_B_max'})]; %#ok<AGROW>
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
            mean(RTIME(:),'omitnan'), ...
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
            'degen_C_frac','degen_B_frac','rt_mean_s'})]; %#ok<AGROW>
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
t0 = tic;
try
    [th_bpd, r_bpd] = bpd_baseline(X, P);
    th_bpd = th_bpd(1); r_bpd = r_bpd(1);
catch
    th_bpd = (P.theta_lo + P.theta_hi)/2;
    r_bpd  = 1/(0.5*(P.u_min + P.u_max));
end
R.rt_bpd = toc(t0);
p0 = ones(P.d,1)/P.d;

% ---- ARM A ---------------------------------------------------------------
t0 = tic;
[thA, rA, ~, niA, cvA] = loc_safe_clkl(Rh, W, th_bpd, r_bpd, p0, P);
R.rt_A = toc(t0);

r_hat = struct(); th_hat = struct(); n_it = struct(); conv = struct();
r_hat.A = rA; th_hat.A = thA; n_it.A = niA; conv.A = cvA;

if only_A
    R.rt_C = 0; R.rt_B = 0;
    R.disp_omega = NaN; R.disp_kappa = NaN; R.disp_omega_B = NaN;
    R.degen_C = NaN; R.degen_B = NaN; R.k_argmin_L = NaN; R.cond_B_max = NaN;
    R = loc_finish(R, r_hat, th_hat, n_it, conv, P, CFG);
    return
end

% ---- ARMS C -------------------------------------------------------------
t0 = tic;
omC  = nan(K_s,1); kaC  = nan(K_s,1); LkC = nan(K_s,1);
niC  = nan(K_s,1); cvC  = nan(K_s,1);
omCp = nan(K_s,1); kaCp = nan(K_s,1); niCp = nan(K_s,1); cvCp = nan(K_s,1);
for k = 1:K_s
    Pk = loc_Pk(P, k);
    [thk, rk, Lk, nik, cvk] = loc_safe_clkl(Rh(k), W, th_bpd, r_bpd, p0, Pk);
    [omC(k), kaC(k)] = loc_to_natural(thk, rk, P);
    LkC(k) = Lk; niC(k) = nik; cvC(k) = cvk;
    if CFG.do_C_percarrier
        try
            [tb, rb] = bpd_baseline(X(:,:,k), Pk);
            tb = tb(1); rb = rb(1);
        catch
            tb = th_bpd; rb = r_bpd;
        end
        [thk2, rk2, ~, ni2, cv2] = loc_safe_clkl(Rh(k), W, tb, rb, p0, Pk);
        [omCp(k), kaCp(k)] = loc_to_natural(thk2, rk2, P);
        niCp(k) = ni2; cvCp(k) = cv2;
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
[~, R.k_argmin_L] = min(LkC);
R.cond_B_max = max(condB);

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

function [th, r, L, ni, cv] = loc_safe_clkl(Rc, W, th0, r0, p0, P)
th = NaN; r = NaN; L = Inf; ni = NaN; cv = 0;
try
    [th, r, ~, ~, info] = wb_clkl_driver_pc(Rc, W, th0, r0, p0, P);
    th = th(1); r = r(1);
    if isfield(info,'L_hist') && ~isempty(info.L_hist); L = info.L_hist(end); end
    if isfield(info,'n_iter');    ni = info.n_iter;              end
    if isfield(info,'converged'); cv = double(info.converged);   end
catch
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
