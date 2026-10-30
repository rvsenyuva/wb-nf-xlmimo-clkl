function out = ctA_sweep_paperC(P, snr_vec, N_MC, ref_dump, out_dir, stamp, cfg)
%CTA_SWEEP_PAPERC  The CT-A estimator on trials of Paper C's production model.
%
%   OUT = CTA_SWEEP_PAPERC(P, SNR_VEC, N_MC, REF_DUMP, OUT_DIR, STAMP, CFG)
%
%   Adapted from matlab/ctA2/ctA2_sweep.m of rvsenyuva/wb-nf-crb-globecom26
%   at commit ba1862cb39aa638c0dbc700f8f310e3e7dd12dcf, blob
%   4b5deb7241cbabba2ded2e0ba15ff175aad63c88. Each hunk against that script
%   and its reason is in CTA_EDITS.md, and the registration in CTA_SPEC.md,
%   of sessions/2026-10-29_CTA-spec in paperC-lab.
%
%   DRAFT. This file was written as text. It has not been parsed or run.
%
%   The estimator and its helpers are the unchanged copies in this folder.
%   A trial is the trial of the production chain with the same seed and SNR:
%   rng(mc_seed, 'twister'), then one two-argument call of
%   wb_channel_gen_ofdm_nf, as run_one_realisation does at L449 and L464-465
%   of run_monte_carlo_paperC.m (commit 1e5f9506).
%
%   INPUTS
%     P         production struct, built as the run sheet builds it.
%     snr_vec   SNR values [dB].
%     N_MC      seeds per SNR: mc_seed = P.mc_seed_base + (1:N_MC).
%     ref_dump  per-trial CSV the scene gate G-03 compares with.
%     out_dir   folder of the two CSV files.
%     stamp     text in the two file names.
%     cfg       no defaults; the run sheet sets every field.
%       .stage         'S0' (G-07 and G-06 only) or 'RUN'
%       .npar          workers of the parfor over seeds
%       .N_CAND        candidates carried to stage 2
%       .K_COARSE      detection subcarriers
%       .setup_expect  1 x 4 cell of text: r_lo, r_hi, u_min, u_max
%
%   OUTPUTS, in out_dir
%     ctA_trials_<stamp>.csv   one row per (seed, SNR), appended after every
%                              chunk, so a fault keeps what was written
%     ctA_summary_<stamp>.csv  one row per SNR, written at the end
%
%   The console carries gate lines, set-up lines and elapsed times. It carries
%   no statistic of the run.
%
%   No bound is computed here. No error is compared with a bound. No value of
%   a trial's truth reaches the estimate path: eta_true is read by
%   wb_nf_music_est only after the estimate is fixed (its L51-55, L240-250).

narginchk(7, 7);
t_all = tic;
out   = struct('ok', false, 'fn_trials', '', 'fn_summary', '');

try

assert(P.d == 1, 'ctA_sweep_paperC: the estimator takes one path');
assert(~P.use_quadriga_channel && ~P.use_riviello_snr_axis, ...
    'ctA_sweep_paperC: only the primary generator path is registered');
assert(~P.wb_gen_write_csv, ...
    'ctA_sweep_paperC: the generator must not write its CSV from a worker');
assert(all(isfield(cfg, {'stage', 'npar', 'N_CAND', 'K_COARSE', 'setup_expect'})), ...
    'ctA_sweep_paperC: cfg has no defaults');
assert(any(strcmp(cfg.stage, {'S0', 'RUN'})), ...
    'ctA_sweep_paperC: cfg.stage is S0 or RUN');

L = ctA1_locals();

%% ========================================================================
%  SETUP -- the constants of ctA2_sweep.m L136-138 at Paper C's geometry
% =========================================================================
q_scale = 2*pi*P.d_ant/P.lambda_c;
s0      = pi*P.d_ant^2/P.lambda_c;
m_bar   = ((0:P.M-1) - (P.M-1)/2).';

[GRID, gi] = loc_grid(P, s0, m_bar, cfg);

NBOOT = 2000;
MAXIT = 400;

OPT0 = struct('s0', s0, 'q', q_scale, 'npar', 0, ...
              'batch', numel(snr_vec), 'assert_internal', true, ...
              'tangent', true, ...
              'tol_dec', 2e-12, 'tol_ml', 1e-9, 'maxit', MAXIT);
assert(OPT0.tangent, 'OPT.tangent must be true for every number in the figure');
assert(OPT0.maxit == MAXIT, 'maxit must be the uniform production budget');

fprintf('CTA-CONFIG M = %d, N_RF = %d, K_s = %d, N = %d, d = %d\n', ...
        P.M, P.N_RF, P.K_s, P.N, P.d);
fprintf('CTA-CONFIG grid %d x %d, steps %.6g and %.6g mainlobes\n', ...
        gi.n_om, gi.n_ka, gi.step_om_ml, gi.step_ka_ml);
fprintf('CTA-CONFIG K_COARSE = %d, N_CAND = %d, candidate share of the grid %.6g, maxit = %d\n', ...
        cfg.K_COARSE, cfg.N_CAND, gi.cand_share, MAXIT);
fprintf('CTA-CONFIG SNR values %d, N_MC %d, seed base %d, workers %d, stage %s\n', ...
        numel(snr_vec), N_MC, P.mc_seed_base, cfg.npar, cfg.stage);

%% ========================================================================
%  G-07 -- the set-up against the sheet, text domain
% =========================================================================
got07 = {sprintf('%.4f', P.r_lo_fac*P.r_RD), sprintf('%.4f', P.r_hi_fac*P.r_RD), ...
         sprintf('%.6g', P.u_min), sprintf('%.6g', P.u_max)};
ok07  = isequal(got07, cfg.setup_expect);
fprintf('CTA-G07 [text of %%.4f and %%.6g; r_lo r_hi u_min u_max against the sheet]: %s %s %s %s: %s\n', ...
        got07{1}, got07{2}, got07{3}, got07{4}, loc_vs(ok07));
if ~ok07
    error('ctA:G07', 'G-07 FAIL');
end

%% ========================================================================
%  G-06 -- self-test before any trial
% =========================================================================
loc_selftest(P, snr_vec(1), GRID, OPT0, q_scale, s0, m_bar, L);

if strcmp(cfg.stage, 'S0')
    out.ok = true;
    fprintf('CTA-S0-DONE elapsed %.1f s\n', toc(t_all));
    return;
end

%% ========================================================================
%  SWEEP -- one estimator call per seed, on the stack of its SNR values
% =========================================================================
ref       = loc_read_ref(ref_dump);
n_snr     = numel(snr_vec);
seed_base = P.mc_seed_base;
npar      = cfg.npar;
CHUNK     = max(1, npar);
n_chunk   = ceil(N_MC/CHUNK);

ctx = struct('GRID', GRID, 'OPT0', OPT0, 'q_scale', q_scale, 's0', s0, ...
             'm_bar', m_bar);

fn_t = fullfile(out_dir, ['ctA_trials_' stamp '.csv']);
fn_s = fullfile(out_dir, ['ctA_summary_' stamp '.csv']);
loc_trials_header(fn_t);
out.fn_trials  = fn_t;
out.fn_summary = fn_s;

R     = cell(n_snr, 1);       % per-SNR accumulated results
n_cmp = 0;
n_dif = 0;

for ic = 1:n_chunk
    t_c = tic;
    idx = ((ic-1)*CHUNK + 1) : min(ic*CHUNK, N_MC);
    res = cell(numel(idx), 1);

    parfor (j = 1:numel(idx), npar)
        res{j} = loc_one_seed(P, snr_vec, seed_base, idx(j), ctx);
    end

    fid = fopen(fn_t, 'a');
    assert(fid > 0, 'ctA_sweep_paperC: cannot open %s', fn_t);
    for jr = 1:numel(idx)
        T = res{jr};
        for s = 1:n_snr
            want = ref(sprintf('%.6g|%d', snr_vec(s), T.mc_seed));
            got  = {sprintf('%.10g', T.th0*180/pi), sprintf('%.10g', T.r0)};
            n_cmp = n_cmp + 2;
            n_dif = n_dif + ~strcmp(got{1}, want{1}) + ~strcmp(got{2}, want{2});

            Q = loc_take(T.pack, s);
            Q.th0     = T.th0;
            Q.r0      = T.r0;
            Q.mc_seed = T.mc_seed;
            R{s} = loc_cat(R{s}, Q);
            loc_trial_row(fid, snr_vec(s), seed_base, T, Q, P.M, q_scale);
        end
    end
    fclose(fid);

    if n_dif > 0
        fprintf('CTA-G03 [text of %%.10g; theta_true_deg and r_true_m against the dump]: cells compared %d, cells that differ %d: FAIL\n', ...
                n_cmp, n_dif);
        error('ctA:G03', 'G-03 FAIL at chunk %d', ic);
    end
    fprintf('CTA-TIME chunk %d of %d, %.1f s\n', ic, n_chunk, toc(t_c));
end

fprintf('CTA-G03 [text of %%.10g; theta_true_deg and r_true_m against the dump]: cells compared %d, cells that differ %d: PASS\n', ...
        n_cmp, n_dif);
fprintf('CTA-G04 [bit identity of theta, r and W across the SNR values of a seed]: seeds %d, SNR values %d, assertion raised 0: PASS\n', ...
        N_MC, n_snr);
fprintf('CTA-G05 [wb_nf_music_est assertions, stage-2 against stage-3 cost below 1e-10]: calls %d, assertion raised 0: PASS\n', ...
        N_MC);

%% ========================================================================
%  G-02 -- seeds and trial counts, integer domain
% =========================================================================
ok02 = true;
for s = 1:n_snr
    ok02 = ok02 && isequal(sort(R{s}.mc_seed), seed_base + (1:N_MC).');
end
fprintf('CTA-G02 [integer; each SNR holds exactly the seeds base+1 to base+N_MC]: SNR values %d, seeds %d, base %d: %s\n', ...
        n_snr, N_MC, seed_base, loc_vs(ok02));
if ~ok02
    error('ctA:G02', 'G-02 FAIL');
end

%% ========================================================================
%  SUMMARY
% =========================================================================
KV = {};
for s = 1:n_snr
    stb = RandStream('mt19937ar', 'Seed', 70000 + s);
    ix  = randi(stb, N_MC, N_MC, NBOOT);
    S   = loc_stats(R{s}, P.M, q_scale, ix, L);
    KV  = loc_row(KV, snr_vec(s), S, cfg.K_COARSE, cfg.N_CAND, MAXIT, gi);
end
loc_flush(KV, fn_s);

out.ok = true;
fprintf('CTA-DONE files written, elapsed %.1f s\n', toc(t_all));

catch ME
    fprintf('CTA-HALT: %s\n', ME.message);
    rethrow(ME);
end

end

%% ========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function T = loc_one_seed(P, snr_vec, seed_base, mc_idx, ctx)
% One seed: its trial at every SNR, one estimator call on the stack. W and the
% scene are the same at every SNR of a seed (Ledger F-134), which G-04 checks.
L       = ctA1_locals();
n_snr   = numel(snr_vec);
mc_seed = seed_base + mc_idx;
Ustack  = zeros(P.N_RF, P.K_s, n_snr);
th0 = [];  r0 = [];  W = [];  Tw = [];

for s = 1:n_snr
    rng(mc_seed, 'twister');
    [~, Y_full, ~, theta_gen, r_gen, ~, ~, W_comb] = ...
        wb_channel_gen_ofdm_nf(P, snr_vec(s));
    if s == 1
        th0 = theta_gen;
        r0  = r_gen;
        W   = W_comb;
        [~, Tw] = L.whiten(W);
    else
        assert(isequal(theta_gen, th0) && isequal(r_gen, r0) && isequal(W_comb, W), ...
            'G-04 FAIL at seed %d', mc_seed);
    end
    Ustack(:,:,s) = loc_eigvecs(Y_full, Tw, P.N);
end

om0 = ctx.q_scale*cos(th0);
ka0 = ctx.s0*sin(th0)^2/r0;

O = ctx.OPT0;
O.eta_true = [om0 ka0];
[eta, info] = wb_nf_music_est(Ustack, W, Tw, P.alpha_k_vec, ctx.m_bar, ctx.GRID, O);

T = struct();
T.pack    = loc_pack(eta, info);
T.mc_idx  = mc_idx;
T.mc_seed = mc_seed;
T.th0     = th0;
T.r0      = r0;
T.om0     = om0;
T.ka0     = ka0;
end

function U = loc_eigvecs(Y_full, Tw, N_snap)
% The lines of wb_nf_gen_eigvecs.m L103-109 on compressed data of the
% production chain, after the whitener of ctA1_locals.m L108-114.
Ks = size(Y_full, 3);
U  = zeros(size(Tw, 1), Ks);
for k = 1:Ks
    Y = Tw*Y_full(:,:,k);
    R = (Y*Y')/N_snap;                  % NO sample-mean removal (9.4.2)
    [Ve, De] = eig((R+R')/2);
    dv = real(diag(De));
    [~, ord] = sort(dv, 'descend');

    uh     = Ve(:, ord(1));
    U(:,k) = uh / norm(uh);
end
end

function [GRID, gi] = loc_grid(P, s0, m_bar, cfg)
% The grids of ctA2_sweep.m L155-161 at Paper C's geometry. The steps of that
% script, in mainlobes of its array, are kept; the half-width of the capture
% window at wb_nf_music_est.m L244 is in cells, so the steps set its width in
% mainlobes. The kappa box comes from P, never from a trial's truth.
M_G       = 256;                                 % ctA2_sweep.m L133
c_chirp_G = s0*sin(40*pi/180)^2;                 % L141, L143; s0 is the same
d_om_G    = 2*pi/(1536 - 1);                     % L157
d_ka_G    = (c_chirp_G/0.5 - c_chirp_G/200)/(288 - 1);    % L158-159

step_om_ml = d_om_G/(2*pi/M_G);                  % wb_nf_music_est.m L127
step_ka_ml = d_ka_G/(pi/((M_G - 1)/2)^2);        % wb_nf_music_est.m L128

ml_om  = 2*pi/P.M;
ml_ka  = pi/max(m_bar.^2);
ka_min = s0*sin(P.theta_lo)^2*P.u_min;
ka_max = s0*sin(P.theta_hi)^2*P.u_max;

n_om_all = ceil(2*pi/(step_om_ml*ml_om)) + 1;
n_ka     = ceil((ka_max - ka_min)/(step_ka_ml*ml_ka)) + 1;

om_grid = linspace(-pi, pi, n_om_all).';  om_grid = om_grid(2:end-1);
ka_grid = linspace(ka_min, ka_max, n_ka).';
GRID    = struct('om_grid', om_grid, 'ka_grid', ka_grid, ...
                 'idx_c', round(linspace(1, P.K_s, cfg.K_COARSE)), ...
                 'N_CAND', cfg.N_CAND);

gi = struct();
gi.n_om       = numel(om_grid);
gi.n_ka       = numel(ka_grid);
gi.step_om_ml = (om_grid(2) - om_grid(1))/ml_om;
gi.step_ka_ml = (ka_grid(2) - ka_grid(1))/ml_ka;
gi.cand_share = cfg.N_CAND/(gi.n_om*gi.n_ka);
end

function loc_selftest(P, snr, GRID, OPT0, q_scale, s0, m_bar, L)
% G-06. Noise-free, one scene, no trial of the run. (a) loc_eigvecs against
% wb_nf_gen_eigvecs on the same vectors. (b) wb_nf_music_est flags capture at
% the mid-box pair of run_monte_carlo_paperC.m L217-219 on those vectors.
rng(P.mc_seed_base + 1, 'twister');
[~, ~, ~, ~, ~, ~, ~, W] = wb_channel_gen_ofdm_nf(P, snr);
[~, Tw] = L.whiten(W);

th_nom = (P.theta_lo + P.theta_hi)/2;
r_nom  = P.r_lo_fac*P.r_RD + (P.r_hi_fac - P.r_lo_fac)*P.r_RD/2;
om_nom = q_scale*cos(th_nom);
ka_nom = s0*sin(th_nom)^2/r_nom;

b     = L.bvec(W, Tw, P.alpha_k_vec, m_bar, m_bar.^2, om_nom, ka_nom);
U_ref = wb_nf_gen_eigvecs(b, 1, 0, P.N);

s_row = exp(1j*2*pi*(0:P.N-1)/P.N);              % wb_nf_gen_eigvecs.m L73
Y_c   = zeros(P.N_RF, P.N, P.K_s);
for k = 1:P.K_s
    Y_c(:,:,k) = Tw \ (b(:,k)*s_row);
end
U_new = loc_eigvecs(Y_c, Tw, P.N);

err_u = max(abs(1 - abs(sum(conj(U_ref).*U_new, 1))));
ok_a  = err_u < 1e-10;                           % wb_nf_music_est.m L231

O = OPT0;
O.batch    = 1;
O.eta_true = [om_nom ka_nom];
[~, info]  = wb_nf_music_est(U_new, W, Tw, P.alpha_k_vec, m_bar, GRID, O);
ok_b = info.captured;

fprintf('CTA-G06 [real, tolerance 1e-10, and boolean; noise-free]: eigenvector difference %.3e, capture flag %d: %s\n', ...
        err_u, ok_b, loc_vs(ok_a && ok_b));
if ~(ok_a && ok_b)
    error('ctA:G06', 'G-06 FAIL');
end
end

function ref = loc_read_ref(fn)
% Reads, from the per-trial dump, the two columns of the scene gate and the
% two of its key, as text. No other column is kept.
txt = fileread(fn);
ln  = regexp(txt, '\r?\n', 'split');
if isempty(ln{end}); ln(end) = []; end
hdr = strsplit(ln{1}, ',');
c_s = find(strcmp(hdr, 'SNR_dB'));
c_m = find(strcmp(hdr, 'mc_seed'));
c_t = find(strcmp(hdr, 'theta_true_deg'));
c_r = find(strcmp(hdr, 'r_true_m'));
assert(isscalar(c_s) && isscalar(c_m) && isscalar(c_t) && isscalar(c_r), ...
    'ctA_sweep_paperC: the reference dump lacks a column of the gate');
ref = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 2:numel(ln)
    f = strsplit(ln{i}, ',');
    ref(sprintf('%s|%s', f{c_s}, f{c_m})) = {f{c_t}, f{c_r}};
end
end

function P = loc_pack(eta, info)
% Fields are assigned INDIVIDUALLY, never through the struct() constructor.
% struct('f', c) with a cell-valued c builds a STRUCT ARRAY of size(c)
% rather than a scalar struct holding the cell -- which would silently
% shred stop_reason and verdict into 600 one-element structs.
P = struct();
P.eta           = eta;
P.r_hat         = info.r_hat;
P.theta_hat_deg = info.theta_hat_deg;
P.n_iter        = info.n_iter;
P.lambda2       = info.lambda2;
P.g_inf         = info.g_inf;
P.cond_min      = info.cond_Hs_min;
P.cond_max      = info.cond_Hs_max;
P.step_ml       = info.step_ml;
P.acc_frac      = info.acc_frac;
P.cert_pass     = info.cert_pass;
P.boundary      = info.boundary;
P.illcond       = info.illcond;
P.nonstationary = info.nonstationary;
P.captured      = info.captured;
P.rank          = info.rank;
P.Chat_coarse   = info.Chat_coarse;
P.Chat_final    = info.Chat_final;
P.stop_reason   = info.stop_reason;      % n x 1 cell, assigned directly
P.verdict       = info.verdict;          % n x 1 cell, assigned directly
end

function Q = loc_take(A, s)
% Row s of every field of a pack.
Q = struct();
f = fieldnames(A);
for i = 1:numel(f)
    Q.(f{i}) = A.(f{i})(s,:);
end
end

function A = loc_cat(A, B)
if isempty(A); A = B; return; end
f = fieldnames(B);
for i = 1:numel(f)
    A.(f{i}) = [A.(f{i}); B.(f{i})];
end
end

function [e_r, e_th, gross] = loc_err(r_hat, th_hat_deg, r0, th0, M, q_scale)
% Errors with the sign of run_monte_carlo_paperC.m L730-731. Gross-error flag
% of ctA2_sweep.m L421, on each trial's own truth.
e_r   = r0 - r_hat;
e_th  = th0*180/pi - th_hat_deg;
g_th  = -1./(q_scale*sin(th0));                  % ctA2_sweep.m L168
ml_th_deg = (2*pi/M)*abs(g_th)*180/pi;           % L169-170
gross = (abs(e_th) > ml_th_deg/2) | (r_hat < r0/2) | (r_hat > 2*r0);
end

function S = loc_stats(R, M, q_scale, ix, L)
% All statistics at one SNR. NOTHING here tunes anything.
n = numel(R.r_hat);
S = struct();

[e_r, e_th, gross] = loc_err(R.r_hat, R.theta_hat_deg, R.r0, R.th0, M, q_scale);

S.RMSE_r     = sqrt(mean(e_r.^2));
S.RMSE_theta = sqrt(mean(e_th.^2));
S.ci_RMSE_r     = loc_ci_rmse(e_r, ix, L);
S.ci_RMSE_theta = loc_ci_rmse(e_th, ix, L);

% ---- Given no gross error (registered in CTA_SPEC.md Sec. 9, MS-05) ------
keep = ~gross;
S.n_not_gross = sum(keep);
if any(keep)
    S.RMSE_r_ng     = sqrt(mean(e_r(keep).^2));
    S.RMSE_theta_ng = sqrt(mean(e_th(keep).^2));
else
    S.RMSE_r_ng = NaN;  S.RMSE_theta_ng = NaN;
end

% ---- Health, taxonomy, detection ---------------------------------------
S.capture       = mean(R.captured);
S.cert_pass     = mean(R.cert_pass);
S.nonstat_rate  = mean(R.nonstationary);
S.boundary_rate = mean(R.boundary);
S.illcond_rate  = mean(R.illcond);
S.gross_rate    = mean(gross);
S.n_iter_med    = median(R.n_iter);
S.n_iter_p95    = L.prctile(R.n_iter, 95);
S.n_iter_max    = max(R.n_iter);
S.n_CONVERGED   = sum(strcmp(R.stop_reason, 'CONVERGED'));
S.n_MAXIT       = sum(strcmp(R.stop_reason, 'MAXIT'));
S.n_STALL       = sum(strcmp(R.stop_reason, 'STALL'));
S.n_OK          = sum(strcmp(R.verdict, 'OK'));
S.n_v_MAXIT     = sum(strcmp(R.verdict, 'MAXIT'));
S.n_v_BOUNDARY  = sum(strcmp(R.verdict, 'BOUNDARY'));
S.n_v_STALL     = sum(strcmp(R.verdict, 'STALL'));
S.n_v_ILLCOND   = sum(strcmp(R.verdict, 'ILLCOND'));
S.cond_max      = max(R.cond_max);

rk = R.rank(R.captured);
if isempty(rk)
    S.rank_p50 = NaN; S.rank_p95 = NaN; S.rank_p99 = NaN; S.rank_max = NaN;
else
    S.rank_p50 = median(rk);        S.rank_p95 = L.prctile(rk, 95);
    S.rank_p99 = L.prctile(rk, 99); S.rank_max = max(rk);
end

% Open-marker flag of ctA2_sweep.m L478: a diagnostic column, no reading.
S.open_marker = double(S.capture < 0.99);
S.N_MC = n;
end

function ci = loc_ci_rmse(err, ix, L)
e2 = err(:).^2;
rb = sqrt(mean(reshape(e2(ix), size(ix,1), size(ix,2)), 1));
ci = [L.prctile(rb, 2.5), L.prctile(rb, 97.5)];
end

function KV = loc_row(KV, SNR, S, KC, NC, mx, gi)
KV(end+1,:) = {SNR, S.N_MC, KC, NC, mx, gi.n_om, gi.n_ka, ...
    S.RMSE_theta, S.ci_RMSE_theta(1), S.ci_RMSE_theta(2), ...
    S.RMSE_r, S.ci_RMSE_r(1), S.ci_RMSE_r(2), ...
    S.n_not_gross, S.RMSE_theta_ng, S.RMSE_r_ng, ...
    S.capture, S.rank_p50, S.rank_p95, S.rank_p99, S.rank_max, ...
    S.n_CONVERGED, S.n_MAXIT, S.n_STALL, ...
    S.n_OK, S.n_v_MAXIT, S.n_v_BOUNDARY, S.n_v_STALL, S.n_v_ILLCOND, ...
    S.cert_pass, S.nonstat_rate, S.boundary_rate, S.illcond_rate, S.gross_rate, ...
    S.n_iter_med, S.n_iter_p95, S.n_iter_max, S.cond_max, S.open_marker};
end

function loc_trials_header(fn)
vn = {'SNR_dB','mc_seed_base','mc_idx','mc_seed','theta_true_deg','r_true_m', ...
    'omega_true','kappa_true','omega_hat','kappa_hat','r_hat_m','theta_hat_deg', ...
    'err_r_m','err_theta_deg','Chat_coarse','Chat_final','n_iter','lambda2', ...
    'g_inf','cond_Hs_min','cond_Hs_max','step_ml','acc_frac','cert_pass', ...
    'boundary','illcond','nonstationary','captured','rank','gross', ...
    'stop_reason','verdict'};
fid = fopen(fn, 'w');
assert(fid > 0, 'ctA_sweep_paperC: cannot open %s', fn);
fprintf(fid, '%s\n', strjoin(vn, ','));
fclose(fid);
end

function loc_trial_row(fid, snr, seed_base, T, Q, M, q_scale)
[e_r, e_th, gross] = loc_err(Q.r_hat, Q.theta_hat_deg, T.r0, T.th0, M, q_scale);
fprintf(fid, '%.6g,%d,%d,%d,%.10g,%.10g,', ...
        snr, seed_base, T.mc_idx, T.mc_seed, T.th0*180/pi, T.r0);
fprintf(fid, '%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,', ...
        T.om0, T.ka0, Q.eta(1), Q.eta(2), Q.r_hat, Q.theta_hat_deg, e_r, e_th);
fprintf(fid, '%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,', ...
        Q.Chat_coarse, Q.Chat_final, Q.n_iter, Q.lambda2, Q.g_inf, ...
        Q.cond_min, Q.cond_max, Q.step_ml, Q.acc_frac);
fprintf(fid, '%d,%d,%d,%d,%d,%g,%d,%s,%s\n', ...
        Q.cert_pass, Q.boundary, Q.illcond, Q.nonstationary, Q.captured, ...
        Q.rank, gross, Q.stop_reason{1}, Q.verdict{1});
end

function loc_flush(KV, fn)
vn = {'SNR_dB','N_MC','K_COARSE','N_CAND','maxit','n_om','n_ka', ...
  'RMSE_theta_deg','RMSE_theta_ci_lo','RMSE_theta_ci_hi', ...
  'RMSE_r_m','RMSE_r_ci_lo','RMSE_r_ci_hi', ...
  'n_not_gross','RMSE_theta_not_gross_deg','RMSE_r_not_gross_m', ...
  'capture','rank_p50','rank_p95','rank_p99','rank_max', ...
  'n_CONVERGED','n_MAXIT','n_STALL', ...
  'v_OK','v_MAXIT','v_BOUNDARY','v_STALL','v_ILLCOND', ...
  'cert_pass_rate','nonstationary_rate','boundary_rate','illcond_rate', ...
  'gross_rate','n_iter_med','n_iter_p95','n_iter_max','cond_Hs_max', ...
  'open_marker'};
writetable(cell2table(KV, 'VariableNames', vn), fn);
end

function s = loc_vs(ok)
if ok; s = 'PASS'; else; s = 'FAIL'; end
end
