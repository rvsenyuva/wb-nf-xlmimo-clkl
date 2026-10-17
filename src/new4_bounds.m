function new4_bounds(stage, results_dir, stamp, r_hi_fac_m1m3, r_hi_fac_m8)
%NEW4_BOUNDS  The bound computation of the NEW-4 run: every gate, every bound.
%
%  new4_bounds('failfast', results_dir, stamp, r_hi_fac_m1m3, r_hi_fac_m8)
%  new4_bounds('full',     results_dir, stamp, r_hi_fac_m1m3, r_hi_fac_m8)
%
%  Written by the NEW-4 code session of paperC-lab (sessions/2026-10-18_NEW4-CODE)
%  from sessions/2026-10-17_NEW4-SPEC/NEW4_SPEC.md as decided in DECISIONS_N4G.md,
%  with DECISIONS_N4C.md and DECISIONS_N4CB.md. Written as text; never parsed or
%  executed by its author (N6). Held at a commit of the code repository only after
%  the guide's check and the author's block (N4C-D2).
%
%  The function is a bound computation. It calls no estimator and no channel
%  generator, reads no estimator column, and writes or prints no error statistic,
%  no estimator output and no quotient or difference of an estimator statistic
%  and a bound (inherited-finding block of Ledger N-195).
%
%  Called functions: setup_production_P_v4, wb_crb_compressed (unchanged) and
%  built-ins. The only random-number call is inside the copied local function
%  loc_recon (W2).
%
%  stage 'failfast': the config print and the gate N4-G2 at one row (SNR 10 dB of
%                    the M1/M3 file, both columns); writes no file (W4, W15).
%  stage 'full'    : N4-G2 at every row, N4-G3, the per-trial bounds with N4-G4 and
%                    N4-G5, the aggregation, N4-G6, the quantiles of N-34, the
%                    single-geometry bound of M8, then the files (W15).
%
%  results_dir, stamp and the two box values are fixed by the final run sheet and
%  passed in (W1; roadmap P-9; N4G-D8 (iii)). Every file written is in results_dir
%  and is named new4_<kind>_<stamp>.csv.
%
%  Gates N4-G0, N4-G1 and N4-G7 are the run sheet's Git Bash block, not this
%  function (N4G-D11; OPEN-1). The function prints nothing about its own hash.
%
%  Labels in comments: Wn = requirement of the starter of the code session;
%  N4-Gn, N4-In, On = items of NEW4_SPEC.md Sec. f and g; N4G-Dn, N4C-Dn,
%  N4CB-Dn = author decisions; OPEN-n = a choice no record fixes, listed in
%  OPEN_NEW4.md.

% =========================================================================
%  W1; N4-I1; N4G-D11, N4G-D17: arguments, stage, results folder
% =========================================================================
if nargin ~= 5
    error('new4_bounds:usage', ...
        'new4_bounds(stage, results_dir, stamp, r_hi_fac_m1m3, r_hi_fac_m8): five arguments are required.');
end
if ~ischar(stage) || ~(strcmp(stage, 'failfast') || strcmp(stage, 'full'))
    error('new4_bounds:stage', 'stage must be ''failfast'' or ''full''.');
end
if ~ischar(results_dir) || isempty(results_dir) || ~ischar(stamp) || isempty(stamp)
    error('new4_bounds:args', 'results_dir and stamp must be non-empty character vectors.');
end
box_ok = @(v) isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v > 0;
if ~box_ok(r_hi_fac_m1m3) || ~box_ok(r_hi_fac_m8)
    error('new4_bounds:r_hi_fac', ...
        'r_hi_fac_m1m3 and r_hi_fac_m8 must be positive finite real scalars (the run sheet gives them with their source, F-130).');
end
is_full = strcmp(stage, 'full');

% =========================================================================
%  Named constants: every input path is relative to the root of the code
%  repository (D-20; N4G-D17), so the function runs from any current folder
% =========================================================================
FILE_PREFIX = 'new4_';                      % W1: every file written begins with this and carries the stamp
CSV_M1M3    = fullfile('results', 'regime_probe_rhi020', 'mc_snr_sweep_20260528_204613.csv');                % SPEC c.1, CR [R01]: the committed CSV of M1 and M3
CSV_M8      = fullfile('results', 'riviello_snr', 'mc_snr_sweep_20260610_203626.csv');                       % SPEC c.5: the committed CSV of M8
CSV_TRIALS  = fullfile('results', 'paperC_A3_20260907', 'mc_snr_trials_20260908_145855.csv');                % N4CB-D1: the A3 trials of the CTL block (not committed in the code repository)
CSV_A4SUM   = fullfile('results', 'paperC_A4_t50', 't50_a4_summary_20260919_200409.csv');                    % N4CB-D1, N4G-D12: A4's summary (not committed in the code repository)
A4_LEG       = 'snr';                       % N4-G6: the value of the column leg of A4's summary rows for the SNR sweep (A4 blob, writer lines of the summary)
A4_SOURCE_TS = '20260908_145855';           % N4-G6: the CTL block, the value of the column source_ts (A4 blob CTL_TS)
N_SNR_BLOCK = 13;                           % W5: SNR rows per trial in the trials file (A4 blob S_READ); checked against the M1/M3 file; OPEN-19
FF_SNR_TEXT = '10';                         % W4: the fail-fast row, SNR 10 dB, by the text of its SNR_dB cell
QUANT_SNR_TEXT = '10';                      % W11, N4G-D19: the SNR of the scene quantiles, by the text of its SNR_dB cell
QLEVELS     = [50 95];                      % W11: quantile levels in percent; OPEN-4
TOP_PCT     = 5;                            % W11, N-34: the top group is this percent of the trials; OPEN-6
COND_MAX    = 4.5e9;   % COPIED from kb/PaperC_B4_T50_T57_BoundSpec.md L495 (G9: max cond_eq threshold; W8)
LABEL_FA_OBS  = 'full-array CRB (observation x)';   % COPIED from kb/PaperC_Fact_Ledger.md L1728 (N-32 (b), confirmed 2026-09-23)
LABEL_CM_A    = 'compressed-domain CRB';   % COPIED from kb/PaperC_Fact_Ledger.md L1728 (N-32 (b), first part)
LABEL_CM_B    = '(observation y)';   % COPIED from kb/PaperC_Fact_Ledger.md L1729 (N-32 (b), second part)
LABEL_TM_A    = 'the square root of the';   % COPIED from kb/PaperC_Fact_Ledger.md L1724 (N-32 (a), 2026-09-23 row, first part)
LABEL_TM_B    = 'mean, over the plotted trials, of the per-trial CRB at each trial''s own';   % COPIED from kb/PaperC_Fact_Ledger.md L1725 (N-32 (a), second part)
LABEL_TM_C    = 'geometry and combiner';   % COPIED from kb/PaperC_Fact_Ledger.md L1726 (N-32 (a), third part)
LABEL_NB_A    = 'not a bound for WB-CL-KL as';   % COPIED from kb/PaperC_Fact_Ledger.md L1110 (N-32 (b), the note of the compressed curve, first part)
LABEL_NB_B    = 'implemented';   % COPIED from kb/PaperC_Fact_Ledger.md L1111 (N-32 (b), the note of the compressed curve, second part)
LABEL_SG      = 'single-geometry';   % COPIED from kb/PaperC_Fact_Ledger.md L1112 (N-32 (c), the word for a one-geometry curve)
LABEL_CM_OBS  = [LABEL_CM_A ' ' LABEL_CM_B];                          % N-32 (b) label of the compressed curve
LABEL_TM_AGG  = [LABEL_TM_A ' ' LABEL_TM_B ' ' LABEL_TM_C];           % N-32 (a) label of the trial-matched aggregation; OPEN-9
LABEL_NOTBOUND = [LABEL_NB_A ' ' LABEL_NB_B];                         % N-32 (b) note of a compressed curve drawn near WB-CL-KL
GEOM_TRIAL  = 'per-trial geometry (T-50 form) of the CTL block';      % N4G-D13: geometry named in the table; OPEN-10
GEOM_SINGLE = 'nominal geometry of loc_nominal (mid-angle and mid-range of the box)';   % N4G-D6, N-32 (c): the geometry of M8's curve; OPEN-10
COL_R_BOUND  = 'range_bound.sqrt_of_mean_over_trials_of_per_trial_variance_bound_full_Jacobian_equilibrated_exact_inverse_info_raw_eq.m';   % O1; N4-I5; OPEN-28: the family is a value of the column family, not part of this name
COL_TH_BOUND = 'angle_bound.sqrt_of_mean_over_trials_of_per_trial_variance_bound_theta_form_equilibrated_exact_inverse_info_raw_eq.deg';     % O1; N4-I5
COL_R_SINGLE = 'range_bound.sqrt_of_single_geometry_variance_bound_full_Jacobian_equilibrated_exact_inverse_info_raw_eq.m';                % O1 of M8; N4-I5
COL_COND_MAX = 'cond_Jn_equilibrated.max_over_trials.loc_inv_equilibrated_matlab_cond_info_raw_cond_eq_all.dimless';                       % O1; W8, C-12, C-13
COL_COND_ONE = 'cond_Jn_equilibrated.single_call.loc_inv_equilibrated_matlab_cond_info_raw_cond_eq_all.dimless';                           % O1 of M8; W8, C-12, C-13

here_dir = fileparts(mfilename('fullpath'));    % W1: the file stands in src/ of the code repository
root_dir = fileparts(here_dir);
f_m1m3   = fullfile(root_dir, CSV_M1M3);
f_m8     = fullfile(root_dir, CSV_M8);
f_trials = fullfile(root_dir, CSV_TRIALS);
f_a4sum  = fullfile(root_dir, CSV_A4SUM);

if is_full
    if exist(results_dir, 'dir') ~= 7
        fprintf('N4-FOLDER results folder [domain: results_dir before the run]: FAIL (the folder does not exist)\n');
        error('new4_bounds:folder', 'results_dir does not exist.');
    end
    d_old = dir(fullfile(results_dir, [FILE_PREFIX '*']));
    if ~isempty(d_old)
        fprintf('N4-FOLDER results folder [domain: results_dir before the run]: FAIL (%d file(s) begin with %s)\n', numel(d_old), FILE_PREFIX);
        error('new4_bounds:folder', 'results_dir already holds files that begin with %s.', FILE_PREFIX);
    end
end

fprintf('=== new4_bounds: stage %s, start %s ===\n', stage, datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('MATLAB release: %s\n', version);
fprintf('results_dir: %s\n', results_dir);
fprintf('stamp: %s\n', stamp);
fprintf('root of the code repository (from the path of this file): %s\n', root_dir);
fprintf('r_hi_fac_m1m3 = %.17g, r_hi_fac_m8 = %.17g (given by the run sheet with their source, F-130; N4G-D8 (iii))\n', ...
    r_hi_fac_m1m3, r_hi_fac_m8);

in_files = {f_m1m3, f_m8};
if is_full
    in_files = [in_files, {f_trials, f_a4sum}];
end
for i_f = 1:numel(in_files)
    if exist(in_files{i_f}, 'file') ~= 2
        fprintf('N4-INPUT: FAIL (input file not found: %s)\n', in_files{i_f});
        error('new4_bounds:input', 'input file not found: %s', in_files{i_f});
    end
    fprintf('input file: %s\n', in_files{i_f});
end

% =========================================================================
%  W3; N4-I2; N4G-D8 (iii), SPEC c.1, c.5: the config of the two curves
% =========================================================================
P  = loc_make_P(r_hi_fac_m1m3);     % M1 and M3: lines L37-38 of A4's blob, then r_hi_fac from the run sheet
P8 = loc_make_P(r_hi_fac_m8);       % M8: the same, into a separate struct
loc_print_P('N4-I2 M1/M3 struct', P);
loc_print_P('N4-I2 M8 struct', P8);
n_cfg_diff = loc_cfg_compare('M1/M3', f_m1m3, P) + loc_cfg_compare('M8', f_m8, P8);   % OPEN-2
fprintf('N4-I2 config columns of the committed CSVs against the structs [not a gate: N4-G2 below confirms the config, N4G-D8 (iii)]: columns that differ %d\n', ...
    n_cfg_diff);

[snr_s, pth_s, pr_s]  = loc_pub_crb_strings(f_m1m3);     % W17: bound columns and SNR_dB of the WB-CL-KL rows only
[snr8_s, ~, pr8_s]    = loc_pub_crb_strings(f_m8);
fprintf('N4-I3 published rows read: M1/M3 file %d WB-CL-KL rows, M8 file %d WB-CL-KL rows\n', numel(snr_s), numel(snr8_s));
S = numel(snr_s);                                         % OPEN-19: the SNR points of the plotted block are the SNR_dB cells of the M1/M3 file
snr_val = zeros(S, 1);
for s = 1:S
    snr_val(s) = str2double(snr_s(s));                    % converted once from the cell's own text
end
if is_full && S ~= N_SNR_BLOCK
    fprintf('N4-INPUT: FAIL (the M1/M3 file has %d WB-CL-KL rows, the trials file is read with %d SNR rows per trial)\n', S, N_SNR_BLOCK);
    error('new4_bounds:input', 'SNR row count of the M1/M3 file differs from N_SNR_BLOCK.');
end

% =========================================================================
%  W4; N4-G2; N4-I3; N4G-D11, N4G-D12, N4G-D8 (iii): reproduction of the
%  published bound columns through the default call (BoundSpec G1a)
% =========================================================================
repro = cell(0, 9);
g2_cells = 0;
g2_diff  = 0;
if is_full
    jj1 = 1:S;
    jj8 = 1:numel(snr8_s);
    g2_domain = sprintf('every WB-CL-KL row: M1 crb_r_m %d rows, M3 crb_theta_deg %d rows, M8 crb_r_m %d rows', S, S, numel(snr8_s));
else
    jj1 = find(snr_s == FF_SNR_TEXT);
    jj8 = [];
    g2_domain = 'one row, SNR 10 dB of the M1/M3 file, columns crb_r_m and crb_theta_deg (fail-fast stage)';
    if numel(jj1) ~= 1
        fprintf('N4-G2 reproduction [domain: %s]: cells compared 0, cells that differ 0: FAIL (the file has %d rows with SNR text %s)\n', ...
            g2_domain, numel(jj1), FF_SNR_TEXT);
        error('new4_bounds:G2', 'N4-G2 FAIL: the fail-fast row is not unique.');
    end
end
for j = jj1
    snr_dB_s = str2double(snr_s(j));
    [crb_r_val, crb_theta_val] = loc_default_call(P, snr_dB_s);
    cc = {'M1', 'crb_r_m', crb_r_val, pr_s(j); 'M3', 'crb_theta_deg', crb_theta_val, pth_s(j)};
    for k = 1:2
        t_csv = char(cc{k, 4});
        t_new = sprintf('%.6g', cc{k, 3});
        eq_ok = strcmp(t_new, t_csv);
        g2_cells = g2_cells + 1;
        if ~eq_ok
            g2_diff = g2_diff + 1;
            fprintf('N4-G2 DIFF: file %s, row SNR %s, column %s: csv %s, recomputed %s\n', CSV_M1M3, char(snr_s(j)), cc{k, 2}, t_csv, t_new);
        end
        repro(end + 1, :) = {'N4-G2', cc{k, 1}, CSV_M1M3, ['SNR_dB=' char(snr_s(j))], cc{k, 2}, t_csv, t_csv, t_new, double(eq_ok)}; %#ok<AGROW>
    end
end
for j = jj8
    snr_dB_s = str2double(snr8_s(j));
    [crb_r_val, ~] = loc_default_call(P8, snr_dB_s);
    t_csv = char(pr8_s(j));
    t_new = sprintf('%.6g', crb_r_val);
    eq_ok = strcmp(t_new, t_csv);
    g2_cells = g2_cells + 1;
    if ~eq_ok
        g2_diff = g2_diff + 1;
        fprintf('N4-G2 DIFF: file %s, row SNR %s, column %s: csv %s, recomputed %s\n', CSV_M8, char(snr8_s(j)), 'crb_r_m', t_csv, t_new);
    end
    repro(end + 1, :) = {'N4-G2', 'M8', CSV_M8, ['SNR_dB=' char(snr8_s(j))], 'crb_r_m', t_csv, t_csv, t_new, double(eq_ok)}; %#ok<AGROW>
end
if g2_diff == 0, g2_verdict = 'PASS'; else, g2_verdict = 'FAIL'; end
% OPEN-16: the exact text of every verdict line; OPEN-3: resolution and halting of N4-G2 and N4-G6
fprintf('N4-G2 reproduction [domain: %s; default call without opts as the driver L224-229; compared as text at %%.6g]: cells compared %d, cells that differ %d: %s\n', ...
    g2_domain, g2_cells, g2_diff, g2_verdict);
if g2_diff ~= 0
    fprintf('NEW4 HALT: N4-G2 FAIL. Run no further stage.\n');
    error('new4_bounds:G2', 'N4-G2 FAIL: %d cell(s) differ.', g2_diff);
end

if ~is_full
    fprintf('NEW4 FAILFAST COMPLETE: N4-G2 PASS at the one row; no file written; N4-G3 and the later stages were not run.\n');
    return
end

% =========================================================================
%  W5; N4-G3; N4-I3; N4G-D11, BoundSpec G0a, F-134 (i): the geometry of every
%  trial reconstructed from its seed and compared with the trials file
% =========================================================================
[seeds, base, th_str, r_str, g0b] = loc_read_trials(f_trials, N_SNR_BLOCK);
n = numel(seeds);
recon = false(n, 1);
th_i  = nan(n, 1);
r_i   = nan(n, 1);
Wc    = cell(n, 1);
for i = 1:n
        [th, r, Wi] = loc_recon(seeds(i), P);
        recon(i) = strcmp(sprintf('%.10g', th*180/pi), th_str(i)) && strcmp(sprintf('%.10g', r), r_str(i));
        th_i(i) = th; r_i(i) = r;
    Wc{i} = Wi;                                   % the trial's combiner, kept for W6
end
g3_diff = sum(~recon);
for i = 1:n
    t_th = sprintf('%.10g', th_i(i) * 180 / pi);
    t_r  = sprintf('%.10g', r_i(i));
    repro(end + 1, :) = {'N4-G3', 'M1/M3', CSV_TRIALS, sprintf('mc_seed=%d', seeds(i)), 'theta_true_deg', char(th_str(i)), char(th_str(i)), t_th, double(strcmp(t_th, th_str(i)))}; %#ok<AGROW>
    repro(end + 1, :) = {'N4-G3', 'M1/M3', CSV_TRIALS, sprintf('mc_seed=%d', seeds(i)), 'r_true_m', char(r_str(i)), char(r_str(i)), t_r, double(strcmp(t_r, r_str(i)))}; %#ok<AGROW>
end
if g3_diff == 0 && g0b, g3_verdict = 'PASS'; else, g3_verdict = 'FAIL'; end
fprintf('N4-G3 geometry reconstruction [domain: every trial of %s, mc_seed_base %d; theta in deg and r in m compared as text at %%.10g; reader flag (geometry equal at every SNR, trials in order) %d]: trials compared %d, trials that differ %d: %s\n', ...
    CSV_TRIALS, base, double(g0b), n, g3_diff, g3_verdict);
if ~strcmp(g3_verdict, 'PASS')
    fprintf('NEW4 HALT: N4-G3 FAIL. Run no further stage.\n');
    error('new4_bounds:G3', 'N4-G3 FAIL: %d trial(s) differ, reader flag %d.', g3_diff, double(g0b));
end

% =========================================================================
%  W6, W7, W8; N4-I4, N4-G4, N4-G5; N4G-D1, D3, D9, D11, D12: two calls per
%  trial and SNR, family full array (W = I) and family trial-matched. All SNR
%  points of the block are computed (OPEN-18); no parallel pool is used (OPEN-22).
% =========================================================================
kinds = {'call_error', 'nonfinite', 'nonpositive'};
vr_fa  = nan(n, S);  vth_fa = nan(n, S);  cd_fa = nan(n, S);
vr_tm  = nan(n, S);  vth_tm = nan(n, S);  cd_tm = nan(n, S);
n_checked = 0;
for s = 1:S
    N0 = 10^(-snr_val(s) / 10);                   % p = 1 and N0 = 10^(-SNR/10), as the driver L220-221 and A4's loc_vars
    for i = 1:n
        n_checked = n_checked + 1;
        R = loc_checked_call(th_i(i), r_i(i), N0, P, {eye(P.M)});
        if ~isempty(R.kind)
            loc_g4_fail(kinds, R.kind, n_checked, 'per-trial bounds', ...
                sprintf('SNR %s, trial %d, seed %d, family full array', char(snr_s(s)), i, seeds(i)), R.msg);
        end
        vr_fa(i, s) = R.vr;  vth_fa(i, s) = R.vth;  cd_fa(i, s) = R.cnd;
        n_checked = n_checked + 1;
        R = loc_checked_call(th_i(i), r_i(i), N0, P, {Wc{i}});
        if ~isempty(R.kind)
            loc_g4_fail(kinds, R.kind, n_checked, 'per-trial bounds', ...
                sprintf('SNR %s, trial %d, seed %d, family trial-matched', char(snr_s(s)), i, seeds(i)), R.msg);
        end
        vr_tm(i, s) = R.vr;  vth_tm(i, s) = R.vth;  cd_tm(i, s) = R.cnd;
    end
    fprintf('N4-I4 SNR %s dB (%d of %d): %d trials x 2 families computed\n', char(snr_s(s)), s, S, n);
end
fprintf('N4-G4 validity [domain: every call of the per-trial bounds, %d trials x %d SNR points x 2 families; finite and positive range and angle variance bounds, finite cond]: calls checked %d, calls invalid 0 (by kind: %s): PASS\n', ...
    n, S, n_checked, loc_kind_zero(kinds));

maxc_fa = max(cd_fa, [], 1);                      % 1 x S, max over the trials
maxc_tm = max(cd_tm, [], 1);
n_over  = sum(maxc_fa > COND_MAX) + sum(maxc_tm > COND_MAX);   % OPEN-21: the number itself is compared
if n_over == 0, g5_verdict = 'PASS'; else, g5_verdict = 'FAIL'; end
fprintf('N4-G5 reliability of the equilibrated inverse [domain: per SNR point and family over the %d trials; quantity: cond of the equilibrated matrix Jn of loc_inv_equilibrated, helper MATLAB cond, from info.raw.cond_eq_all (rulings C-12, C-13); not the retained-subspace number of loc_pinv_svd]: largest value over points and families %s, threshold %s (BoundSpec G9; N4G-D11), points and families over the threshold %d of %d: %s\n', ...
    n, sprintf('%.6g', max([maxc_fa, maxc_tm])), sprintf('%.6g', COND_MAX), n_over, 2 * S, g5_verdict);
if n_over ~= 0
    fprintf('NEW4 HALT: N4-G5 FAIL. Run no further stage.\n');
    error('new4_bounds:G5', 'N4-G5 FAIL: %d point(s) over the threshold.', n_over);
end

% =========================================================================
%  W9; N4-I5; N4G-D1, D3, D9: aggregation in the variance domain (T-50 form)
% =========================================================================
agg_r_fa  = sqrt(mean(vr_fa,  1));                % 1 x S, m
agg_th_fa = sqrt(mean(vth_fa, 1));                % 1 x S, deg
agg_r_tm  = sqrt(mean(vr_tm,  1));
agg_th_tm = sqrt(mean(vth_tm, 1));

% =========================================================================
%  W10; N4-G6; N4-I3; N4G-D12: A4's aggregates at the CTL block, recomputed
%  here and compared with A4's summary rows (no USW-model leg)
% =========================================================================
A4 = loc_read_a4_summary(f_a4sum);
fams  = {'fa', 'cm'};                             % A4's family tags: full array, trial-matched
forms = {'r_fj', 'theta'};                        % A4's form tags: full-Jacobian range, angle
agg_all = {agg_r_fa, agg_th_fa; agg_r_tm, agg_th_tm};
g6_cells = 0;
g6_diff  = 0;
for s = 1:S
    for fi = 1:2
        for fo = 1:2
            t_new = sprintf('%.6g', agg_all{fi, fo}(s));
            key = sprintf('SNR_dB=%s;leg=%s;source_ts=%s;family=%s;inverse=eq;form=%s;agg=T50', char(snr_s(s)), A4_LEG, A4_SOURCE_TS, fams{fi}, forms{fo});
            kk = find(A4.leg == A4_LEG & A4.source_ts == A4_SOURCE_TS & A4.family == fams{fi} & A4.inverse == "eq" ...
                & A4.form == forms{fo} & A4.agg == "T50" & abs(str2double(A4.SNR_dB) - snr_val(s)) < 1e-9);   % OPEN-24
            g6_cells = g6_cells + 1;
            if numel(kk) ~= 1
                g6_diff = g6_diff + 1;
                fprintf('N4-G6 DIFF: %s: the summary file holds %d matching rows\n', key, numel(kk));
                repro(end + 1, :) = {'N4-G6', 'M1/M3', CSV_A4SUM, key, 'value', sprintf('%d matching rows', numel(kk)), '', t_new, 0}; %#ok<AGROW>
            else
                t_written = char(A4.value(kk));
                t_csv6 = sprintf('%.6g', str2double(t_written));
                eq_ok = strcmp(t_csv6, t_new);
                if ~eq_ok
                    g6_diff = g6_diff + 1;
                    fprintf('N4-G6 DIFF: %s: A4 %s, recomputed %s\n', key, t_csv6, t_new);
                end
                repro(end + 1, :) = {'N4-G6', 'M1/M3', CSV_A4SUM, key, 'value', t_written, t_csv6, t_new, double(eq_ok)}; %#ok<AGROW>
            end
        end
    end
end
if g6_diff == 0, g6_verdict = 'PASS'; else, g6_verdict = 'FAIL'; end
fprintf('N4-G6 A4 reproduction [domain: the %d SNR points of the CTL block, families full array and trial-matched, range (full-Jacobian) and angle, T50 aggregate, equilibrated inverse; compared as text at %%.6g; no USW-model leg]: cells compared %d, cells that differ %d: %s\n', ...
    S, g6_cells, g6_diff, g6_verdict);
if g6_diff ~= 0
    fprintf('NEW4 HALT: N4-G6 FAIL. Run no further stage.\n');
    error('new4_bounds:G6', 'N4-G6 FAIL: %d cell(s) differ.', g6_diff);
end

% =========================================================================
%  W11; N4G-D19; N-34; SPEC O3: scene quantiles of sqrt(CRB(g_i, W_i)) at one
%  stated SNR, and the top-group share of the mean variance (no band, no interval)
% =========================================================================
sq = find(snr_s == QUANT_SNR_TEXT);
if numel(sq) ~= 1
    fprintf('N-34 quantiles [domain: SNR text %s]: FAIL (the M1/M3 file has %d rows with that SNR text)\n', QUANT_SNR_TEXT, numel(sq));
    error('new4_bounds:quantiles', 'N-34 quantile SNR is not unique.');
end
v_q = vr_tm(:, sq);                               % per-trial range variances, trial-matched family, m^2; OPEN-7
x_q = sqrt(v_q);                                  % sqrt(CRB(g_i, W_i)), m
qrows = cell(0, 8);
for iq = 1:numel(QLEVELS)
    qrows(end + 1, :) = {char(snr_s(sq)), snr_val(sq), 'trial_matched', n, ...
        'quantile_nearest_rank_of_sqrt_of_per_trial_range_variance_bound', QLEVELS(iq), loc_quantile_nr(x_q, QLEVELS(iq)), 'm'}; %#ok<AGROW>
end
vs = sort(v_q, 'descend');
k_top = ceil(n * TOP_PCT / 100);                  % OPEN-6: the group size is rounded up; integer arithmetic
qrows(end + 1, :) = {char(snr_s(sq)), snr_val(sq), 'trial_matched', n, 'top_group_size_in_trials', TOP_PCT, k_top, 'count'};
qrows(end + 1, :) = {char(snr_s(sq)), snr_val(sq), 'trial_matched', n, ...
    'share_of_the_mean_range_variance_bound_in_the_top_group_of_trials', TOP_PCT, sum(vs(1:k_top)) / sum(vs), 'dimensionless'};
fprintf('N-34 scene quantiles [domain: SNR text %s, %d trials, family trial-matched; levels in percent %s; top group %d percent = %d trials; nearest-rank quantiles]: computed, written to the quantile table\n', ...
    QUANT_SNR_TEXT, n, mat2str(QLEVELS), TOP_PCT, k_top);

% =========================================================================
%  W12; N4G-D6, SPEC c.5; N-32 (c): the bound of M8 at one stated geometry
% =========================================================================
[th8, r8] = loc_nominal(P8);
S8 = numel(snr8_s);
vr8 = nan(S8, 1);
cd8 = nan(S8, 1);
snr8_val = zeros(S8, 1);
n_checked8 = 0;
for j = 1:S8
    snr8_val(j) = str2double(snr8_s(j));
    N0 = 10^(-snr8_val(j) / 10);
    n_checked8 = n_checked8 + 1;
    R = loc_checked_call(th8, r8, N0, P8, {eye(P8.M)});
    if ~isempty(R.kind)
        loc_g4_fail(kinds, R.kind, n_checked8, 'M8 single geometry', sprintf('SNR %s', char(snr8_s(j))), R.msg);
    end
    vr8(j) = R.vr;
    cd8(j) = R.cnd;
end
agg_r8 = sqrt(vr8);                               % m, one geometry, one call: no aggregation
fprintf('N4-G4 validity [domain: every call of the M8 single-geometry bound, %d SNR points; finite and positive variance bounds, finite cond]: calls checked %d, calls invalid 0 (by kind: %s): PASS\n', ...
    S8, n_checked8, loc_kind_zero(kinds));
n_over8 = sum(cd8 > COND_MAX);
if n_over8 == 0, g5_verdict8 = 'PASS'; else, g5_verdict8 = 'FAIL'; end
fprintf('N4-G5 reliability of the equilibrated inverse [domain: the %d SNR points of M8; quantity: cond of the equilibrated matrix Jn of loc_inv_equilibrated, helper MATLAB cond, from info.raw.cond_eq_all (rulings C-12, C-13)]: largest value %s, threshold %s, points over the threshold %d of %d: %s\n', ...
    S8, sprintf('%.6g', max(cd8)), sprintf('%.6g', COND_MAX), n_over8, S8, g5_verdict8);
if n_over8 ~= 0
    fprintf('NEW4 HALT: N4-G5 FAIL (M8). Run no further stage.\n');
    error('new4_bounds:G5', 'N4-G5 FAIL (M8): %d point(s) over the threshold.', n_over8);
end

% =========================================================================
%  W13; N4-I7; N4G-D17: files (O1 to O3 and the quantile table); one header
%  line naming every column; values at %.17g; fprintf, not writetable. The files
%  are written here, after every gate; a halt before this point writes none (OPEN-17).
% =========================================================================
f_out = @(kind) fullfile(results_dir, [FILE_PREFIX kind '_' stamp '.csv']);
f1 = f_out('bound_m1');           % OPEN-8: one table per regenerated curve
f3 = f_out('bound_m3');
f8 = f_out('bound_m8');
f2 = f_out('repro');
fT = f_out('pertrial');
fQ = f_out('quantiles');

hdr_cfg = {'M', 'N_RF', 'N', 'K_s', 'bandwidth_Hz', 'd', 'r_RD_m', 'r_lo_fac', 'r_hi_fac', 'theta_lo_deg', 'theta_hi_deg'};
cfg1 = loc_cfg_cols(P);
cfg8 = loc_cfg_cols(P8);
hdr_lab = {'map_row', 'curve_id', 'family', 'geometry', 'snr_dB_text', 'snr_dB', ...
    'label_observation_N32b', 'label_aggregation_N32a', 'label_not_a_bound_N32b'};      % OPEN-9, OPEN-10

hdr_m1 = [hdr_lab, hdr_cfg, {'n_trials', COL_R_BOUND, COL_COND_MAX}];
hdr_m3 = [hdr_lab, hdr_cfg, {'n_trials', COL_TH_BOUND, COL_COND_MAX}];
rows_m1 = cell(0, numel(hdr_m1));
rows_m3 = cell(0, numel(hdr_m3));
for s = 1:S
    pre_fa = {'full_array', GEOM_TRIAL, char(snr_s(s)), snr_val(s), LABEL_FA_OBS, '', ''};
    pre_tm = {'trial_matched', GEOM_TRIAL, char(snr_s(s)), snr_val(s), LABEL_CM_OBS, LABEL_TM_AGG, LABEL_NOTBOUND};
    rows_m1(end + 1, :) = [{'M1', 'M1_full_array'},    pre_fa, cfg1, {n, agg_r_fa(s),  maxc_fa(s)}]; %#ok<AGROW>
    rows_m1(end + 1, :) = [{'M1', 'M1_trial_matched'}, pre_tm, cfg1, {n, agg_r_tm(s),  maxc_tm(s)}]; %#ok<AGROW>
    rows_m3(end + 1, :) = [{'M3', 'M3_full_array'},    pre_fa, cfg1, {n, agg_th_fa(s), maxc_fa(s)}]; %#ok<AGROW>
    rows_m3(end + 1, :) = [{'M3', 'M3_trial_matched'}, pre_tm, cfg1, {n, agg_th_tm(s), maxc_tm(s)}]; %#ok<AGROW>
end

hdr_m8 = [{'map_row', 'curve_id', 'family', 'geometry', 'snr_dB_text', 'snr_dB', 'label_observation_N32b', ...
    'label_single_geometry_N32c', 'theta_nom_deg', 'r_nom_m'}, hdr_cfg, {'n_geometries', COL_R_SINGLE, COL_COND_ONE}];
rows_m8 = cell(0, numel(hdr_m8));
for j = 1:S8
    rows_m8(end + 1, :) = [{'M8', 'M8_full_array_single_geometry', 'full_array', GEOM_SINGLE, char(snr8_s(j)), snr8_val(j), ...
        LABEL_FA_OBS, LABEL_SG, th8 * 180 / pi, r8}, cfg8, {1, agg_r8(j), cd8(j)}]; %#ok<AGROW>
end

hdr_rep = {'gate', 'map_row', 'file', 'row_key', 'column', 'csv_text_as_written', ...
    'csv_text_at_comparison_resolution', 'recomputed_text_at_comparison_resolution', 'equal_1yes_0no'};

hdr_pt = {'trial_index', 'mc_seed_base', 'mc_seed', 'theta_deg_reconstructed', 'r_m_reconstructed', 'snr_dB_text', 'snr_dB', ...
    'full_array.range_variance_bound_full_Jacobian_equilibrated_info_raw_eq.m2', ...
    'full_array.angle_variance_bound_theta_form_equilibrated_info_raw_eq.deg2', ...
    'full_array.cond_Jn_equilibrated.loc_inv_equilibrated_matlab_cond_info_raw_cond_eq_all.dimless', ...
    'trial_matched.range_variance_bound_full_Jacobian_equilibrated_info_raw_eq.m2', ...
    'trial_matched.angle_variance_bound_theta_form_equilibrated_info_raw_eq.deg2', ...
    'trial_matched.cond_Jn_equilibrated.loc_inv_equilibrated_matlab_cond_info_raw_cond_eq_all.dimless'};
rows_pt = cell(n * S, numel(hdr_pt));
ir = 0;
for s = 1:S
    for i = 1:n
        ir = ir + 1;
        rows_pt(ir, :) = {i, base, seeds(i), th_i(i) * 180 / pi, r_i(i), char(snr_s(s)), snr_val(s), ...
            vr_fa(i, s), vth_fa(i, s), cd_fa(i, s), vr_tm(i, s), vth_tm(i, s), cd_tm(i, s)};
    end
end

hdr_q = {'snr_dB_text', 'snr_dB', 'family', 'n_trials', 'statistic', 'parameter', 'value', 'unit'};

loc_write_csv(f1, hdr_m1, rows_m1);
loc_write_csv(f3, hdr_m3, rows_m3);
loc_write_csv(f8, hdr_m8, rows_m8);
loc_write_csv(f2, hdr_rep, repro);
loc_write_csv(fT, hdr_pt, rows_pt);
loc_write_csv(fQ, hdr_q, qrows);

% =========================================================================
%  W14; N4-I8: completion line naming every file written
% =========================================================================
fprintf('NEW4 COMPLETE: stage full, %d trials x %d SNR points, M8 %d SNR points; files written: %s ; %s ; %s ; %s ; %s ; %s\n', ...
    n, S, S8, f1, f3, f8, f2, fT, fQ);
end


% =========================================================================
%  LOCAL FUNCTIONS COPIED LINE FOR LINE FROM THE BLOB OF diagnostics/t50_a4_bound.m
%  (W2; N4CB-G1 item 2). loc_recon is the only place of this file that sets
%  the generator state (OPEN-23).
% =========================================================================
% COPIED from diagnostics/t50_a4_bound.m@0ffe3d3a458d L347-351 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
function [th0, r0] = loc_nominal(Q)
% run_monte_carlo_paperC.m L217-219 (live)
th0 = (Q.theta_lo + Q.theta_hi) / 2;
r0  = Q.r_lo_fac * Q.r_RD + (Q.r_hi_fac - Q.r_lo_fac) * Q.r_RD / 2;
end

% COPIED from diagnostics/t50_a4_bound.m@0ffe3d3a458d L363-375 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
function [th, r, W] = loc_recon(seed, Q)
% Trial geometry and combiner from mc_seed, in the generator's draw order:
% run_monte_carlo_paperC.m L429 -> wb_channel_gen_ofdm_nf.m L132-135, L151-152, L159-160
rng(seed, 'twister');
D     = (Q.M - 1) * Q.d_ant;
r_RD  = 2 * D^2 / Q.lambda_c;
r_min = Q.r_lo_fac * r_RD;
r_max = Q.r_hi_fac * r_RD;
th    = Q.theta_lo + (Q.theta_hi - Q.theta_lo) * rand(Q.d, 1);
r     = r_min      + (r_max - r_min)            * rand(Q.d, 1);
phi   = 2*pi * rand(Q.M, Q.N_RF);
W     = (1/sqrt(Q.M)) * exp(1j * phi);
end

% COPIED from diagnostics/t50_a4_bound.m@0ffe3d3a458d L459-471 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
function [snr_s, pth_s, pr_s] = loc_pub_crb_strings(csv_path)
L = splitlines(string(fileread(csv_path))); L = L(strlength(L) > 0);
hdr = split(L(1), ',');
im = find(hdr == "method"); is = find(hdr == "SNR_dB");
it = find(hdr == "crb_theta_deg"); ir = find(hdr == "crb_r_m");
snr_s = strings(0,1); pth_s = strings(0,1); pr_s = strings(0,1);
for j = 2:numel(L)
    f = split(L(j), ',');
    if f(im) == "WB-CL-KL"
        snr_s(end+1,1) = f(is); pth_s(end+1,1) = f(it); pr_s(end+1,1) = f(ir); %#ok<AGROW>
    end
end
end

% COPIED from diagnostics/t50_a4_bound.m@0ffe3d3a458d L473-485 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
function [seeds, base, th_str, r_str, ok] = loc_read_trials(csv_path, S)
L = splitlines(string(fileread(csv_path))); L = L(strlength(L) > 0);
hdr = split(L(1), ','); C = split(L(2:end), ',');
ic = @(nm) find(hdr == nm);
n = size(C, 1) / S;
idx  = str2double(C(:, ic("mc_idx")));
thS  = C(:, ic("theta_true_deg")); rS = C(:, ic("r_true_m"));
ok = mod(size(C,1), S) == 0 && isequal(idx, repmat((1:n).', S, 1)) && ...
     isequal(thS, repmat(thS(1:n), S, 1)) && isequal(rS, repmat(rS(1:n), S, 1));
seeds  = str2double(C(1:n, ic("mc_seed")));
base   = str2double(C(1, ic("mc_seed_base")));
th_str = thS(1:n); r_str = rS(1:n);
end


% =========================================================================
%  LOCAL FUNCTIONS WRITTEN FOR THIS RUN
% =========================================================================
function P = loc_make_P(r_hi_fac)
% W3; N4-I2; N4G-D8 (iii): the effective struct of a curve. Lines L37-38 of A4's
% blob as they stand, then r_hi_fac from the run sheet. The blob's L39 sets the
% same field from a constant of A4's own run; the blob's L40-41 are not copied
% (OPEN-20).
% COPIED from diagnostics/t50_a4_bound.m@0ffe3d3a458d L37-38 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
P = setup_production_P_v4('snr', 'full');
if ~isfield(P, 'r_RD'), P.r_RD = 2*((P.M-1)*P.d_ant)^2 / P.lambda_c; end
P.r_hi_fac = r_hi_fac;
end

function [crb_r_val, crb_theta_val] = loc_default_call(P_s, SNR_dB_s)
% W4; N4-G2: the default call exactly as the driver forms it, no opts: the
% nominal geometry, equal powers, N0 from the SNR, the call, the mean over the
% draws of the returned values.
% COPIED from src/run_monte_carlo_paperC.m@ab795eb3bfa4 L217-221 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
    theta_nom = (P_s.theta_lo + P_s.theta_hi) / 2;
    r_nom     = P_s.r_lo_fac * P_s.r_RD + ...
                (P_s.r_hi_fac - P_s.r_lo_fac) * P_s.r_RD / 2;
    p_nom     = ones(P_s.d, 1) / P_s.d;
    N0_nom    = 10^(-SNR_dB_s / 10);
% COPIED from src/run_monte_carlo_paperC.m@ab795eb3bfa4 L224-229 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
        [crb_r_s, crb_theta_s, ~] = wb_crb_compressed( ...
            repmat(theta_nom, P_s.d, 1), ...
            repmat(r_nom,     P_s.d, 1), ...
            p_nom, N0_nom, P_s);
        crb_r_val     = mean(crb_r_s);       % scalar [m],  sqrt-CRB
        crb_theta_val = mean(crb_theta_s);   % scalar [deg], sqrt-CRB
end

function [vr, vth, cnd] = loc_bound_call(theta, r, N0, P, Wl)
% W6; N4-I4; N4G-D1: one call through the D1 interface at the trial's geometry
% (theta in rad, r in m), p = 1, the combiner list Wl; the variances are formed
% from info.raw.eq (the equilibrated exact inverse) by the mapping of BoundSpec
% L313-318, as A4's loc_vars forms them. Nothing is taken from info.raw.pinv.
o = struct(); o.W_list = Wl; o.return_raw = true;
[~, ~, info] = wb_crb_compressed(theta, r, 1, N0, P, o);
% COPIED from diagnostics/t50_a4_bound.m@0ffe3d3a458d L65-66 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
c1 = 2*pi*P.d_ant/P.lambda_c;                 % omega = c1*cos(theta)
c0 = pi*P.d_ant^2/P.lambda_c;                 % kappa = c0*sin^2(theta)/r
% COPIED from diagnostics/t50_a4_bound.m@0ffe3d3a458d L382-383 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
kap  = c0 * sin(theta)^2 / r;
g    = [-2*r*cos(theta) / (c1*sin(theta)^2); -r/kap];   % d r / d(omega, kappa)
% COPIED from diagnostics/t50_a4_bound.m@0ffe3d3a458d L385 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
dwdt = c1 * sin(theta);                                  % |d omega / d theta|
Coo = info.raw.eq.C_oo(1, 1);
Ckk = info.raw.eq.C_kk(1, 1);
Cok = info.raw.eq.C_ok(1, 1);
% COPIED from diagnostics/t50_a4_bound.m@0ffe3d3a458d L391-392 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
    vr_f = g(1)^2*Coo + 2*g(1)*g(2)*Cok + g(2)^2*Ckk;    % full-Jacobian
    vt   = max(Coo, 0) / dwdt^2 * (180/pi)^2;            % L202, L209 convention
vr  = vr_f;
vth = vt;
cnd = info.raw.cond_eq_all(1);
end

function R = loc_checked_call(theta, r, N0, P, Wl)
% W7; N4-G4: a call of loc_bound_call inside try/catch (the equilibrated inverse
% asserts a positive finite diagonal, CR [R102]); R.kind is empty for a valid
% call, else the kind of the first invalid value.
R.vr = NaN; R.vth = NaN; R.cnd = NaN; R.kind = ''; R.msg = '';
try
    [R.vr, R.vth, R.cnd] = loc_bound_call(theta, r, N0, P, Wl);
catch ME
    R.kind = 'call_error';
    R.msg  = ME.message;
    return
end
if ~(isreal(R.vr) && isfinite(R.vr)) || ~(isreal(R.vth) && isfinite(R.vth)) || ~(isreal(R.cnd) && isfinite(R.cnd))
    R.kind = 'nonfinite';
elseif R.vr <= 0 || R.vth <= 0
    R.kind = 'nonpositive';
end
end

function loc_g4_fail(kinds, kind, n_checked, domain, ctx, msg)
% W7; N4-G4: the first invalid value prints the FAIL line and halts
parts = cell(1, numel(kinds));
for i = 1:numel(kinds)
    parts{i} = sprintf('%s %d', kinds{i}, double(strcmp(kinds{i}, kind)));
end
fprintf('N4-G4 validity [domain: %s]: calls checked %d, calls invalid 1 (by kind: %s): FAIL\n', domain, n_checked, strjoin(parts, ', '));
fprintf('N4-G4 first invalid value: %s, kind %s, message: %s\n', ctx, kind, msg);
fprintf('NEW4 HALT: N4-G4 FAIL. Run no further stage; no file written.\n');
error('new4_bounds:G4', 'N4-G4 FAIL (%s, %s).', ctx, kind);
end

function s = loc_kind_zero(kinds)
parts = cell(1, numel(kinds));
for i = 1:numel(kinds)
    parts{i} = sprintf('%s 0', kinds{i});
end
s = strjoin(parts, ', ');
end

function loc_print_P(tag, P)
% W3; N4-I2: every scalar field at %.17g, class and size of every other field
fprintf('--- %s: effective parameter struct (every scalar field at %%.17g) ---\n', tag);
fn = fieldnames(P);
for i = 1:numel(fn)
    v = P.(fn{i});
    if (isnumeric(v) || islogical(v)) && isscalar(v)
        fprintf('P.%s = %.17g\n', fn{i}, double(v));
    else
        fprintf('P.%s : %s, size %s\n', fn{i}, class(v), mat2str(size(v)));
    end
end
fprintf('%s named: N_seed = %.17g, rng_seed_W = %.17g, r_hi_fac = %.17g\n', tag, P.N_seed, P.rng_seed_W, P.r_hi_fac);
fprintf('--- end of %s ---\n', tag);
end

function n_diff = loc_cfg_compare(tag, csv_path, P)
% W3; N4-I2: the config columns of a committed CSV (WB-CL-KL rows) printed beside
% the struct's values, each as text in the form the Monte Carlo driver wrote it
% (integers at %d, B_hz and r_RD_m at %.6g; the driver's fmt_num and L852-862).
% A disagreement is counted and printed; it does not halt (OPEN-2).
cols = {'M', 'N_RF', 'N', 'K_s', 'B_hz', 'd', 'r_RD_m'};
want = {sprintf('%d', P.M), sprintf('%d', P.N_RF), sprintf('%d', P.N), sprintf('%d', P.K_s), ...
        sprintf('%.6g', P.B), sprintf('%d', P.d), sprintf('%.6g', P.r_RD)};
L = splitlines(string(fileread(csv_path))); L = L(strlength(L) > 0);
hdr = split(L(1), ',');
im = find(hdr == "method");
n_diff = 0;
for c = 1:numel(cols)
    ic = find(hdr == string(cols{c}));
    if isempty(ic)
        n_diff = n_diff + 1;
        fprintf('N4-CFG [%s] column %s: the CSV has no such column: DIFFERENT\n', tag, cols{c});
        continue
    end
    seen = strings(0, 1);
    for j = 2:numel(L)
        f = split(L(j), ',');
        if f(im) == "WB-CL-KL"
            seen(end + 1, 1) = f(ic); %#ok<AGROW>
        end
    end
    u = unique(seen);
    same = (numel(u) == 1) && (u(1) == string(want{c}));
    if ~same, n_diff = n_diff + 1; end
    if same, word = 'equal'; else, word = 'DIFFERENT'; end
    fprintf('N4-CFG [%s] column %s: csv text(s) {%s}, struct text %s: %s\n', tag, cols{c}, strjoin(cellstr(u(:).'), ' | '), want{c}, word);
end
end

function c = loc_cfg_cols(P)
% O1: the config columns of a table row, in the order of hdr_cfg
c = {P.M, P.N_RF, P.N, P.K_s, P.B, P.d, P.r_RD, P.r_lo_fac, P.r_hi_fac, P.theta_lo * 180 / pi, P.theta_hi * 180 / pi};
end

function A = loc_read_a4_summary(csv_path)
% W10; N4-G6: A4's summary file read as text, columns found by their names
L = splitlines(string(fileread(csv_path))); L = L(strlength(L) > 0);
hdr = split(L(1), ',');
C = split(L(2:end), ',');
names = {'leg', 'source_ts', 'SNR_dB', 'family', 'inverse', 'form', 'agg', 'value'};
for k = 1:numel(names)
    ic = find(hdr == string(names{k}));
    if numel(ic) ~= 1
        error('new4_bounds:a4summary', 'column %s not found once in %s.', names{k}, csv_path);
    end
    A.(names{k}) = C(:, ic);
end
end

function x = loc_quantile_nr(v, q)
% W11: the nearest-rank quantile, level q in percent: the smallest sorted value
% whose rank is at least q percent of the count (OPEN-5)
vs = sort(v(:));
k = ceil(numel(vs) * q / 100);
k = max(k, 1);
x = vs(k);
end

% The two CSV helpers below are copied from E3's function of record (OPEN-25).
% COPIED from src/wb_crb_d2_e3.m@fc6d9c3502a5 L769-775 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
function s = loc_csvfield(t)
if any(t == ',') || any(t == '"')
    s = ['"' strrep(t, '"', '""') '"'];
else
    s = t;
end
end

% COPIED from src/wb_crb_d2_e3.m@fc6d9c3502a5 L777-797 of the code repository at 69a752cccbc5 (blob, CR removed) by copy_new4.pl
function loc_write_csv(fname, hdr, rows)
% S13: header line, then one line per row; numbers at %.17g; fprintf, not writetable
fid = fopen(fname, 'w');
if fid < 0
    error('wb_crb_d2_e3:file', 'cannot open %s for writing.', fname);
end
cl = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:size(rows, 1)
    parts = cell(1, size(rows, 2));
    for j = 1:size(rows, 2)
        v = rows{i, j};
        if ischar(v)
            parts{j} = loc_csvfield(v);
        else
            parts{j} = sprintf('%.17g', double(v));
        end
    end
    fprintf(fid, '%s\n', strjoin(parts, ','));
end
end
