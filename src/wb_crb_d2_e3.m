function wb_crb_d2_e3(stage, results_dir, stamp)
%WB_CRB_D2_E3  Step E3-1 and the compressed-domain bound at two paths (E3).
%
%  wb_crb_d2_e3('failfast', results_dir, stamp)
%  wb_crb_d2_e3('full',     results_dir, stamp)
%
%  Written by the C4 session of paperC-lab (sessions/2026-10-16_C4) from
%  sessions/2026-10-16_C4-SPEC/C4_SPEC.md as decided in DECISIONS_C4G.md and
%  DECISIONS_C4A.md. Held at a commit of the code repository only after the
%  author's review (C4A-D1). Written as text; never parsed or executed by
%  its author (N6).
%
%  The function is a bound computation. It calls no estimator and writes or
%  prints no error statistic, no estimator output and no ratio of any
%  estimator statistic to any bound.
%
%  Called functions: setup_production_P_v4, wb_crb_compressed (unchanged)
%  and MATLAB built-ins. No random-number call and no combiner is drawn
%  here: the draws come from the own mechanism of wb_crb_compressed
%  (C4G-D3 iii, C4G-D13).
%
%  stage 'failfast': the reproduction gate E3-G3 at one row (SNR 10 dB,
%                    N_RF 8) and nothing else; writes no file (C4G-D8 iv).
%  stage 'full'    : the reproduction gate at all rows, then Step E3-1 and
%                    the bounds at every point, then the files.
%
%  results_dir and stamp are fixed by the final run sheet and passed in
%  (C4G-D19). Every file written is in results_dir and is named
%  wb_crb_d2_e3_<kind>_<stamp>.csv.
%
%  Labels in comments: Sxx = requirement of the C4 starter; Ixx = step of
%  C4_SPEC.md Sec. g; C4G-Dn, C4A-Dn = author decisions; OPEN-n = a choice
%  no record fixes, listed in OPEN_C4.md.

% =========================================================================
%  I1; S1: arguments, stage, results folder, start line
% =========================================================================
if nargin ~= 3
    error('wb_crb_d2_e3:usage', ...
        'wb_crb_d2_e3(stage, results_dir, stamp): three arguments are required.');
end
if ~ischar(stage) || ~(strcmp(stage, 'failfast') || strcmp(stage, 'full'))
    error('wb_crb_d2_e3:stage', 'stage must be ''failfast'' or ''full''.');
end
if ~ischar(results_dir) || isempty(results_dir) || ~ischar(stamp) || isempty(stamp)
    error('wb_crb_d2_e3:args', 'results_dir and stamp must be non-empty character vectors.');
end
is_full = strcmp(stage, 'full');

if is_full
    if exist(results_dir, 'dir') ~= 7
        fprintf('E3-S1 results folder [domain: results_dir before the run]: FAIL (the folder does not exist)\n');
        error('wb_crb_d2_e3:folder', 'results_dir does not exist.');
    end
    d_old = dir(fullfile(results_dir, 'wb_crb_d2_e3_*'));       % OPEN-11
    if ~isempty(d_old)
        fprintf('E3-S1 results folder [domain: results_dir before the run]: FAIL (%d file(s) begin with wb_crb_d2_e3_)\n', numel(d_old));
        error('wb_crb_d2_e3:folder', 'results_dir already holds files that begin with wb_crb_d2_e3_.');
    end
end

fprintf('=== wb_crb_d2_e3: stage %s, start %s ===\n', stage, datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('MATLAB release: %s\n', version);
fprintf('results_dir: %s\n', results_dir);
fprintf('stamp: %s\n', stamp);

% The root of the code repository from this file's own path (S1, C4A-D1).
here_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(here_dir);
anchors_file = fullfile(root_dir, 'results', 'paperC_phase35_fig6_crb_d2_v2', ...
    'wb_crb_d2_sweep_20260526_225015_v3.csv');
anchors_header = 'snr_db,N_RF,crb_r1_m,crb_r2_m,crb_theta1_deg,crb_theta2_deg,crb_r1_std_m,crb_r2_std_m,cond_J_mean,n_singular_mean';   % COPIED from sessions/2026-10-16_C4-SPEC/C4_CODE_READING.md L58 (header line only)
fprintf('anchors file: %s\n', anchors_file);

% =========================================================================
%  Registered constants, each copied by script from the line named (S10)
% =========================================================================
EDGE_LO    = 0.1;   % COPIED from sessions/2026-10-16_C4-SPEC/C4_SPEC.md L202 (band, lower edge, S10)
EDGE_HI    = 0.3;   % COPIED from sessions/2026-10-16_C4-SPEC/C4_SPEC.md L203 (band, upper edge, S10)
SHARE_PCT  = 90;   % COPIED from sessions/2026-10-16_C4-SPEC/DECISIONS_C4G.md L116 (C4G-D3 v: at least this percent of the draws)
COND_MAX   = 4.5e9;   % COPIED from sessions/2026-10-16_C4-SPEC/DECISIONS_C4G.md L142 (C4G-D7 ii: RS when cond(C_J) exceeds it)
LABEL_PART1 = 'the square root of the mean, over the';   % COPIED from sessions/2026-10-16_C4-SPEC/DECISIONS_C4G.md L196 (C4G-D15 label, first part)
LABEL_PART2 = '50 combiner draws, of the per-draw compressed-domain CRB of range';   % COPIED from sessions/2026-10-16_C4-SPEC/DECISIONS_C4G.md L197 (C4G-D15 label, second part)
LABEL_RANGE = [LABEL_PART1 ' ' LABEL_PART2];                  % OPEN-1: where the label is carried
FF_SNR_DB  = 10;     % fail-fast row, C4G-D8 (iv): SNR 10 dB
FF_NRF     = 8;      % fail-fast row, C4G-D8 (iv): N_RF 8
REG_TXT    = {'below the lower edge', 'between the edges', 'above the upper edge'};
READ_TXT   = {'R1', 'R2', 'R3'};
SNR_CONV   = 'B: N0 = p_1 / 10^(SNR_dB/10)';                  % C4G-D6

% =========================================================================
%  I2; S3: configuration, the driver's lines (C4G-D1, C4G-D6)
% =========================================================================
P = setup_production_P_v4('snr', 'full');

P.d     = 2;
P.d_max = 2;

snr_vec  = -5 : 2.5 : 25;        % 13 pts  (unchanged)
nrf_vec  = [8, 16, 32, 64];       % added N_RF=64 full-array reference
n_snr    = numel(snr_vec);
n_nrf    = numel(nrf_vec);

theta_true = [30; 50] * pi / 180;
r_true     = [3.0; 7.0];
p_true     = [1.0; 1.0];

d      = numel(theta_true);
n_par  = 3*d + 1;
n_e    = 2*d;
if d ~= 2 || P.d ~= d
    error('wb_crb_d2_e3:config', 'the two-path configuration is required (d = 2).');
end
n_pts  = n_snr * n_nrf;
n_draw = P.N_seed;
if n_draw < 1 || n_draw ~= round(n_draw)
    error('wb_crb_d2_e3:config', 'P.N_seed must be a positive integer.');
end

loc_print_P(P);
fprintf('E3 P.N_seed = %.17g, P.rng_seed_W = %.17g (C4G-D3 ii, iii)\n', P.N_seed, P.rng_seed_W);
fprintf('E3 grid: %d SNR values x %d N_RF values = %d points; noise reference %s (C4G-D6)\n', ...
    n_snr, n_nrf, n_pts, SNR_CONV);
fprintf('E3 scene: theta_true [rad] = [%.17g %.17g], r_true [m] = [%.17g %.17g], p_true = [%.17g %.17g]\n', ...
    theta_true(1), theta_true(2), r_true(1), r_true(2), p_true(1), p_true(2));
if n_draw ~= 50
    fprintf('E3 NOTE: the label of C4G-D15 names 50 draws; this run has P.N_seed = %d (OPEN-9)\n', n_draw);
end

% =========================================================================
%  I9; S4; C4A-G1 item 3: reproduction gate E3-G3 (C4G-D9), before any
%  per-draw call; FAIL halts
% =========================================================================
[A_anch, anch_err] = loc_read_anchors(anchors_file, anchors_header);      % OPEN-12
if isempty(anch_err) && is_full && size(A_anch, 1) ~= n_snr * n_nrf
    anch_err = sprintf('the file has %d data rows, the gate needs %d', size(A_anch, 1), n_snr * n_nrf);
end
if ~isempty(anch_err)
    fprintf('E3-G3 reproduction [domain: %s]: cells compared 0, cells that differ 0: FAIL (%s)\n', ...
        stage, anch_err);
    fprintf('E3 HALT: E3-G3 FAIL. Run no further stage.\n');
    error('wb_crb_d2_e3:G3', 'E3-G3 FAIL: %s', anch_err);
end

if is_full
    idx_nrf = 1 : n_nrf;
    idx_snr = 1 : n_snr;
    g3_domain = sprintf('all %d rows of the anchors file', n_pts);
else
    idx_nrf = find(nrf_vec == FF_NRF);
    idx_snr = find(abs(snr_vec - FF_SNR_DB) < 1e-9);
    g3_domain = 'one row, SNR 10 dB and N_RF 8 only (fail-fast stage)';
end
if isempty(idx_nrf) || isempty(idx_snr)
    error('wb_crb_d2_e3:config', 'the fail-fast row is not on the grid.');
end

g3_cols  = {'crb_r1_m', 'crb_r2_m', 'crb_theta1_deg', 'crb_theta2_deg'};   % file columns 3 to 6
g3_cells = 0;
g3_diff  = 0;
for i_nrf = idx_nrf                                                       % as driver L89 (N_RF outer)
    P_local      = P;
    P_local.N_RF = nrf_vec(i_nrf);
    for i_snr = idx_snr                                                   % as driver L95 (SNR inner)
        snr_dB = snr_vec(i_snr);
        N0     = p_true(1) / (10^(snr_dB/10));
        fprintf('E3-G3 row: N_RF = %d, SNR = %.17g dB, N0 = %.17g\n', P_local.N_RF, snr_dB, N0);

        [crb_r, crb_theta, info] = wb_crb_compressed( ...
            theta_true, r_true, p_true, N0, P_local);

        comp_vals = [crb_r(1), crb_r(2), crb_theta(1), crb_theta(2)];
        i_row = find(abs(A_anch(:, 1) - snr_dB) < 1e-9 & A_anch(:, 2) == P_local.N_RF);
        for k = 1:4
            g3_cells = g3_cells + 1;
            if numel(i_row) ~= 1                                          % OPEN-15: no unique matching row counts as a differing cell
                g3_diff = g3_diff + 1;
                fprintf('E3-G3 DIFF: row SNR %.17g dB, N_RF %d, column %s: the file holds %d matching rows\n', ...
                    snr_dB, P_local.N_RF, g3_cols{k}, numel(i_row));
            else
                t_file = sprintf('%.6g', A_anch(i_row, 2 + k));
                t_comp = sprintf('%.6g', comp_vals(k));
                if ~strcmp(t_file, t_comp)
                    g3_diff = g3_diff + 1;
                    fprintf('E3-G3 DIFF: row SNR %.17g dB, N_RF %d, column %s: file %s, computed %s\n', ...
                        snr_dB, P_local.N_RF, g3_cols{k}, t_file, t_comp);
                end
            end
        end
    end
end
if g3_diff == 0, g3_verdict = 'PASS'; else, g3_verdict = 'FAIL'; end
% OPEN-3: the exact text of every verdict line (E3-G3, E3-G1, E3-1, E3-G2)
fprintf('E3-G3 reproduction [domain: %s; default call without opts; columns crb_r1_m crb_r2_m crb_theta1_deg crb_theta2_deg compared as text at %%.6g]: cells compared %d, cells that differ %d: %s\n', ...
    g3_domain, g3_cells, g3_diff, g3_verdict);
if g3_diff ~= 0
    fprintf('E3 HALT: E3-G3 FAIL. Run no further stage.\n');
    error('wb_crb_d2_e3:G3', 'E3-G3 FAIL: %d cell(s) differ.', g3_diff);
end

if ~is_full
    % OPEN-6: the fail-fast stage runs E3-G3 at one row and nothing else (C4G-D8 iv)
    fprintf('E3 FAILFAST COMPLETE: E3-G3 PASS at the one row; no file written; Step E3-1 and the bounds were not run.\n');
    return
end

% =========================================================================
%  I3 to I8; S5 to S9: per point and per draw (C4G-D13)
% =========================================================================
kinds  = {'call_error', 'J_nonfinite', 'J_asymmetric', 'J_diag_nonpositive', ...
          'Jeff_diag_nonpositive', 'return_nonfinite'};                   % E3-G1 kinds, OPEN-8
ju_idx = find(triu(true(n_par)));                                         % OPEN-16: column-major order of the upper triangle (F3)
[ju_r, ju_c] = ind2sub([n_par n_par], ju_idx);
n_ju   = numel(ju_idx);

pt_nrf = zeros(n_pts, 1);
pt_snr = zeros(n_pts, 1);
pt_N0  = zeros(n_pts, 1);
st_seed   = nan(n_pts, n_draw);
st_J      = nan(n_pts, n_draw, n_ju);
st_cj_mall = nan(n_pts, n_draw);
st_cj_moff = nan(n_pts, n_draw);
st_cj_eig  = nan(n_pts, n_draw, n_par);
st_cj_cond = nan(n_pts, n_draw);
st_cj_ent  = nan(n_pts, n_draw, d);
st_cj_cross = nan(n_pts, n_draw, 4);
st_ce_mall = nan(n_pts, n_draw);
st_ce_moff = nan(n_pts, n_draw);
st_ce_eig  = nan(n_pts, n_draw, n_e);
st_ce_cond = nan(n_pts, n_draw);
st_ce_ent  = nan(n_pts, n_draw, d);
st_ce_cross = nan(n_pts, n_draw, 4);
st_Vr      = nan(n_pts, n_draw, d);
st_Vth     = nan(n_pts, n_draw, d);
st_cond_eq = nan(n_pts, n_draw);

n_checked = 0;
ip = 0;
for i_nrf = 1 : n_nrf                                                     % as driver L89 (N_RF outer)
    P_local      = P;
    P_local.N_RF = nrf_vec(i_nrf);
    for i_snr = 1 : n_snr                                                 % as driver L95 (SNR inner)
        snr_dB = snr_vec(i_snr);
        N0     = p_true(1) / (10^(snr_dB/10));
        ip = ip + 1;
        pt_nrf(ip) = P_local.N_RF;
        pt_snr(ip) = snr_dB;
        pt_N0(ip)  = N0;
        fprintf('E3-1 point %d of %d: N_RF = %d, SNR = %.17g dB, N0 = %.17g, draws = %d\n', ...
            ip, n_pts, P_local.N_RF, snr_dB, N0, n_draw);

        for i_draw = 1 : n_draw
            % S5: one call per draw, N_seed = 1, seed = base + draw - 1 (C4G-D13)
            seed_i = P_local.rng_seed_W + i_draw - 1;
            P_draw = P_local;
            P_draw.N_seed     = 1;
            P_draw.rng_seed_W = seed_i;
            opts = struct();
            opts.return_J   = true;
            opts.return_raw = true;
            n_checked = n_checked + 1;
            kind = '';
            msg  = '';
            info = [];
            try
                [~, ~, info] = wb_crb_compressed(theta_true, r_true, p_true, N0, P_draw, opts);
            catch ME
                kind = 'call_error';
                msg  = ME.message;
            end
            if isempty(kind)
                J = info.J_WB;
                kind = loc_check_J(J, n_par);                              % S6
            end
            if isempty(kind)
                R = loc_corr_returns(J, d);                                % S7
                if ~R.ok_e
                    kind = 'Jeff_diag_nonpositive';
                end
            end
            if isempty(kind)
                Bd = loc_bounds(P_local, theta_true, r_true, info);        % S9
                c_eq = info.raw.cond_eq_all(1);
                if ~loc_all_finite(R, Bd, c_eq)
                    kind = 'return_nonfinite';
                end
            end
            if ~isempty(kind)
                loc_g1_fail(kinds, kind, n_checked, ip, i_draw, seed_i, msg);
            end

            st_seed(ip, i_draw) = seed_i;
            st_J(ip, i_draw, :) = J(ju_idx).';
            st_cj_mall(ip, i_draw)   = R.CJ.maxabs_all;
            st_cj_moff(ip, i_draw)   = R.CJ.maxabs_off;
            st_cj_eig(ip, i_draw, :) = reshape(R.CJ.eig, 1, 1, []);
            st_cj_cond(ip, i_draw)   = R.CJ.cond;
            st_cj_ent(ip, i_draw, :) = reshape(R.CJ.ent, 1, 1, []);
            st_cj_cross(ip, i_draw, :) = reshape(R.CJ.cross, 1, 1, []);
            st_ce_mall(ip, i_draw)   = R.CJe.maxabs_all;
            st_ce_moff(ip, i_draw)   = R.CJe.maxabs_off;
            st_ce_eig(ip, i_draw, :) = reshape(R.CJe.eig, 1, 1, []);
            st_ce_cond(ip, i_draw)   = R.CJe.cond;
            st_ce_ent(ip, i_draw, :) = reshape(R.CJe.ent, 1, 1, []);
            st_ce_cross(ip, i_draw, :) = reshape(R.CJe.cross, 1, 1, []);
            st_Vr(ip, i_draw, :)  = reshape(Bd.Vr, 1, 1, []);
            st_Vth(ip, i_draw, :) = reshape(Bd.Vth, 1, 1, []);
            st_cond_eq(ip, i_draw) = c_eq;
        end
    end
end
fprintf('E3-G1 validity [domain: every draw at every point, %d points x %d draws]: draws checked %d, draws invalid 0 (by kind: %s): PASS\n', ...
    n_pts, n_draw, n_checked, loc_kind_zero(kinds));

% =========================================================================
%  I7; S10: Step E3-1 readings per point and overall (C4G-D2, C4G-D3 iv-v,
%  C4G-D4, C4G-D7, C4G-D17, C4G-D18)
% =========================================================================
q_mat  = {'CJ7full_FIM', 'CJ7full_FIM', 'CJe4_nuisance_eliminated', 'CJe4_nuisance_eliminated'};
q_ent  = {'omega1_kappa1', 'omega2_kappa2', 'omega1_kappa1', 'omega2_kappa2'};
med_v  = nan(n_pts, 4);
med_rg = zeros(n_pts, 4);
n_in_v = zeros(n_pts, 4);
sh_ok  = false(n_pts, 4);
maxc_cj = nan(n_pts, 1);
maxc_ce = nan(n_pts, 1);
rs_pt   = false(n_pts, 1);
reading = cell(n_pts, 1);
n_rs = 0; n_r1 = 0; n_r2 = 0; n_r3 = 0; n_ru = 0;
for ip = 1 : n_pts
    vals = zeros(4, n_draw);
    vals(1, :) = abs(st_cj_ent(ip, :, 1));
    vals(2, :) = abs(st_cj_ent(ip, :, 2));
    vals(3, :) = abs(st_ce_ent(ip, :, 1));
    vals(4, :) = abs(st_ce_ent(ip, :, 2));
    for q = 1:4
        med_v(ip, q)  = median(vals(q, :));                                % C4G-D3 iv: median; OPEN-13
        med_rg(ip, q) = loc_region(med_v(ip, q), EDGE_LO, EDGE_HI);        % C4G-D17
        n_in = 0;
        for i_draw = 1:n_draw
            if loc_region(vals(q, i_draw), EDGE_LO, EDGE_HI) == med_rg(ip, q)
                n_in = n_in + 1;
            end
        end
        n_in_v(ip, q) = n_in;
        sh_ok(ip, q)  = (n_in * 100 >= SHARE_PCT * n_draw);                % C4G-D3 v
    end
    maxc_cj(ip) = max(st_cj_cond(ip, :));
    maxc_ce(ip) = max(st_ce_cond(ip, :));
    rs_pt(ip)   = loc_over(maxc_cj(ip), COND_MAX) || loc_over(maxc_ce(ip), COND_MAX);   % C4G-D7 i-iii
    if rs_pt(ip)                                                           % OPEN-10: RS read first and alone
        reading{ip} = 'RS';
        n_rs = n_rs + 1;
    elseif all(med_rg(ip, :) == med_rg(ip, 1)) && all(sh_ok(ip, :))
        reading{ip} = READ_TXT{med_rg(ip, 1)};
        if med_rg(ip, 1) == 1, n_r1 = n_r1 + 1; elseif med_rg(ip, 1) == 2, n_r2 = n_r2 + 1; else, n_r3 = n_r3 + 1; end
    else
        reading{ip} = 'RU';
        n_ru = n_ru + 1;
    end

    for q = 1:4
        fprintf('E3-1 point %d [N_RF = %d, SNR = %.17g dB; %d draws] matrix %s, entry (%s), summary median over draws of |entry| = %s, region of the median: %s, draws in that region %d of %d (share rule %s)\n', ...
            ip, pt_nrf(ip), pt_snr(ip), n_draw, q_mat{q}, q_ent{q}, sprintf('%.6g', med_v(ip, q)), ...
            REG_TXT{med_rg(ip, q)}, n_in_v(ip, q), n_draw, loc_yesno(sh_ok(ip, q), 'met', 'not met'));
    end
    fprintf('E3-1 point %d stop-rule quantity cond(C_J) [MATLAB cond; max over draws]: matrix CJ7full_FIM %s, matrix CJe4_nuisance_eliminated %s, threshold %s: %s\n', ...
        ip, sprintf('%.6g', maxc_cj(ip)), sprintf('%.6g', maxc_ce(ip)), sprintf('%.6g', COND_MAX), ...
        loc_yesno(rs_pt(ip), 'NEAR-SINGULAR', 'not near-singular'));
    fprintf('E3-1 point %d reading: %s\n', ip, reading{ip});
end
fprintf('E3-1 verdict overall [domain: %d points, each point read on its own (C4G-D2), %d draws per point]: points near-singular (RS) %d, R1 %d, R2 %d, R3 %d, RU %d\n', ...
    n_pts, n_draw, n_rs, n_r1, n_r2, n_r3, n_ru);
if n_rs > 0
    fprintf('E3-1 STOP (RS): C_J is near-singular at %d point(s). E3 stops and returns to the author; no region is read; the run continues and writes everything (C4G-D7 iv).\n', n_rs);
end

% =========================================================================
%  I9; S11: E3-G2 inversion reliability; reports, never halts (C4A-D2)
% =========================================================================
maxc_eq = nan(n_pts, 1);
g2_over = false(n_pts, 1);
for ip = 1 : n_pts
    maxc_eq(ip) = max(st_cond_eq(ip, :));
    g2_over(ip) = loc_over(maxc_eq(ip), COND_MAX);
end
n_g2_fail = sum(g2_over);                                                  % OPEN-14: PASS when no point is over the threshold
if n_g2_fail == 0, g2_verdict = 'PASS'; else, g2_verdict = 'FAIL'; end
fprintf('E3-G2 inversion reliability [domain: per point, all draws; quantity: cond of the equilibrated matrix Jn of loc_inv_equilibrated, helper MATLAB cond, from info.raw.cond_eq_all]: max over points of the max over draws %s, threshold %s, points over the threshold %d of %d: %s (reports, never halts, C4A-D2)\n', ...
    sprintf('%.6g', max(maxc_eq)), sprintf('%.6g', COND_MAX), n_g2_fail, n_pts, g2_verdict);

% =========================================================================
%  I9; S12: aggregation over the draws, variance domain (C4G-D15)
% =========================================================================
agg_r  = nan(n_pts, d);
agg_th = nan(n_pts, d);
for ip = 1 : n_pts
    for ell = 1:d
        mv = mean(st_Vr(ip, :, ell));
        if mv >= 0, agg_r(ip, ell) = sqrt(mv); end                         % OPEN-8
        mt = mean(st_Vth(ip, :, ell));
        if mt >= 0, agg_th(ip, ell) = sqrt(mt); end
    end
end
fprintf('E3 range aggregate label (C4G-D15): %s\n', LABEL_RANGE);

% =========================================================================
%  I10; S13: files F1, F2, F3 (C4G-D12, C4G-D19)
% =========================================================================
f1_name = fullfile(results_dir, ['wb_crb_d2_e3_perdraw_'  stamp '.csv']);
f2_name = fullfile(results_dir, ['wb_crb_d2_e3_perpoint_' stamp '.csv']);
f3_name = fullfile(results_dir, ['wb_crb_d2_e3_fimdump_'  stamp '.csv']);

% ---- F1: one row per point and draw (SPEC e, columns of F1) ----
% OPEN-4: how C4G and C4A narrow the columns of SPEC e (no pseudoinverse
% columns, both matrices, F3 written, angle recorded, names carry matrix and helper)
rows1 = {};
hdr1  = {};
r1 = 0;
for ip = 1 : n_pts
    for i_draw = 1 : n_draw
        nm = {}; vv = {};
        [nm, vv] = loc_put(nm, vv, 'point_id', ip);
        [nm, vv] = loc_put(nm, vv, 'd', d);
        [nm, vv] = loc_put(nm, vv, 'M', P.M);
        [nm, vv] = loc_put(nm, vv, 'N_RF', pt_nrf(ip));
        [nm, vv] = loc_put(nm, vv, 'K_s', P.K_s);
        [nm, vv] = loc_put(nm, vv, 'N_snapshots', P.N);
        [nm, vv] = loc_put(nm, vv, 'bandwidth_Hz', P.B);
        [nm, vv] = loc_put(nm, vv, 'snr_dB', pt_snr(ip));
        [nm, vv] = loc_put(nm, vv, 'snr_convention', SNR_CONV);
        [nm, vv] = loc_put(nm, vv, 'N0_linear', pt_N0(ip));
        [nm, vv] = loc_put(nm, vv, 'p1_linear', p_true(1));
        [nm, vv] = loc_put(nm, vv, 'p2_linear', p_true(2));
        [nm, vv] = loc_put(nm, vv, 'theta1_deg', theta_true(1) * 180 / pi);
        [nm, vv] = loc_put(nm, vv, 'r1_m', r_true(1));
        [nm, vv] = loc_put(nm, vv, 'theta2_deg', theta_true(2) * 180 / pi);
        [nm, vv] = loc_put(nm, vv, 'r2_m', r_true(2));
        [nm, vv] = loc_put(nm, vv, 'draw_index', i_draw);
        [nm, vv] = loc_put(nm, vv, 'seed', st_seed(ip, i_draw));
        [nm, vv] = loc_put(nm, vv, 'validity_flag_E3G1_1_valid', 1);
        % matrix CJ7full_FIM (S7 a)
        mt = 'CJ7full_FIM';
        [nm, vv] = loc_put(nm, vv, [mt '.maxabs_CJ_minus_I_all_entries.matlab_max_abs.dimless'], st_cj_mall(ip, i_draw));
        [nm, vv] = loc_put(nm, vv, [mt '.maxabs_CJ_minus_I_offdiagonal_entries.matlab_max_abs.dimless'], st_cj_moff(ip, i_draw));
        for k = 1:n_par
            [nm, vv] = loc_put(nm, vv, sprintf('%s.eig%02d_ascending.matlab_eig_sorted.dimless', mt, k), st_cj_eig(ip, i_draw, k));
        end
        [nm, vv] = loc_put(nm, vv, [mt '.cond.matlab_cond.dimless'], st_cj_cond(ip, i_draw));
        [nm, vv] = loc_put(nm, vv, [mt '.stoprule_cond_over_threshold_flag.text_at_6g.flag'], loc_over(st_cj_cond(ip, i_draw), COND_MAX));
        for l = 1:d
            [nm, vv] = loc_put(nm, vv, sprintf('%s.entry_omega%d_kappa%d_signed.code_sign_convention.dimless', mt, l, l), st_cj_ent(ip, i_draw, l));
            [nm, vv] = loc_put(nm, vv, sprintf('%s.entry_omega%d_kappa%d_abs.abs_of_signed.dimless', mt, l, l), abs(st_cj_ent(ip, i_draw, l)));
        end
        cr_nm = {'omega1_omega2', 'omega1_kappa2', 'omega2_kappa1', 'kappa1_kappa2'};
        for k = 1:4
            [nm, vv] = loc_put(nm, vv, sprintf('%s.cross_path_entry_%s_signed_recorded_only.code_sign_convention.dimless', mt, cr_nm{k}), st_cj_cross(ip, i_draw, k));
        end
        % matrix CJe4_nuisance_eliminated (S7 b)
        mt = 'CJe4_nuisance_eliminated';
        [nm, vv] = loc_put(nm, vv, [mt '.maxabs_CJ_minus_I_all_entries.matlab_max_abs.dimless'], st_ce_mall(ip, i_draw));
        [nm, vv] = loc_put(nm, vv, [mt '.maxabs_CJ_minus_I_offdiagonal_entries.matlab_max_abs.dimless'], st_ce_moff(ip, i_draw));
        for k = 1:n_e
            [nm, vv] = loc_put(nm, vv, sprintf('%s.eig%02d_ascending.matlab_eig_sorted.dimless', mt, k), st_ce_eig(ip, i_draw, k));
        end
        [nm, vv] = loc_put(nm, vv, [mt '.cond.matlab_cond.dimless'], st_ce_cond(ip, i_draw));
        [nm, vv] = loc_put(nm, vv, [mt '.stoprule_cond_over_threshold_flag.text_at_6g.flag'], loc_over(st_ce_cond(ip, i_draw), COND_MAX));
        for l = 1:d
            [nm, vv] = loc_put(nm, vv, sprintf('%s.entry_omega%d_kappa%d_signed.code_sign_convention.dimless', mt, l, l), st_ce_ent(ip, i_draw, l));
            [nm, vv] = loc_put(nm, vv, sprintf('%s.entry_omega%d_kappa%d_abs.abs_of_signed.dimless', mt, l, l), abs(st_ce_ent(ip, i_draw, l)));
        end
        for k = 1:4
            [nm, vv] = loc_put(nm, vv, sprintf('%s.cross_path_entry_%s_signed_recorded_only.code_sign_convention.dimless', mt, cr_nm{k}), st_ce_cross(ip, i_draw, k));
        end
        % bounds (S9) and the quantity of E3-G2
        for l = 1:d
            [nm, vv] = loc_put(nm, vv, sprintf('range_variance_bound_path%d.full_Jacobian_equilibrated_exact_inverse_info_raw_eq.m2', l), st_Vr(ip, i_draw, l));
        end
        for l = 1:d
            [nm, vv] = loc_put(nm, vv, sprintf('angle_variance_bound_path%d.equilibrated_exact_inverse_info_raw_eq_recorded_not_reported.deg2', l), st_Vth(ip, i_draw, l));
        end
        [nm, vv] = loc_put(nm, vv, 'cond_Jn_equilibrated.loc_inv_equilibrated_matlab_cond_info_raw_cond_eq_all.dimless', st_cond_eq(ip, i_draw));
        r1 = r1 + 1;
        if r1 == 1
            hdr1 = nm;
            rows1 = cell(n_pts * n_draw, numel(vv));
        elseif numel(nm) ~= numel(hdr1)
            error('wb_crb_d2_e3:csv', 'F1 row %d has %d fields, the header has %d.', r1, numel(nm), numel(hdr1));
        end
        rows1(r1, :) = vv;
    end
end

% ---- F2: one row per point ----
rows2 = {};
hdr2  = {};
for ip = 1 : n_pts
    nm = {}; vv = {};
    [nm, vv] = loc_put(nm, vv, 'point_id', ip);
    [nm, vv] = loc_put(nm, vv, 'd', d);
    [nm, vv] = loc_put(nm, vv, 'M', P.M);
    [nm, vv] = loc_put(nm, vv, 'N_RF', pt_nrf(ip));
    [nm, vv] = loc_put(nm, vv, 'K_s', P.K_s);
    [nm, vv] = loc_put(nm, vv, 'N_snapshots', P.N);
    [nm, vv] = loc_put(nm, vv, 'bandwidth_Hz', P.B);
    [nm, vv] = loc_put(nm, vv, 'snr_dB', pt_snr(ip));
    [nm, vv] = loc_put(nm, vv, 'snr_convention', SNR_CONV);
    [nm, vv] = loc_put(nm, vv, 'N0_linear', pt_N0(ip));
    [nm, vv] = loc_put(nm, vv, 'n_draws', n_draw);
    for q = 1:4
        bs = sprintf('%s.entry_%s_abs', q_mat{q}, q_ent{q});
        [nm, vv] = loc_put(nm, vv, [bs '.median_over_draws.dimless'], med_v(ip, q));
        [nm, vv] = loc_put(nm, vv, [bs '.median_text_at_6g.text'], sprintf('%.6g', med_v(ip, q)));
        [nm, vv] = loc_put(nm, vv, [bs '.region_of_median_1below_2between_3above.code'], med_rg(ip, q));
        [nm, vv] = loc_put(nm, vv, [bs '.draws_in_region_of_median.count'], n_in_v(ip, q));
        [nm, vv] = loc_put(nm, vv, [bs '.share_rule_met_C4G_D3_v.flag'], sh_ok(ip, q));
    end
    [nm, vv] = loc_put(nm, vv, 'CJ7full_FIM.cond.max_over_draws.matlab_cond.dimless', maxc_cj(ip));
    [nm, vv] = loc_put(nm, vv, 'CJe4_nuisance_eliminated.cond.max_over_draws.matlab_cond.dimless', maxc_ce(ip));
    [nm, vv] = loc_put(nm, vv, 'stoprule_threshold.C4G_D7_ii.dimless', COND_MAX);
    [nm, vv] = loc_put(nm, vv, 'stoprule_RS_point_flag.any_draw_either_matrix.flag', rs_pt(ip));
    [nm, vv] = loc_put(nm, vv, 'reading_R1_R2_R3_RU_RS.C4_SPEC_b5.text', reading{ip});
    [nm, vv] = loc_put(nm, vv, 'cond_Jn_equilibrated.max_over_draws.loc_inv_equilibrated_matlab_cond.dimless', maxc_eq(ip));
    [nm, vv] = loc_put(nm, vv, 'E3G2_over_threshold_flag.C4G_D7_ii.flag', g2_over(ip));
    for l = 1:d
        [nm, vv] = loc_put(nm, vv, sprintf('range_aggregate_path%d.sqrt_of_mean_over_draws_of_per_draw_variance_bound_full_Jacobian_equilibrated.m', l), agg_r(ip, l));
    end
    for l = 1:d
        [nm, vv] = loc_put(nm, vv, sprintf('angle_aggregate_path%d.sqrt_of_mean_over_draws_of_per_draw_variance_bound_recorded_not_reported.deg', l), agg_th(ip, l));
    end
    [nm, vv] = loc_put(nm, vv, 'range_aggregate_label.C4G_D15.text', LABEL_RANGE);
    if ip == 1
        hdr2 = nm;
        rows2 = cell(n_pts, numel(vv));
    elseif numel(nm) ~= numel(hdr2)
        error('wb_crb_d2_e3:csv', 'F2 row %d has %d fields, the header has %d.', ip, numel(nm), numel(hdr2));
    end
    rows2(ip, :) = vv;
end

loc_write_csv(f1_name, hdr1, rows1);
loc_write_csv(f2_name, hdr2, rows2);

% ---- F3: every entry of the upper triangle of each draw's J_WB ----
fid3 = fopen(f3_name, 'w');
if fid3 < 0
    error('wb_crb_d2_e3:file', 'cannot open %s for writing.', f3_name);
end
fprintf(fid3, 'point_id,N_RF,snr_dB,draw_index,seed,row,col,J_WB_entry.upper_triangle_column_major.raw_FIM_omega_kappa_p_N0.dimension_of_the_FIM\n');
for ip = 1 : n_pts
    for i_draw = 1 : n_draw
        for k = 1 : n_ju
            fprintf(fid3, '%d,%d,%.17g,%d,%.17g,%d,%d,%.17g\n', ip, pt_nrf(ip), pt_snr(ip), i_draw, ...
                st_seed(ip, i_draw), ju_r(k), ju_c(k), st_J(ip, i_draw, k));
        end
    end
end
fclose(fid3);

% =========================================================================
%  I10; S14: completion line
% =========================================================================
fprintf('E3 COMPLETE: stage full, %d points x %d draws; files written: %s ; %s ; %s\n', ...
    n_pts, n_draw, f1_name, f2_name, f3_name);
end


% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function loc_print_P(P)
% I2: every scalar field at %.17g, class and size of every other field
fprintf('--- E3 effective parameter struct P (every scalar field at %%.17g) ---\n');
fn = fieldnames(P);
for i = 1:numel(fn)
    v = P.(fn{i});
    if (isnumeric(v) || islogical(v)) && isscalar(v)
        fprintf('P.%s = %.17g\n', fn{i}, double(v));
    else
        fprintf('P.%s : %s, size %s\n', fn{i}, class(v), mat2str(size(v)));
    end
end
fprintf('--- end of P ---\n');
end

function [A, err] = loc_read_anchors(fname, hdr_expected)
% S4, M4: the anchors file read as text; the header must be the expected one
A = [];
err = '';
fid = fopen(fname, 'r');
if fid < 0
    err = 'cannot open the anchors file';
    return
end
cl = onCleanup(@() fclose(fid));
ln = fgetl(fid);
if ~ischar(ln) || ~strcmp(strtrim(ln), hdr_expected)
    err = 'the header line differs from the expected header';
    return
end
n_cols = numel(strsplit(hdr_expected, ','));
rows = zeros(0, n_cols);
ln = fgetl(fid);
while ischar(ln)
    ln = strtrim(ln);
    if ~isempty(ln)
        v = str2double(strsplit(ln, ','));
        if numel(v) ~= n_cols || any(~isfinite(v))
            err = 'a data row does not have the header''s number of finite fields';
            return
        end
        rows(end + 1, :) = v; %#ok<AGROW>
    end
    ln = fgetl(fid);
end
A = rows;
end

function r = loc_region(v, lo, hi)
% C4G-D17: judged on the text at %.6g; the edges belong to the middle region
v6 = str2double(sprintf('%.6g', v));
if v6 < lo
    r = 1;
elseif v6 > hi
    r = 3;
else
    r = 2;
end
end

function tf = loc_over(v, thr)
% C4G-D7 ii: judged on the text at %.6g (OPEN-7)
tf = (str2double(sprintf('%.6g', v)) > thr);
end

function s = loc_yesno(tf, a, b)
if tf, s = a; else, s = b; end
end

function kind = loc_check_J(J, n_par)
% S6: J_WB finite, equal to its transpose, positive diagonal
kind = '';
if ~isnumeric(J) || ~isequal(size(J), [n_par n_par]) || ~isreal(J) || any(~isfinite(J(:)))
    kind = 'J_nonfinite';
elseif ~isequal(J, J.')
    kind = 'J_asymmetric';
elseif ~all(diag(J) > 0)
    kind = 'J_diag_nonpositive';
end
end

function R = loc_corr_returns(J, d)
% S7, M2: the correlation forms as the GLOBECOM report function loc_CJ_report
% forms them (C4_CODE_READING.md K46, L610-614 and L635-640 of its source):
%   (a) CJ  = correlation form of the full FIM J_WB, 3d+1 by 3d+1;
%   (b) CJe = correlation form of the block of (omega, kappa) left after
%       (p_1..p_d, N0) are eliminated by a Schur complement, 2d by 2d.
% Parameter order as wb_crb_compressed.m L19: omega_1..omega_d, kappa_1..kappa_d,
% p_1..p_d, N0.
n  = size(J, 1);
J  = (J + J.') / 2;
dg = sqrt(diag(J));
CJ = J ./ (dg * dg.');
CJ = (CJ + CJ.') / 2;
R.CJ = loc_matrix_returns(CJ, d);

xi = 1:(2*d);
nu = (2*d+1):n;
Jeff = J(xi, xi) - J(xi, nu) * (J(nu, nu) \ J(nu, xi));
Jeff = (Jeff + Jeff.') / 2;
de = sqrt(diag(Jeff));
R.ok_e = isreal(de) && all(isfinite(de)) && all(de > 0);
if R.ok_e
    CJe = Jeff ./ (de * de.');
    R.CJe = loc_matrix_returns(CJe, d);
else
    R.CJe = loc_matrix_returns(nan(2*d), d);
end
end

function S = loc_matrix_returns(C, d)
% S7: D-08's three returns and the entries, for one correlation-form matrix
n = size(C, 1);
S.maxabs_all = NaN;
S.maxabs_off = NaN;
S.eig        = nan(n, 1);
S.cond       = NaN;
S.ent        = nan(1, d);
S.cross      = nan(1, 4);
if ~all(isfinite(C(:)))
    return
end
S.maxabs_all = max(abs(C - eye(n)), [], 'all');       % as K46 (all entries); OPEN-2
E = abs(C - eye(n));
E(1:n+1:end) = 0;
S.maxabs_off = max(E, [], 'all');                     % off-diagonal entries only; OPEN-2
ev = eig(C);
if ~isreal(ev)
    return
end
S.eig  = sort(ev, 'ascend');
S.cond = cond(C);                                     % MATLAB cond (2-norm), named in the output
for l = 1:d
    S.ent(l) = C(l, d + l);                           % (omega_l, kappa_l), code sign convention
end
S.cross = [C(1, 2), C(1, d + 2), C(2, d + 1), C(d + 1, d + 2)];   % cross-path entries, recorded only
end

function B = loc_bounds(P, theta_true, r_true, info)
% S9, M3: variance bounds from info.raw.eq (diagonally equilibrated exact
% inverse of the draw's J_WB) by the mapping of BoundSpec Sec. 3 (L313-318):
%   V_r,FJ  = g' C_(omega,kappa) g, g = [-2 r cos(theta)/(c1 sin^2 theta); -r/kappa]   [m^2]
%   V_theta = max(C_oo, 0) / (c1 sin theta)^2 * (180/pi)^2                            [deg^2]
%   c1 = 2 pi d_ant / lambda_c, c0 = pi d_ant^2 / lambda_c, kappa = c0 sin^2(theta) / r
d  = numel(theta_true);
c1 = 2 * pi * P.d_ant / P.lambda_c;
c0 = pi * P.d_ant^2 / P.lambda_c;
B.Vr  = nan(1, d);
B.Vth = nan(1, d);
for ell = 1:d
    th  = theta_true(ell);
    r   = r_true(ell);
    kap = c0 * sin(th)^2 / r;
    g   = [-2 * r * cos(th) / (c1 * sin(th)^2); -r / kap];
    Coo = info.raw.eq.C_oo(ell, 1);
    Ckk = info.raw.eq.C_kk(ell, 1);
    Cok = info.raw.eq.C_ok(ell, 1);
    B.Vr(ell)  = g.' * [Coo Cok; Cok Ckk] * g;
    B.Vth(ell) = max(Coo, 0) / (c1 * sin(th))^2 * (180 / pi)^2;
end
end

function tf = loc_all_finite(R, B, c_eq)
% S6: every return of S7 and S9 is finite
v = [R.CJ.maxabs_all, R.CJ.maxabs_off, R.CJ.eig(:).', R.CJ.cond, R.CJ.ent, R.CJ.cross, ...
     R.CJe.maxabs_all, R.CJe.maxabs_off, R.CJe.eig(:).', R.CJe.cond, R.CJe.ent, R.CJe.cross, ...
     B.Vr, B.Vth, c_eq];
tf = isreal(v) && all(isfinite(v));
end

function s = loc_kind_zero(kinds)
parts = cell(1, numel(kinds));
for i = 1:numel(kinds)
    parts{i} = sprintf('%s 0', kinds{i});
end
s = strjoin(parts, ', ');
end

function loc_g1_fail(kinds, kind, n_checked, ip, i_draw, seed_i, msg)
% E3-G1: any invalid draw is FAIL and halts (C4G-D8 ii)
parts = cell(1, numel(kinds));
for i = 1:numel(kinds)
    parts{i} = sprintf('%s %d', kinds{i}, double(strcmp(kinds{i}, kind)));
end
fprintf('E3-G1 validity [domain: every draw at every point]: draws checked %d, draws invalid 1 (by kind: %s): FAIL\n', ...
    n_checked, strjoin(parts, ', '));
fprintf('E3-G1 first invalid draw: point %d, draw %d, seed %.17g, kind %s, message: %s\n', ip, i_draw, seed_i, kind, msg);
% OPEN-5: when E3-G1 halts, no file is written (the files are written at the end)
fprintf('E3 HALT: E3-G1 FAIL. Run no further stage; no file written.\n');
error('wb_crb_d2_e3:G1', 'E3-G1 FAIL at point %d, draw %d (%s).', ip, i_draw, kind);
end

function [nm, vv] = loc_put(nm, vv, name, val)
nm{end + 1} = name;
vv{end + 1} = val;
end

function s = loc_csvfield(t)
if any(t == ',') || any(t == '"')
    s = ['"' strrep(t, '"', '""') '"'];
else
    s = t;
end
end

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
