function t53_randw_residual(stage)
%T53_RANDW_RESIDUAL  T-53 (Lane A item A5): random-W full-Jacobian residual,
%  with its distribution, at Paper C's own configuration.
%
%  Spec section  : PaperC_T53_T54_Specs.md Sec. 3 (T-53), AS AMENDED BY
%                  PaperC_T53_T54_Specs_AddendumB.md (N-35, approved
%                  2026-09-23). Where the two disagree the amended spec
%                  wins.
%  Gate IDs       : PC-T1, PC-T2, G-R1a, G-R1b, G-R1c, G-C1, G-R5, G-CND,
%                  G-NS, G-DST (Sec. 3.7); G-SEED (Addendum B.1); R4-SQ
%                  (Addendum B.2); R5c's reproduction check (Addendum B.3).
%  Blocks         :
%    stage dispatch          -- this header switch
%    PRE-FLIGHT (R0)         -- Sec. 3.7 PC-T1, PC-T2 (line counts, config)
%    R1  REPRO               -- Sec. 3.2-3.3; gates G-R1a/b/c
%  G-R1a/b compare in absolute dB and G-R1c compares integer counts per N-38 (PaperC_A5_R1_GateFailure_Adjudication.md); original verdicts printed as G-R1a-orig, G-R1b-orig, G-R1c-orig.
%    R2  PRIMARY              -- Sec. 3.2, 3.4; the blind measurement
%    R3  GRID                 -- Sec. 3.2 (T-49 grid)
%    R4  NRF                  -- Sec. 3.2, Addendum B.2 (P2b fit set, R4-SQ)
%    R5  SELFCHK               -- Sec. 3.2, 3.7 (G-R5, G-C1); production grid
%    R5c CENTRED SELFCHK      -- Addendum B.3 item 4 (new)
%    R6  SNR                  -- Sec. 3.2; reported, never gated
%    R7  SEEDCHK              -- Addendum B.7 (new); reported, never gated
%
%  CONVENTIONS (Spec Sec. 2.1, non-negotiable):
%    - inverse: equilibrated only. info.raw.eq is the sole reported source;
%      the function's crb_r/crb_theta returns and info.raw.pinv are NEVER
%      used for a reported quantity (the pinv-minus-eq column at R2 is the
%      one exception, at the summary level, labelled).
%    - call wb_crb_compressed via opts.W_list and opts.return_raw only.
%      NO edit is made to wb_crb_compressed.m.
%    - c0 = 3e8 (not 299792458); d = 1; p_true = 1; N0 = 10^(-SNR/10).
%    - alpha_sub for every Paper C leg is the PRODUCTION grid
%      (Addendum B.3): k_idx = (-(K_s/2):(K_s/2-1)).',
%      alpha_k = 1 + k_idx*Delta_f/f_c, i.e. setup_production_P_v4('snr','full')
%      L161-163. The GLOBECOM leg (R1) uses its own native grid via the
%      VERBATIM loc_subcarrier_grid copy (Spec Sec. 3.3), not this one.
%
%  Author-side execution. MATLAB R2025b. Provenance P-RUN.
%  Encoding: 7-bit ASCII only.
% =========================================================================

if nargin < 1 || isempty(stage)
    error('t53_randw_residual: stage argument required. Usage: t53_randw_residual(stage), stage in {R0,R1,R2,R3,R4,R5,R5c,R6,R7,ALL}.');
end
valid_stages = {'R0','R1','R2','R3','R4','R5','R5c','R6','R7','ALL'};
assert(any(strcmp(stage, valid_stages)), ...
    't53_randw_residual: stage must be one of R0,R1,R2,R3,R4,R5,R5c,R6,R7,ALL.');

OUT_DIR = 'results/paperC_A5_t53';
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end
ts = datestr(now, 'yyyymmdd_HHMMSS');
diary_file = fullfile(OUT_DIR, sprintf('A5_console_%s.txt', ts));
diary(diary_file); diary on;
t_start = tic;

fprintf('=============================================================\n');
fprintf('  T-53 (A5) -- random-W full-Jacobian residual, with distribution\n');
fprintf('  stage = %s\n', stage);
fprintf('  Addendum B (N-35, approved 2026-09-23) IS BINDING\n');
fprintf('=============================================================\n\n');

% ---- persistent output accumulators (shared across stages in 'ALL') ------
DRAW_ROWS    = {};   % -> t53_draws_<ts>.csv
SUMMARY_ROWS = {};   % -> t53_summary_<ts>.csv
GATE_ROWS    = {};   % -> t53_gates_<ts>.csv
REPRO_ROWS   = {};   % -> t53_repro_<ts>.csv  (R1 only)

% ---- physical constants, locked (Spec Sec. 3.2) ---------------------------
c0_const = 3e8;                 % NOT 299792458
fc_const = 28e9;
lam_c    = c0_const / fc_const;
d_ant    = lam_c / 2;

% =========================================================================
%  R0 -- PRE-FLIGHT: PC-T1 (line counts), PC-T2 (config block dry stamp)
% =========================================================================
if any(strcmp(stage, {'R0','ALL'}))
    fprintf('### R0 -- PRE-FLIGHT (PC-T1, PC-T2) ###\n\n');

    pct1_files = struct( ...
        'name', {'wb_crb_compressed.m', 't49_f032_fulljac_invariance.m', 'stageA_s0_audit.m'}, ...
        'expected', {431, NaN, 649});
    pc_t1_pass = true;
    for i = 1:numel(pct1_files)
        fp = which(pct1_files(i).name);
        if isempty(fp)
            fprintf('PC-T1 %-32s NOT ON PATH\n', pct1_files(i).name);
            pc_t1_pass = false;
            continue
        end
        txt = fileread(fp);
        % Line count matches wc -l convention (newline-terminated lines
        % only): a trailing unterminated line is NOT counted as an extra
        % line. This matches the convention PC-T1's targets (431, 649)
        % were set from -- do not add 1 for a missing trailing newline.
        nlines = sum(txt == newline);
        if isnan(pct1_files(i).expected)
            fprintf('PC-T1 %-32s lines=%d  (as uploaded, no fixed target)  %s\n', ...
                pct1_files(i).name, nlines, fp);
        else
            ok = (nlines == pct1_files(i).expected);
            fprintf('PC-T1 %-32s lines=%d  expected=%d  %s\n', ...
                pct1_files(i).name, nlines, pct1_files(i).expected, ...
                loc_pf(ok));
            pc_t1_pass = pc_t1_pass && ok;
        end
    end
    fprintf('PC-T1 VERDICT: %s\n\n', loc_pf(pc_t1_pass));
    GATE_ROWS(end+1,:) = {'PC-T1','R0','ALL',NaN,NaN,loc_pf(pc_t1_pass)}; %#ok<SAGROW>

    % PC-T2 dry stamp: print the production-grid config block so the
    % printed-config discipline is visible even before any leg runs.
    Pchk = setup_production_P_v4('snr','full');
    k_idx_chk = (-(Pchk.K_s/2):(Pchk.K_s/2-1)).';
    mean_delta1_chk = mean(k_idx_chk * Pchk.Delta_f / Pchk.fc);
    mean_delta2_chk = mean((k_idx_chk * Pchk.Delta_f / Pchk.fc).^2);
    fprintf('PC-T2 DRY STAMP (production grid, Addendum B.3):\n');
    fprintf('  M=%d N_RF=%d Delta_f=%.6g Hz K_s=%d fc=%.6g Hz\n', ...
        Pchk.M, Pchk.N_RF, Pchk.Delta_f, Pchk.K_s, Pchk.fc);
    fprintf('  mean_delta1 = %.9e   mean_delta2 = %.9e\n', ...
        mean_delta1_chk, mean_delta2_chk);
    fprintf('PC-T2 VERDICT: %s (full block re-printed inside every leg below)\n\n', loc_pf(true));
    GATE_ROWS(end+1,:) = {'PC-T2','R0','ALL',NaN,NaN,'PASS'}; %#ok<SAGROW>

    if ~pc_t1_pass
        fprintf('HALT: PC-T1 FAIL. Line-count drift is a finding about the\n');
        fprintf('repository, not a nuisance. Stopping before any leg runs.\n');
        loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
        diary off;
        return
    end
end

if strcmp(stage, 'R0')
    loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
    diary off;
    return
end

% =========================================================================
%  Shared setup: Paper C production P-struct and grid (Addendum B.3)
% =========================================================================
P_PC = setup_production_P_v4('snr', 'full');   % M=64, N_RF=8, K_s=16 @ B=400MHz default
% Note: P_PC.alpha_k_vec/K_s are the DEFAULT-bandwidth (400 MHz) grid;
% every leg below re-derives its own alpha_sub per bandwidth via
% loc_prodgrid_alpha, consistent with Addendum B.3.

% cell_id counter for the amended seed rule (Addendum B.1): base 5300000,
% stride 10000 per cell, n <= 9999 within a cell.
CELL_ID = 0;

% =========================================================================
%  R1 -- REPRO: reproduce stageA_A2_ensemble.csv exactly (gate, NOT BLIND)
% =========================================================================
if any(strcmp(stage, {'R1','ALL'}))
    fprintf('### R1 -- REPRO (GLOBECOM leg, gate against stageA_A2_ensemble.csv) ###\n\n');
    fprintf('R1 reuses Stage A''s seeds and is exempt from the amended seed rule (Addendum B.1).\n');

    M1 = 256; N_RF1 = 16; Delta_f1 = 120e3;
    theta1_deg = 40; r1 = 5; SNR1_dB = 10; N_snap1 = 64; p1 = 1;
    N0_1 = p1 / (10^(SNR1_dB/10));
    theta1 = theta1_deg * pi/180;

    B1_list = [100e6, 400e6, 800e6];
    N_DRAWS1 = 385;
    seeds1 = [42, 10001:(10000 + N_DRAWS1 - 1)];   % VERBATIM, Sec. 3.3 / stageA_s0_audit.m L261

    mean_delta1_r1 = nan(1,numel(B1_list));
    mean_delta2_r1 = nan(1,numel(B1_list));

    fprintf('PC-T2 CONFIG [R1]: M=%d N_RF=%d Delta_f=%.6g fc=%.6g theta=%g deg r=%g m SNR=%g dB N=%d p=%g N0=%.6g seedbase=STAGE_A draws=%d\n', ...
        M1, N_RF1, Delta_f1, fc_const, theta1_deg, r1, SNR1_dB, N_snap1, p1, N0_1, N_DRAWS1);

    R1_ALL = cell(numel(B1_list), 1);
    for ib = 1:numel(B1_list)
        Bv = B1_list(ib);
        [al_sub1, Ks1] = loc_subcarrier_grid(Bv, fc_const, Delta_f1);   % VERBATIM
        mean_delta1_r1(ib) = mean(al_sub1(:) - 1);
        mean_delta2_r1(ib) = mean((al_sub1(:) - 1).^2);
        fprintf('  B=%.0f MHz K_s=%d mean_delta1=%.9e mean_delta2=%.9e\n', ...
            Bv/1e6, Ks1, mean_delta1_r1(ib), mean_delta2_r1(ib));

        rows_b = cell(N_DRAWS1, 1);
        for n = 1:N_DRAWS1
            s = RandStream('mt19937ar', 'Seed', seeds1(n));               %#ok<RAND>
            Wn = (1/sqrt(M1)) * exp(1j * 2*pi * rand(s, M1, N_RF1));      % VERBATIM

            [Bfj, Bdiag, Btheta, raw] = loc_B_of_gW(theta1, r1, p1, N0_1, ...
                M1, N_RF1, N_snap1, al_sub1, lam_c, d_ant, Wn);

            rows_b{n} = {'GLOBECOM', ib, 'R1', M1, N_RF1, Delta_f1, Bv, Ks1, ...
                mean_delta2_r1(ib), theta1_deg, r1, SNR1_dB, n, seeds1(n), ...
                Bfj, Bdiag, Btheta, raw.V_r_nb, raw.V_r_wb, raw.V_r_nb_diag, raw.V_r_wb_diag, ...
                raw.C_oo_nb, raw.C_kk_nb, raw.C_ok_nb, raw.C_oo_wb, raw.C_kk_wb, raw.C_ok_wb, ...
                raw.cond_eq_nb, raw.cond_eq_wb, raw.nsing_nb, raw.nsing_wb};
        end
        R1_ALL{ib} = rows_b;
        for n = 1:N_DRAWS1
            DRAW_ROWS(end+1,:) = R1_ALL{ib}{n}; %#ok<SAGROW>
        end
    end

    % ---- G-R1a, G-R1b, G-R1c: reproduce stageA_A2_ensemble.csv ------------
    ens_path = loc_find_file('stageA_A2_ensemble.csv');
    summ_path = loc_find_file('stageA_A2_summary.csv');
    assert(~isempty(ens_path), 't53_randw_residual: stageA_A2_ensemble.csv not found on path.');
    assert(~isempty(summ_path), 't53_randw_residual: stageA_A2_summary.csv not found on path.');
    Tens = readtable(ens_path);
    Tsum = readtable(summ_path);

    g_r1a_ok = true(N_DRAWS1, numel(B1_list));
    g_r1b_ok = true(N_DRAWS1, numel(B1_list));
    g_r1a_orig = true(N_DRAWS1, numel(B1_list)); g_r1b_orig = g_r1a_orig;
    sgn_ok_ex = g_r1a_orig; sgn_ok_fj = g_r1a_orig;
    absd_ex = nan(N_DRAWS1, numel(B1_list)); absd_fj = absd_ex;
    for ib = 1:numel(B1_list)
        Bv_MHz = B1_list(ib)/1e6;
        mask = abs(Tens.B_MHz - Bv_MHz) < 1e-6;
        Tb = Tens(mask, :);
        [~, ord] = sort(Tb.draw);
        Tb = Tb(ord, :);
        assert(height(Tb) == N_DRAWS1, ...
            't53_randw_residual: stageA_A2_ensemble.csv does not have %d rows at B=%g MHz.', ...
            N_DRAWS1, Bv_MHz);
        for n = 1:N_DRAWS1
            Bdiag_t53 = R1_ALL{ib}{n}{16};   % B_diag_t53 column position
            Bfj_t53   = R1_ALL{ib}{n}{15};   % B_fj_t53 column position
            G_exact_stageA = Tb.G_exact_dB(n);
            G_full_stageA  = Tb.G_fullJac_dB(n);

            reldiff_exact = abs(Bdiag_t53 - G_exact_stageA) / max(abs(G_exact_stageA), eps);
            reldiff_fj    = abs(Bfj_t53   - G_full_stageA)  / max(abs(G_full_stageA),  eps);

            absdiff_exact = abs(Bdiag_t53 - G_exact_stageA);
            absdiff_fj    = abs(Bfj_t53   - G_full_stageA);
            absd_ex(n, ib) = absdiff_exact; absd_fj(n, ib) = absdiff_fj;
            g_r1a_orig(n, ib) = (reldiff_exact <= 1e-8);   % registered rule, reported
            g_r1b_orig(n, ib) = (reldiff_fj    <= 1e-8);   % registered rule, reported
            g_r1a_ok(n, ib)   = (absdiff_exact <= 1e-8);   % N-38 rule [dB], halting
            g_r1b_ok(n, ib)   = (absdiff_fj    <= 1e-8);   % N-38 rule [dB], halting
            sgn_ok_ex(n, ib)  = (sign(Bdiag_t53) == sign(G_exact_stageA));
            sgn_ok_fj(n, ib)  = (sign(Bfj_t53)   == sign(G_full_stageA));

            REPRO_ROWS(end+1,:) = {n, seeds1(n), Bv_MHz, G_exact_stageA, Bdiag_t53, reldiff_exact, ...
                G_full_stageA, Bfj_t53, reldiff_fj, absdiff_exact, absdiff_fj}; %#ok<SAGROW>
        end
    end
    n_r1a = sum(g_r1a_ok(:)); n_r1b = sum(g_r1b_ok(:)); n_tot1 = numel(g_r1a_ok);
    n_r1a_orig = sum(g_r1a_orig(:)); n_r1b_orig = sum(g_r1b_orig(:));
    max_absd_ex = max(absd_ex(:)); max_absd_fj = max(absd_fj(:));
    n_sgn_ex = sum(sgn_ok_ex(:)); n_sgn_fj = sum(sgn_ok_fj(:));
    fprintf('\nG-R1a (B_diag vs G_exact_dB, abs <= 1e-8 dB, N-38): %d/%d  %s  max abs diff = %.6e dB\n', ...
        n_r1a, n_tot1, loc_pf(n_r1a == n_tot1), max_absd_ex);
    fprintf('G-R1b (B_fj vs G_fullJac_dB, abs <= 1e-8 dB, N-38): %d/%d  %s  max abs diff = %.6e dB\n', ...
        n_r1b, n_tot1, loc_pf(n_r1b == n_tot1), max_absd_fj);
    fprintf('G-R1a-orig (registered, rel <= 1e-8, reported only): %d/%d  %s\n', ...
        n_r1a_orig, n_tot1, loc_pf(n_r1a_orig == n_tot1));
    fprintf('G-R1b-orig (registered, rel <= 1e-8, reported only): %d/%d  %s\n', ...
        n_r1b_orig, n_tot1, loc_pf(n_r1b_orig == n_tot1));
    fprintf('R1 sign agreement: diag %d/%d, fj %d/%d\n', n_sgn_ex, n_tot1, n_sgn_fj, n_tot1);
    GATE_ROWS(end+1,:) = {'G-R1a','R1','ALL_1155',max_absd_ex,1e-8,loc_pf(n_r1a==n_tot1)}; %#ok<SAGROW>
    GATE_ROWS(end+1,:) = {'G-R1b','R1','ALL_1155',max_absd_fj,1e-8,loc_pf(n_r1b==n_tot1)}; %#ok<SAGROW>
    GATE_ROWS(end+1,:) = {'G-R1a-orig','R1','ALL_1155',n_r1a_orig/n_tot1,1e-8,loc_pf(n_r1a_orig==n_tot1)}; %#ok<SAGROW>
    GATE_ROWS(end+1,:) = {'G-R1b-orig','R1','ALL_1155',n_r1b_orig/n_tot1,1e-8,loc_pf(n_r1b_orig==n_tot1)}; %#ok<SAGROW>

    r1a_pass = (n_r1a == n_tot1); r1b_pass = (n_r1b == n_tot1);

    % G-R1c: cell summary at B=400 MHz, n=385 vs stageA_A2_summary.csv
    ib400 = find(B1_list == 400e6, 1);
    Bfj_all_400   = cell2mat(cellfun(@(c) c{15}, R1_ALL{ib400}, 'UniformOutput', false));
    Bdiag_all_400 = cell2mat(cellfun(@(c) c{16}, R1_ALL{ib400}, 'UniformOutput', false));
    row400 = Tsum(Tsum.B_MHz == 400 & Tsum.n == 385, :);
    assert(height(row400) == 1, 't53_randw_residual: stageA_A2_summary.csv row (B=400,n=385) not unique.');

    mean_t53 = mean(Bdiag_all_400); median_t53 = median(Bdiag_all_400); sd_t53 = std(Bdiag_all_400);
    fracpos_t53 = mean(Bdiag_all_400 > 0);
    k_t53 = sum(Bdiag_all_400 > 0);
    mask400 = abs(Tens.B_MHz - 400) < 1e-6;
    k_ens = sum(Tens.G_exact_dB(mask400) > 0);
    k_sum = round(row400.frac_positive * row400.n);
    reldiff_mean   = abs(mean_t53 - row400.mean_dB) / abs(row400.mean_dB);
    reldiff_median = abs(median_t53 - row400.median_dB) / abs(row400.median_dB);
    reldiff_sd     = abs(sd_t53 - row400.sd_dB) / abs(row400.sd_dB);
    g_r1c_orig = (reldiff_mean <= 1e-6) && (reldiff_median <= 1e-6) && (reldiff_sd <= 1e-6) && ...
                 (fracpos_t53 == row400.frac_positive);          % registered rule, reported
    g_r1c_ok   = (reldiff_mean <= 1e-6) && (reldiff_median <= 1e-6) && (reldiff_sd <= 1e-6) && ...
                 (k_t53 == k_ens);                                % N-38 repair
    fprintf('G-R1c (cell summary at B=400,n=385 vs stageA_A2_summary.csv): %s\n', loc_pf(g_r1c_ok));
    fprintf('  mean   t53=%.9f  stageA=%.9f  reldiff=%.3e\n', mean_t53, row400.mean_dB, reldiff_mean);
    fprintf('  median t53=%.9f  stageA=%.9f  reldiff=%.3e\n', median_t53, row400.median_dB, reldiff_median);
    fprintf('  sd     t53=%.9f  stageA=%.9f  reldiff=%.3e\n', sd_t53, row400.sd_dB, reldiff_sd);
    fprintf('  frac_positive t53=%.9f  stageA=%.9f\n', fracpos_t53, row400.frac_positive);
    fprintf('  count B>0: t53=%d  stageA_ensemble=%d  stageA_summary(round f*n)=%d\n', k_t53, k_ens, k_sum);
    fprintf('G-R1c-orig (registered, == on frac_positive, reported only): %s\n', loc_pf(g_r1c_orig));
    GATE_ROWS(end+1,:) = {'G-R1c','R1','B400_n385',NaN,1e-6,loc_pf(g_r1c_ok)}; %#ok<SAGROW>
    GATE_ROWS(end+1,:) = {'G-R1c-orig','R1','B400_n385',NaN,NaN,loc_pf(g_r1c_orig)}; %#ok<SAGROW>

    % ---- per-cell summary rows for R1 (reported with distribution, F-030) --
    for ib = 1:numel(B1_list)
        Bfj_ib   = cell2mat(cellfun(@(c) c{15}, R1_ALL{ib}, 'UniformOutput', false));
        Bdiag_ib = cell2mat(cellfun(@(c) c{16}, R1_ALL{ib}, 'UniformOutput', false));
        Btheta_ib = cell2mat(cellfun(@(c) c{17}, R1_ALL{ib}, 'UniformOutput', false));
        Vr_nb_ib = cell2mat(cellfun(@(c) c{18}, R1_ALL{ib}, 'UniformOutput', false));
        Vr_wb_ib = cell2mat(cellfun(@(c) c{19}, R1_ALL{ib}, 'UniformOutput', false));
        cond_eq_nb_ib = cell2mat(cellfun(@(c) c{28}, R1_ALL{ib}, 'UniformOutput', false));
        cond_eq_wb_ib = cell2mat(cellfun(@(c) c{29}, R1_ALL{ib}, 'UniformOutput', false));
        nsing_nb_ib = cell2mat(cellfun(@(c) c{30}, R1_ALL{ib}, 'UniformOutput', false));
        nsing_wb_ib = cell2mat(cellfun(@(c) c{31}, R1_ALL{ib}, 'UniformOutput', false));
        [~, Ks_ib] = loc_subcarrier_grid(B1_list(ib), fc_const, Delta_f1);
        SUMMARY_ROWS = loc_report_cell_all(SUMMARY_ROWS, 'GLOBECOM', ib, 'R1', M1, N_RF1, B1_list(ib), ...
            Ks_ib, theta1_deg, r1, SNR1_dB, N_DRAWS1, Bfj_ib, Bdiag_ib, Btheta_ib, ...
            max(cond_eq_nb_ib(:)), max(cond_eq_wb_ib(:)), mean_delta2_r1(ib), Vr_nb_ib, Vr_wb_ib, ...
            max(max(nsing_nb_ib(:)), max(nsing_wb_ib(:))));
    end

    fprintf('\nR1 SUMMARY: G-R1a=%s G-R1b=%s G-R1c=%s | registered: G-R1a-orig=%s G-R1b-orig=%s G-R1c-orig=%s\n\n', ...
        loc_pf(r1a_pass), loc_pf(r1b_pass), loc_pf(g_r1c_ok), ...
        loc_pf(n_r1a_orig == n_tot1), loc_pf(n_r1b_orig == n_tot1), loc_pf(g_r1c_orig));

    if ~(r1a_pass && r1b_pass)
        fprintf('HALT: G-R1a or G-R1b FAIL. The harness is wrong and nothing\n');
        fprintf('downstream is readable. Do not run R2..R6.\n');
        loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
        diary off;
        return
    end
end

if strcmp(stage, 'R1')
    loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
    diary off;
    return
end

% =========================================================================
%  R2 -- PRIMARY: the blind measurement, Paper C configuration, 2213 draws
% =========================================================================
SEED_REGISTRY = [];   % running list of all seeds used from R2 onward (G-SEED)

if any(strcmp(stage, {'R2','ALL'}))
    fprintf('### R2 -- PRIMARY (Paper C configuration, blind, n=2213) ###\n\n');
    CELL_ID = CELL_ID + 1;
    theta2_deg = 40; r2 = 2.6578125; SNR2_dB = 10; N_draws2 = 2213;
    B2 = 400e6;

    [SUMMARY_ROWS, DRAW_ROWS, SEED_REGISTRY, GATE_ROWS] = loc_run_cell( ...
        'PaperC', CELL_ID, 'R2', P_PC, B2, theta2_deg, r2, SNR2_dB, N_draws2, ...
        lam_c, d_ant, fc_const, SEED_REGISTRY, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, true);
end

if strcmp(stage, 'R2')
    loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
    diary off;
    return
end

% =========================================================================
%  R3 -- GRID: geometry x bandwidth dependence, Paper C array, n=385/cell
% =========================================================================
if any(strcmp(stage, {'R3','ALL'}))
    fprintf('### R3 -- GRID (geometry x bandwidth, Paper C array, n=385/cell) ###\n\n');
    B3_list = [100 200 400 600 800] * 1e6;
    R3_list = [1.50 2.13 3.00 5.00 10.00];
    THETA3_list = [20 40 60];
    N_draws3 = 385; SNR3_dB = 10;

    for iB = 1:numel(B3_list)
        for iT = 1:numel(THETA3_list)
            for iR = 1:numel(R3_list)
                CELL_ID = CELL_ID + 1;
                [SUMMARY_ROWS, DRAW_ROWS, SEED_REGISTRY, GATE_ROWS] = loc_run_cell( ...
                    'PaperC', CELL_ID, 'R3', P_PC, B3_list(iB), THETA3_list(iT), R3_list(iR), ...
                    SNR3_dB, N_draws3, lam_c, d_ant, fc_const, SEED_REGISTRY, ...
                    DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, false);
            end
        end
    end
    fprintf('R3 complete: %d cells (5 B x 3 theta x 5 r).\n\n', numel(B3_list)*numel(THETA3_list)*numel(R3_list));
end

if strcmp(stage, 'R3')
    loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
    diary off;
    return
end

% =========================================================================
%  R4 -- NRF: N_RF ladder at fixed geometry, Addendum B.2 (P2b fit set, R4-SQ)
% =========================================================================
if any(strcmp(stage, {'R4','ALL'}))
    fprintf('### R4 -- NRF LADDER (Addendum B.2: fit set {4,8,16,32}; N_RF=64 -> R4-SQ) ###\n\n');
    NRF4_list = [4 8 16 32 64];
    theta4_deg = 40; r4 = 2.6578125; SNR4_dB = 10; N_draws4 = 385; B4 = 400e6;

    R4_SD = nan(size(NRF4_list));
    R4_NRF = NRF4_list;
    R4_BFJ = cell(size(NRF4_list));   % per-setting full draw vectors, for P2b's draw-level bootstrap
    for i = 1:numel(NRF4_list)
        CELL_ID = CELL_ID + 1;
        Pi = P_PC; Pi.N_RF = NRF4_list(i);
        is_sq = (NRF4_list(i) == Pi.M);   % N_RF == M -> harness check R4-SQ

        [SUMMARY_ROWS, DRAW_ROWS, SEED_REGISTRY, GATE_ROWS, Bfj_here, sd_here, meanB_here, maxcond_here, maxnsing_here] = ...
            loc_run_cell_nrf( ...
            'PaperC', CELL_ID, 'R4', Pi, B4, theta4_deg, r4, SNR4_dB, N_draws4, ...
            lam_c, d_ant, fc_const, SEED_REGISTRY, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS);
        R4_SD(i) = sd_here;
        R4_BFJ{i} = Bfj_here;

        if is_sq
            % R4-SQ: SD <= 1e-6 dB and mean within 1e-6 dB of B_ref (Addendum B.2, B.3)
            [al_sq, ~] = loc_prodgrid_alpha(B4, Pi.fc, Pi.K_s);
            B_ref_sq = 10*log10(mean(al_sq.^2));
            sq_sd_ok = (sd_here <= 1e-6);
            sq_mean_ok = (abs(meanB_here - B_ref_sq) <= 1e-6);
            fprintf('R4-SQ (N_RF=M=%d, harness self-check, non-halting): SD=%.6e dB (<=1e-6? %s)  mean=%.9e dB vs B_ref=%.9e dB (within 1e-6? %s)\n', ...
                NRF4_list(i), sd_here, loc_pf(sq_sd_ok), meanB_here, B_ref_sq, loc_pf(sq_mean_ok));
            fprintf('R4-SQ max cond_eq=%.6e  n_singular=%d\n', maxcond_here, maxnsing_here);
            GATE_ROWS(end+1,:) = {'R4-SQ','R4',CELL_ID,sd_here,1e-6,loc_pf(sq_sd_ok && sq_mean_ok)}; %#ok<SAGROW>
        end
        if NRF4_list(i) == 4
            fprintf('R4 N_RF=4 (parameter count 3d+1=4): max cond_eq=%.6e  n_singular=%d\n', ...
                maxcond_here, maxnsing_here);
        end
    end

    % ---- P2b: log-log fit of SD vs N_RF over {4,8,16,32} only (drop 64) ---
    fit_mask = ismember(R4_NRF, [4 8 16 32]);
    x = log(R4_NRF(fit_mask)); y = log(R4_SD(fit_mask));
    pfit = polyfit(x, y, 1);
    slope_hat = pfit(1);
    fprintf('\nP2b FIT (N_RF in {4,8,16,32}): log-log slope = %.6f\n', slope_hat);

    % Bootstrap over the four cells' DRAWS (Addendum B.2): resample each
    % cell's 385 draws with replacement, recompute that cell's SD from the
    % resample, refit the slope; repeat 2000 times. This resamples the
    % underlying data the SD is estimated from, not the four (N_RF,SD)
    % summary points themselves (resampling only 4 points would be
    % statistically degenerate -- with replacement over 4 items, repeated
    % picks collapse the x-spread the log-log fit needs).
    NBOOT = 2000; rng(20260922, 'twister');
    fit_idx = find(fit_mask);
    xi_nrf = R4_NRF(fit_idx);
    boot_slopes = nan(NBOOT,1);
    for ib = 1:NBOOT
        yi_boot = nan(size(xi_nrf));
        for k = 1:numel(fit_idx)
            draws_k = R4_BFJ{fit_idx(k)};
            nd = numel(draws_k);
            idx = randi(nd, 1, nd);
            yi_boot(k) = std(draws_k(idx));
        end
        pb = polyfit(log(xi_nrf), log(yi_boot), 1);
        boot_slopes(ib) = pb(1);
    end
    ci_lo = prctile(boot_slopes, 2.5); ci_hi = prctile(boot_slopes, 97.5);
    fprintf('P2b bootstrap 95%% interval on slope (2000 resamples of each cell''s draws, seed 20260922): [%.6f, %.6f]\n', ci_lo, ci_hi);
    p2b_ok = (slope_hat >= -0.75) && (slope_hat <= -0.45);
    fprintf('P2b PREDICTION (band [-0.75,-0.45], scored not gated): %s\n\n', loc_pf(p2b_ok));
end

if strcmp(stage, 'R4')
    loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
    diary off;
    return
end

% =========================================================================
%  R5 -- SELFCHK: W = I on the PRODUCTION grid; reproduces F-080 (Paper C
%        leg) and R1's closed form (GLOBECOM leg). Gates G-R5, G-C1.
% =========================================================================
if any(strcmp(stage, {'R5','ALL'}))
    fprintf('### R5 -- SELFCHK (W=I, production grid; G-R5, G-C1) ###\n\n');
    [GATE_ROWS, SUMMARY_ROWS, r5_pass] = loc_run_R5(P_PC, lam_c, d_ant, fc_const, ...
        SUMMARY_ROWS, GATE_ROWS, false);   %#ok<ASGLU>
    if ~r5_pass
        fprintf('HALT: G-R5 FAIL. Stop -- see F-157 note: if this fails but\n');
        fprintf('R5c passes, escalate (harness right, production reference wrong).\n');
        loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
        diary off;
        return
    end
end

if strcmp(stage, 'R5')
    loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
    diary off;
    return
end

% =========================================================================
%  R5c -- CENTRED SELFCHK (Addendum B.3 item 4, new, deterministic, harness
%         check): R5 on the CENTRED grid k_idx + 1/2. Must reproduce
%         F-080/F-081 to 1e-8 dB.
% =========================================================================
if any(strcmp(stage, {'R5c','ALL'}))
    fprintf('### R5c -- CENTRED GRID SELFCHK (Addendum B.3.4; must reproduce F-080/F-081) ###\n\n');
    [GATE_ROWS, SUMMARY_ROWS, r5c_pass] = loc_run_R5(P_PC, lam_c, d_ant, fc_const, ...
        SUMMARY_ROWS, GATE_ROWS, true);   %#ok<ASGLU>
    if ~r5c_pass
        fprintf('HALT: R5c FAIL. The harness does not reproduce T-49''s grid. Stop.\n');
        loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
        diary off;
        return
    end
end

if strcmp(stage, 'R5c')
    loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
    diary off;
    return
end

% =========================================================================
%  R6 -- SNR ladder, Paper C configuration, reported never gated
% =========================================================================
if any(strcmp(stage, {'R6','ALL'}))
    fprintf('### R6 -- SNR LADDER (Paper C configuration, reported never gated) ###\n\n');
    SNR6_list = [-5 0 10 20 25];
    theta6_deg = 40; r6 = 2.6578125; N_draws6 = 385; B6 = 400e6;
    for i = 1:numel(SNR6_list)
        CELL_ID = CELL_ID + 1;
        [SUMMARY_ROWS, DRAW_ROWS, SEED_REGISTRY, GATE_ROWS] = loc_run_cell( ...
            'PaperC', CELL_ID, 'R6', P_PC, B6, theta6_deg, r6, SNR6_list(i), N_draws6, ...
            lam_c, d_ant, fc_const, SEED_REGISTRY, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, false);
    end
    fprintf('R6 complete: %d SNR points.\n\n', numel(SNR6_list));
end

if strcmp(stage, 'R6')
    loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
    diary off;
    return
end

% =========================================================================
%  R7 -- SEEDCHK (Addendum B.7, new, reported NEVER gated). P-RUN
%        confirmation of F-151's seeding test. NOT a gate; NO verdict is
%        formed. VERBATIM combiner draw and CTL parameters from
%        t50_a4_bound.m loc_W_fixed / loc_vars.
% =========================================================================
if any(strcmp(stage, {'R7','ALL'}))
    fprintf('### R7 -- SEEDCHK (Addendum B.7; reported never gated; ~3 min) ###\n\n');

    Q = setup_production_P_v4('snr', 'full');
    Q.r_hi_fac = 0.20;                                          % t50_a4_bound.m L37-41
    D_ap  = (Q.M - 1) * Q.d_ant;
    Q.r_RD = 2 * D_ap^2 / Q.lambda_c;
    Q.u_min = 1/(Q.r_hi_fac*Q.r_RD*Q.u_margin);
    Q.u_max = 1/(Q.r_lo_fac*Q.r_RD/Q.u_margin);

    theta7_deg = 40; r7 = 2.6578125; SNR7_dB = 10;
    theta7 = theta7_deg * pi/180;
    N0_7 = 1/(10^(SNR7_dB/10));

    fprintf('R7 CTL config: theta=%g deg r=%g m r_hi_fac=%.2f SNR=%g dB M=%d N_RF=%d K_s=%d\n', ...
        theta7_deg, r7, Q.r_hi_fac, SNR7_dB, Q.M, Q.N_RF, Q.K_s);

    S_LIST = 0:5049;
    Vr_fj = nan(numel(S_LIST), 1);
    for is = 1:numel(S_LIST)
        s0v = S_LIST(is);
        rng(s0v, 'twister');                                          % VERBATIM t50_a4_bound.m loc_W_fixed
        W = (1/sqrt(Q.M)) * exp(1j * 2*pi * rand(Q.M, Q.N_RF));        % VERBATIM

        vv = loc_vars_local(theta7, r7, N0_7, Q, {W}, lam_c, d_ant);   % local reimplementation of loc_vars, eq path only
        Vr_fj(is) = vv;   % full-Jacobian equilibrated range variance
    end

    T7 = table((S_LIST.'), Vr_fj, 'VariableNames', {'s','V_r_fj'});
    r7csv = fullfile(OUT_DIR, sprintf('t53_R7_seedchk_%s.csv', ts));
    writetable(T7, r7csv);

    mean_window0 = mean(Vr_fj(1:50));            % s = 0..49
    mean_rest    = mean(Vr_fj(51:end));           % s = 50..5049
    delta_dB = 10*log10(mean_window0 / mean_rest);

    n_windows = 100;
    win_means = nan(n_windows,1);
    for iw = 1:n_windows
        idx = (51 + (iw-1)*50) : (50 + iw*50);
        win_means(iw) = mean(Vr_fj(idx));
    end
    sd_windows_dB = std(10*log10(win_means / mean(win_means)));

    all_window_means = [mean_window0; win_means];   % window 0 first, then 1..100
    [~, rank_order] = sort(all_window_means, 'descend');
    rank_window0 = find(rank_order == 1);

    fprintf('\nR7 RESULT: mean(s=0..49) vs mean(s=50..5049) = %+.6f dB\n', delta_dB);
    fprintf('R7 RESULT: SD of 100 disjoint 50-seed window means (s=50..5049) = %.6f dB\n', sd_windows_dB);
    fprintf('R7 RESULT: rank of window 0 among all 101 windows = %d of 101\n', rank_window0);
    fprintf('R7 PASS-THROUGH EXPECTATION (F-151, P-CALC): +0.419 dB, 0.224 dB, rank 4 of 101\n');
    fprintf('R7 is NOT a gate. No verdict is formed from these numbers.\n\n');
end

if strcmp(stage, 'R7')
    loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
    diary off;
    return
end

% =========================================================================
%  ALL -- final write, all gates re-summarised
% =========================================================================
if strcmp(stage, 'ALL')
    loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS);
    diary off;
end

end   % t53_randw_residual


% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function s = loc_pf(tf)
if tf, s = 'PASS'; else, s = 'FAIL'; end
end

function fp = loc_find_file(name)
fp = which(name);
if isempty(fp) && exist(name, 'file')
    fp = name;
end
end

function [al_sub, Ks] = loc_subcarrier_grid(B_val, fc, Delta_f)
% VERBATIM from stageA_s0_audit.m / t49_f032_fulljac_invariance.m /
% wb_crb_globecom2026.m. Do not edit. Used ONLY for the GLOBECOM (R1) leg.
    K_val = max(1, round(B_val / Delta_f));
    Ks    = min(512, K_val);
    if K_val == 1
        al_sub = 1.0;
    else
        k_v    = (0:K_val-1).';
        f_v    = fc + (k_v - (K_val-1)/2) * Delta_f;
        al_v   = f_v / fc;
        idx_s  = round(linspace(1, K_val, Ks));
        al_sub = al_v(idx_s);
    end
end

function [al_sub, Ks] = loc_prodgrid_alpha(B_val, fc, K_s_fixed) %#ok<INUSD>
% Paper C PRODUCTION grid (Addendum B.3): k_idx = (-(K_s/2):(K_s/2-1)).',
% alpha_k = 1 + k_idx * Delta_f / fc (setup_production_P_v4.m L161-163).
% Delta_f is fixed at the array's native 25 MHz (Sec. 6 of that file,
% Lesson L29); K_s is re-derived per bandwidth as K_s(B) = round(B /
% Delta_f_native). This is confirmed, not invented: Spec Sec. 3.10's R3
% cost row states "385 draws x 75 cells, K_s from 4 to 32" for B in
% {100,...,800} MHz, which is exactly round(B/25e6) at the endpoints
% (100/25=4, 800/25=32). K_s_fixed is accepted for interface symmetry
% with callers that already know K_s, but is not used: it is always
% re-derived from B_val here, consistent with every leg.
    Delta_f_native = 25e6;
    K_val = max(1, round(B_val / Delta_f_native));
    Ks = K_val;
    k_idx = (-(Ks/2):(Ks/2-1)).';
    al_sub = 1 + k_idx * (Delta_f_native / fc);
end

function [B_fj, B_diag, B_theta, raw] = loc_B_of_gW(theta, r, p, N0, M, N_RF, N_snap, ...
    alpha_sub, lam_c, d_ant, W)
% Sec. 3.1: form V_r(g,W,alpha=1) [narrowband replicate] and
% V_r(g,W,alpha=alpha_k) [wideband], both from the SAME W and g, via the
% live wb_crb_compressed.m T-50 interface (opts.W_list, opts.return_raw).
% NO edit is made to wb_crb_compressed.m.
K_s = numel(alpha_sub);

P = struct();
P.M = M; P.N = N_snap; P.N_RF = N_RF;
P.lambda_c = lam_c; P.d_ant = d_ant;
P.K_s = K_s;

opts = struct(); opts.W_list = {W}; opts.return_raw = true;

P_nb = P; P_nb.alpha_k_vec = ones(K_s, 1);
P_wb = P; P_wb.alpha_k_vec = alpha_sub(:);

[~, ~, info_nb] = wb_crb_compressed(theta, r, p, N0, P_nb, opts);
[~, ~, info_wb] = wb_crb_compressed(theta, r, p, N0, P_wb, opts);

% g_omega, g_kappa: Sec. 2.1 full-Jacobian range form
q_scale = 2*pi*d_ant/lam_c;
c_ch    = pi*d_ant^2/lam_c * sin(theta)^2;
kappa_v = c_ch / r;
domega_dtheta = -q_scale * sin(theta);
g_om = 2*r*cot(theta) / domega_dtheta;
g_ka = -r / kappa_v;
dkdr = c_ch / r^2;   % |d kappa / d r|

C_oo_nb = info_nb.raw.eq.C_oo(1); C_kk_nb = info_nb.raw.eq.C_kk(1); C_ok_nb = info_nb.raw.eq.C_ok(1);
C_oo_wb = info_wb.raw.eq.C_oo(1); C_kk_wb = info_wb.raw.eq.C_kk(1); C_ok_wb = info_wb.raw.eq.C_ok(1);

V_r_nb_fj = g_om^2*C_oo_nb + 2*g_om*g_ka*C_ok_nb + g_ka^2*C_kk_nb;
V_r_wb_fj = g_om^2*C_oo_wb + 2*g_om*g_ka*C_ok_wb + g_ka^2*C_kk_wb;
V_r_nb_diag = max(C_kk_nb, 0) / dkdr^2;
V_r_wb_diag = max(C_kk_wb, 0) / dkdr^2;
V_theta_nb  = max(C_oo_nb, 0) / (2*pi*d_ant/lam_c*sin(theta))^2 * (180/pi)^2;
V_theta_wb  = max(C_oo_wb, 0) / (2*pi*d_ant/lam_c*sin(theta))^2 * (180/pi)^2;

B_fj    = 10*log10(V_r_nb_fj / V_r_wb_fj);
B_diag  = 10*log10(V_r_nb_diag / V_r_wb_diag);
B_theta = 10*log10(V_theta_nb / V_theta_wb);

raw.V_r_nb = V_r_nb_fj; raw.V_r_wb = V_r_wb_fj;
raw.V_r_nb_diag = V_r_nb_diag; raw.V_r_wb_diag = V_r_wb_diag;
raw.C_oo_nb = C_oo_nb; raw.C_kk_nb = C_kk_nb; raw.C_ok_nb = C_ok_nb;
raw.C_oo_wb = C_oo_wb; raw.C_kk_wb = C_kk_wb; raw.C_ok_wb = C_ok_wb;
raw.cond_eq_nb = info_nb.raw.cond_eq_all(1); raw.cond_eq_wb = info_wb.raw.cond_eq_all(1);
raw.nsing_nb = info_nb.n_singular_all(1); raw.nsing_wb = info_wb.n_singular_all(1);
raw.V_theta_nb = V_theta_nb; raw.V_theta_wb = V_theta_wb;
end

function [SUMMARY_ROWS, DRAW_ROWS, SEED_REGISTRY, GATE_ROWS, Bfj_vec, sd_out, mean_out, maxcond_out, maxnsing_out] = ...
    loc_run_cell_nrf(config, cell_id, leg, P, B_hz, theta_deg, r_m, SNR_dB, N_draws, ...
    lam_c, d_ant, fc_const, SEED_REGISTRY, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS)
% Wraps loc_run_cell and additionally returns the cell's full B_fj draw
% vector (for P2b's draw-level bootstrap) plus its SD/mean/max-cond/
% max-nsing summaries, for R4's P2b fit and R4-SQ check.
[SUMMARY_ROWS, DRAW_ROWS, SEED_REGISTRY, GATE_ROWS, Bfj_vec, cond_eq_all, nsing_all] = ...
    loc_run_cell(config, cell_id, leg, P, B_hz, theta_deg, r_m, SNR_dB, N_draws, ...
    lam_c, d_ant, fc_const, SEED_REGISTRY, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, false); %#ok<ASGLU>
sd_out = std(Bfj_vec);
mean_out = mean(Bfj_vec);
maxcond_out = max(cond_eq_all);
maxnsing_out = max(nsing_all);
end

function [SUMMARY_ROWS, DRAW_ROWS, SEED_REGISTRY, GATE_ROWS, Bfj_vec, cond_eq_all_out, nsing_all_out] = ...
    loc_run_cell(config, cell_id, leg, P, B_hz, theta_deg, r_m, SNR_dB, N_draws, ...
    lam_c, d_ant, fc_const, SEED_REGISTRY, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, print_config)
% Runs one (leg,cell) at the Paper C array on the amended seed rule
% (Addendum B.1: 5300000 + 10000*cell_id + n), production grid
% (Addendum B.3), retaining G-CND/G-NS-flagged draws (Addendum B.4).
theta_r = theta_deg * pi/180;
N0 = 10^(-SNR_dB/10);
[al_sub, K_s] = loc_prodgrid_alpha(B_hz, fc_const, NaN);
mean_delta1 = mean(al_sub(:) - 1);
mean_delta2 = mean((al_sub(:) - 1).^2);

% PC-T2 requires the full config block printed for every leg (Sec. 3.7);
% print_config is accepted for interface symmetry with call sites but the
% config block is unconditional, per PC-T2's own requirement.
fprintf('PC-T2 CONFIG [%s cell %d]: M=%d N_RF=%d Delta_f=%.6g K_s=%d mean_delta1=%.9e mean_delta2=%.9e theta=%g deg r=%.7g m SNR=%g dB N=%d p=1 N0=%.6g seedbase=5300000+10000*%d draws=%d\n', ...
    leg, cell_id, P.M, P.N_RF, B_hz/K_s, K_s, mean_delta1, mean_delta2, theta_deg, r_m, SNR_dB, P.N, N0, cell_id, N_draws);

Bfj_vec = nan(N_draws,1); Bdiag_vec = nan(N_draws,1); Btheta_vec = nan(N_draws,1);
Vr_nb_vec = nan(N_draws,1); Vr_wb_vec = nan(N_draws,1);   % for EW_V_domain_dB (Sec. 3.5)
cond_eq_nb_v = nan(N_draws,1); cond_eq_wb_v = nan(N_draws,1);
nsing_nb_v = nan(N_draws,1); nsing_wb_v = nan(N_draws,1);
flag_cnd_v = false(N_draws,1); flag_ns_v = false(N_draws,1);
seeds_used = nan(N_draws,1);

COND_EQ_LIMIT = 4.5e9;   % G-CND, A4 G9 threshold

for n = 1:N_draws
    seed_n = 5300000 + 10000*cell_id + n;
    seeds_used(n) = seed_n;
    s = RandStream('mt19937ar', 'Seed', seed_n);
    Wn = (1/sqrt(P.M)) * exp(1j * 2*pi * rand(s, P.M, P.N_RF));

    [Bfj, Bdiag, Btheta, raw] = loc_B_of_gW(theta_r, r_m, 1, N0, P.M, P.N_RF, P.N, ...
        al_sub, lam_c, d_ant, Wn);

    Bfj_vec(n) = Bfj; Bdiag_vec(n) = Bdiag; Btheta_vec(n) = Btheta;
    Vr_nb_vec(n) = raw.V_r_nb; Vr_wb_vec(n) = raw.V_r_wb;
    cond_eq_nb_v(n) = raw.cond_eq_nb; cond_eq_wb_v(n) = raw.cond_eq_wb;
    nsing_nb_v(n) = raw.nsing_nb; nsing_wb_v(n) = raw.nsing_wb;
    flag_cnd_v(n) = (raw.cond_eq_nb > COND_EQ_LIMIT) || (raw.cond_eq_wb > COND_EQ_LIMIT);
    flag_ns_v(n)  = (raw.nsing_nb ~= 0) || (raw.nsing_wb ~= 0);

    DRAW_ROWS(end+1,:) = {config, cell_id, leg, P.M, P.N_RF, B_hz/K_s, B_hz, K_s, mean_delta2, ...
        theta_deg, r_m, SNR_dB, n, seed_n, Bfj, Bdiag, Btheta, raw.V_r_nb, raw.V_r_wb, ...
        raw.V_r_nb_diag, raw.V_r_wb_diag, raw.C_oo_nb, raw.C_kk_nb, raw.C_ok_nb, ...
        raw.C_oo_wb, raw.C_kk_wb, raw.C_ok_wb, raw.cond_eq_nb, raw.cond_eq_wb, ...
        raw.nsing_nb, raw.nsing_wb}; %#ok<AGROW,SAGROW>
end

% G-SEED: assert no duplicate seed across the whole run (Addendum B.1)
dup_here = numel(seeds_used) ~= numel(unique(seeds_used));
dup_global = false;
if ~isempty(SEED_REGISTRY)
    combined = [SEED_REGISTRY(:); seeds_used(:)];
    dup_global = numel(combined) ~= numel(unique(combined));
end
gseed_ok = ~dup_here && ~dup_global;
fprintf('G-SEED [%s cell %d]: %s\n', leg, cell_id, loc_pf(gseed_ok));
GATE_ROWS(end+1,:) = {'G-SEED', leg, cell_id, NaN, NaN, loc_pf(gseed_ok)}; %#ok<SAGROW>
if ~gseed_ok
    error('t53_randw_residual: G-SEED FAIL at %s cell %d. Seed allocation is wrong. Halting.', leg, cell_id);
end
SEED_REGISTRY = [SEED_REGISTRY(:); seeds_used(:)];

% G-CND, G-NS (non-halting; retained + flagged, Addendum B.4)
n_flag_cnd = sum(flag_cnd_v); n_flag_ns = sum(flag_ns_v);
fprintf('G-CND [%s cell %d]: max cond_eq = %.6e  (limit %.3e)  flagged=%d/%d\n', ...
    leg, cell_id, max([cond_eq_nb_v; cond_eq_wb_v]), COND_EQ_LIMIT, n_flag_cnd, N_draws);
fprintf('G-NS  [%s cell %d]: flagged=%d/%d\n', leg, cell_id, n_flag_ns, N_draws);
GATE_ROWS(end+1,:) = {'G-CND', leg, cell_id, max([cond_eq_nb_v; cond_eq_wb_v]), COND_EQ_LIMIT, loc_pf(n_flag_cnd==0)}; %#ok<SAGROW>
GATE_ROWS(end+1,:) = {'G-NS', leg, cell_id, n_flag_ns, 0, loc_pf(n_flag_ns==0)}; %#ok<SAGROW>

cond_eq_all_out = [cond_eq_nb_v; cond_eq_wb_v];
nsing_all_out = [nsing_nb_v; nsing_wb_v];

% ---- per-cell report, primary (retained) sample + unflagged subset -------
SUMMARY_ROWS = loc_report_cell_all(SUMMARY_ROWS, config, cell_id, leg, P.M, P.N_RF, B_hz, ...
    K_s, theta_deg, r_m, SNR_dB, N_draws, Bfj_vec, Bdiag_vec, Btheta_vec, ...
    max(cond_eq_nb_v), max(cond_eq_wb_v), mean_delta2, Vr_nb_vec, Vr_wb_vec, ...
    max(max(nsing_nb_v), max(nsing_wb_v)));

any_flag = flag_cnd_v | flag_ns_v;
if any(any_flag)
    keep = ~any_flag;
    fprintf('  UNFLAGGED SUBSET [%s cell %d]: excluded %d of %d draws\n', leg, cell_id, sum(any_flag), N_draws);
    SUMMARY_ROWS = loc_report_cell_all(SUMMARY_ROWS, [config '_UNFLAGGED'], cell_id, leg, P.M, P.N_RF, B_hz, ...
        K_s, theta_deg, r_m, SNR_dB, sum(keep), Bfj_vec(keep), Bdiag_vec(keep), Btheta_vec(keep), ...
        max(cond_eq_nb_v(keep)), max(cond_eq_wb_v(keep)), mean_delta2, Vr_nb_vec(keep), Vr_wb_vec(keep), ...
        max(max(nsing_nb_v(keep)), max(nsing_wb_v(keep))));
end
end

function SUMMARY_ROWS = loc_report_cell_all(SUMMARY_ROWS, config, cell_id, leg, M, N_RF, B_hz, ...
    K_s, theta_deg, r_m, SNR_dB, n, Bfj, Bdiag, Btheta, max_cond_eq_nb, max_cond_eq_wb, ...
    mean_delta2, Vr_nb, Vr_wb, max_nsing)
if nargin < 21 || isempty(max_nsing), max_nsing = NaN; end
% Per-cell reporting, F-030's standard PLUS Addendum B.5 (Wilson intervals,
% distribution-free percentile intervals, bootstrap halfwidth). The
% residual is NEVER printed without SD + 2.5/97.5 percentiles (G-DST).
%
% closed_grid_dB (Sec. 3.5, Addendum B.3 item 3, renamed B_ref_centred_approx):
%   10*log10(1 + mean_delta2), reported never gated (the spec's original
%   "closed form" column; NOT B_ref -- B_ref is computed and gated
%   separately at R5/R5c via loc_run_R5, per Addendum B.3 item 2).
% closed_F081_dB (Sec. 3.5, F-081 exact discrete form):
%   10*log10(1 + (B/fc)^2 * (1 - 1/K_s^2)/12), native-grid form.
quantities = {'B_fj', Bfj; 'B_diag', Bdiag; 'B_theta', Btheta; 'B_fj_minus_diag', Bfj - Bdiag};

if isempty(K_s) || isnan(K_s) || K_s == 0
    closed_grid_dB = NaN;
    closed_F081_dB = NaN;
else
    closed_grid_dB = 10*log10(1 + mean_delta2);
    fc_local = 28e9;
    closed_F081_dB = 10*log10(1 + (B_hz/fc_local)^2 * (1 - 1/K_s^2) / 12);
end

% EW_V_domain_dB (Sec. 3.5): the E_W[V]-domain aggregate,
% 10*log10(mean_n V_r_nb / mean_n V_r_wb), reported never gated.
Vr_nb = Vr_nb(:); Vr_wb = Vr_wb(:);
Vr_nb = Vr_nb(~isnan(Vr_nb)); Vr_wb = Vr_wb(~isnan(Vr_wb));
if isempty(Vr_nb) || isempty(Vr_wb)
    EW_V_domain_dB = NaN;
else
    EW_V_domain_dB = 10*log10(mean(Vr_nb) / mean(Vr_wb));
end

max_cond_eq = max(max_cond_eq_nb, max_cond_eq_wb);

for iq = 1:size(quantities,1)
    qname = quantities{iq,1};
    v = quantities{iq,2};
    v = v(:);
    v = v(~isnan(v));
    if isempty(v)
        continue
    end
    mean_dB = mean(v); median_dB = median(v); sd_dB = std(v);
    se_dB = sd_dB / sqrt(numel(v));
    frac_pos = mean(v > 0);
    pct2p5 = prctile(v, 2.5); pct97p5 = prctile(v, 97.5);

    % Wilson 95% interval on P(B>0) and P(B<0) (Addendum B.5)
    zc = 1.959963984540054;
    nn = numel(v);
    phat_pos = frac_pos;
    wilson_pos = loc_wilson(phat_pos, nn, zc);
    phat_neg = mean(v < 0);
    wilson_neg = loc_wilson(phat_neg, nn, zc);

    % Distribution-free (binomial order-statistic) 95% interval on the
    % 2.5th and 97.5th percentiles (Addendum B.5)
    [lo_idx_025, hi_idx_025] = loc_order_stat_ci(nn, 0.025, zc);
    [lo_idx_975, hi_idx_975] = loc_order_stat_ci(nn, 0.975, zc);
    vs = sort(v);
    df_pct2p5_lo = vs(max(lo_idx_025,1)); df_pct2p5_hi = vs(min(hi_idx_025,nn));
    df_pct97p5_lo = vs(max(lo_idx_975,1)); df_pct97p5_hi = vs(min(hi_idx_975,nn));

    % Bootstrap 95% halfwidth on the mean, 2000 resamples (Addendum B.5)
    NBOOT = 2000; rng(20260922 + cell_id + iq, 'twister');
    bm = nan(NBOOT,1);
    for ib = 1:NBOOT
        idx = randi(nn, 1, nn);
        bm(ib) = mean(v(idx));
    end
    boot_ci = prctile(bm, [2.5 97.5]);
    achieved_halfwidth = (boot_ci(2) - boot_ci(1)) / 2;

    fprintf('  [%s] cell=%d leg=%s n=%d | mean=%+.6f SD=%.6f median=%+.6f fracpos=%.4f p2.5=%+.6f p97.5=%+.6f dB\n', ...
        qname, cell_id, leg, nn, mean_dB, sd_dB, median_dB, frac_pos, pct2p5, pct97p5);
    fprintf('    Wilson95 P(B>0)=[%.4f,%.4f]  P(B<0)=[%.4f,%.4f]  bootstrap halfwidth(mean)=%.6f dB (target 0.05)\n', ...
        wilson_pos(1), wilson_pos(2), wilson_neg(1), wilson_neg(2), achieved_halfwidth);
    fprintf('    distribution-free 95%% CI on p2.5=[%+.6f,%+.6f]  on p97.5=[%+.6f,%+.6f]\n', ...
        df_pct2p5_lo, df_pct2p5_hi, df_pct97p5_lo, df_pct97p5_hi);

    SUMMARY_ROWS(end+1,:) = {config, cell_id, leg, M, N_RF, B_hz, K_s, theta_deg, r_m, SNR_dB, nn, ...
        qname, mean_dB, median_dB, sd_dB, se_dB, frac_pos, pct2p5, pct97p5, EW_V_domain_dB, ...
        closed_grid_dB, closed_F081_dB, max_cond_eq, max_nsing}; %#ok<AGROW,SAGROW>
end
end

function w = loc_wilson(phat, n, z)
denom = 1 + z^2/n;
centre = phat + z^2/(2*n);
half = z * sqrt(phat*(1-phat)/n + z^2/(4*n^2));
w = [(centre - half)/denom, (centre + half)/denom];
w = max(min(w,1),0);
end

function [lo, hi] = loc_order_stat_ci(n, q, z)
% Binomial order-statistic 95% CI on the q-th percentile's rank.
mu = n*q; sg = sqrt(n*q*(1-q));
lo = floor(mu - z*sg);
hi = ceil(mu + z*sg) + 1;
end

function [GATE_ROWS, SUMMARY_ROWS, pass_all] = loc_run_R5(P_PC, lam_c, d_ant, fc_const, SUMMARY_ROWS, GATE_ROWS, centred)
% R5 / R5c: W = I at every R3 (Paper C) cell and every R1 (GLOBECOM)
% bandwidth. Deterministic (one draw per cell). G-R5 compares B_fj with
% B_ref = 10log10(mean_k alpha_k^2) on the grid ACTUALLY USED
% (Addendum B.3 item 2) -- NOT with F-081's closed form directly, except
% R5c which is defined to reproduce F-080/F-081 on the CENTRED grid.
R3_B = [100 200 400 600 800]*1e6;
R3_R = [1.50 2.13 3.00 5.00 10.00];
R3_TH = [20 40 60];

M_pc = P_PC.M; N_RF_pc = P_PC.N_RF; SNR_dB = 10; N0 = 10^(-SNR_dB/10);
W_eye = eye(M_pc);

devs = [];
c1_ratios = [];
n_cells = 0;
for iB = 1:numel(R3_B)
    Delta_f_native = 25e6;
    K_val = max(1, round(R3_B(iB) / Delta_f_native));
    k_idx = (-(K_val/2):(K_val/2-1)).';
    if centred
        k_idx = k_idx + 0.5;
    end
    al_sub = 1 + k_idx * (Delta_f_native / fc_const);
    B_ref = 10*log10(mean(al_sub.^2));

    for iT = 1:numel(R3_TH)
        for iR = 1:numel(R3_R)
            n_cells = n_cells + 1;
            theta_r = R3_TH(iT)*pi/180;
            rr = R3_R(iR);
            [Bfj, Bdiag, Btheta, raw] = loc_B_of_gW(theta_r, rr, 1, N0, M_pc, N_RF_pc, P_PC.N, ...
                al_sub, lam_c, d_ant, W_eye); %#ok<ASGLU>
            dev = abs(Bfj - B_ref);
            devs(end+1) = dev; %#ok<AGROW>

            % G-C1 (Sec. 3.7): cross term vanishes at W=I.
            % abs(2*g_omega*g_kappa*C_ok) / (g_kappa^2*C_kk), threshold
            % 1e-12 (T-49's C1 threshold; F-080 measured max 1.377e-16).
            q_scale = 2*pi*d_ant/lam_c;
            c_ch = pi*d_ant^2/lam_c*sin(theta_r)^2;
            kap = c_ch / rr;
            domega_dtheta = -q_scale*sin(theta_r);
            g_om = 2*rr*cot(theta_r) / domega_dtheta;
            g_ka = -rr / kap;
            c1_ratio = abs(2*g_om*g_ka*raw.C_ok_wb) / max(g_ka^2*abs(raw.C_kk_wb), realmin);
            c1_ratios(end+1) = c1_ratio; %#ok<AGROW>
        end
    end
end

if centred
    max_dev = max(devs);
    r5c_ok = (max_dev <= 1e-8);
    fprintf('R5c: %d cells, max |B_fj - B_ref(centred)| = %.6e dB (<=1e-8? %s)\n', n_cells, max_dev, loc_pf(r5c_ok));
    GATE_ROWS(end+1,:) = {'R5c_reproF080F081','R5c',sprintf('ALL_%dcells',n_cells),max_dev,1e-8,loc_pf(r5c_ok)}; %#ok<SAGROW>
    pass_all = r5c_ok;
else
    max_dev = max(devs);
    gr5_ok = (max_dev <= 1e-8);
    fprintf('G-R5: %d cells, max |B_fj - B_ref| = %.6e dB (<=1e-8? %s)  [HARNESS SELF-CHECK, NOT A PHYSICS TEST]\n', ...
        n_cells, max_dev, loc_pf(gr5_ok));
    GATE_ROWS(end+1,:) = {'G-R5','R5',sprintf('ALL_%dcells',n_cells),max_dev,1e-8,loc_pf(gr5_ok)}; %#ok<SAGROW>

    % G-C1 (cross term at W=I; threshold 1e-12, F-080 measured max 1.377e-16)
    max_c1 = max(c1_ratios);
    gc1_ok = (max_c1 <= 1e-12);
    fprintf('G-C1: %d cells, max cross-term ratio = %.6e (<=1e-12? %s)\n', n_cells, max_c1, loc_pf(gc1_ok));
    GATE_ROWS(end+1,:) = {'G-C1','R5',sprintf('ALL_%dcells',n_cells),max_c1,1e-12,loc_pf(gc1_ok)}; %#ok<SAGROW>
    pass_all = gr5_ok;
end
end

function vv = loc_vars_local(theta, r, N0, Q, Wl, lam_c, d_ant)
% Local, minimal reimplementation of t50_a4_bound.m's loc_vars, EQUILIBRATED
% full-Jacobian range variance ONLY (what R7 needs), via the live
% wb_crb_compressed.m interface (opts.W_list, opts.return_raw). Averaged
% over the combiners in Wl (here always {W}, one seed at a time).
o = struct(); o.W_list = Wl; o.return_raw = true;
[~, ~, info] = wb_crb_compressed(theta, r, 1, N0, Q, o);
q_scale = 2*pi*d_ant/lam_c;
c_ch = pi*d_ant^2/lam_c*sin(theta)^2;
kap = c_ch / r;
g_om = -2*r*cos(theta) / (q_scale*sin(theta)^2);
g_ka = -r/kap;
Coo = info.raw.eq.C_oo; Ckk = info.raw.eq.C_kk; Cok = info.raw.eq.C_ok;
vr_fj = g_om^2*Coo + 2*g_om*g_ka*Cok + g_ka^2*Ckk;
vv = mean(vr_fj);
end

function loc_finish(OUT_DIR, ts, t_start, DRAW_ROWS, SUMMARY_ROWS, GATE_ROWS, REPRO_ROWS)
% Writes the (up to) five Sec. 3.6 output CSVs, ASCII, %.12g or better on
% every dB column, and prints the gate summary block.
draw_cols = {'leg','cell_id','config','M','N_RF','Delta_f_hz','B_hz','K_s','mean_delta2', ...
    'theta_deg','r_m','SNR_dB','draw','seed','B_fj_dB','B_diag_dB','B_theta_dB', ...
    'V_r_nb_fj','V_r_wb_fj','V_r_nb_diag','V_r_wb_diag','C_oo_nb','C_kk_nb','C_ok_nb', ...
    'C_oo_wb','C_kk_wb','C_ok_wb','cond_eq_nb','cond_eq_wb','nsing_nb','nsing_wb'};
summ_cols = {'leg','cell_id','config','M','N_RF','B_hz','K_s','theta_deg','r_m','SNR_dB','n', ...
    'quantity','mean_dB','median_dB','sd_dB','se_dB','frac_positive','pct2p5_dB','pct97p5_dB', ...
    'EW_V_domain_dB','closed_grid_dB','closed_F081_dB','max_cond_eq','max_nsing'};
gate_cols = {'gate','leg','cell_id','statistic','threshold','verdict'};
repro_cols = {'draw','seed','B_MHz','G_exact_stageA','B_diag_t53','reldiff_exact', ...
    'G_fullJac_stageA','B_fj_t53','reldiff_fj','absdiff_exact','absdiff_fj'};

if ~isempty(DRAW_ROWS)
    % Rows were appended as [config,cell_id,leg,...] (Sec. 3.1's natural
    % construction order); draw_cols' header is [leg,cell_id,config,...]
    % (Sec. 3.6 item 2). Reorder once here, same swap as SUMMARY_ROWS below.
    D = DRAW_ROWS;
    D2 = cell(size(D,1), numel(draw_cols));
    for i = 1:size(D,1)
        D2(i,:) = [D(i,3), D(i,2), D(i,1), D(i,4:end)];
    end
    Td = cell2table(D2, 'VariableNames', draw_cols);
    writetable(Td, fullfile(OUT_DIR, sprintf('t53_draws_%s.csv', ts)));
end
if ~isempty(SUMMARY_ROWS)
    % Reorder to match [leg,cell_id,config,...] header above (rows were
    % appended as [config,cell_id,leg,...]; fix column order once here).
    S = SUMMARY_ROWS;
    S2 = cell(size(S,1), numel(summ_cols));
    for i = 1:size(S,1)
        S2(i,:) = {S{i,3}, S{i,2}, S{i,1}, S{i,4}, S{i,5}, S{i,6}, S{i,7}, S{i,8}, S{i,9}, ...
            S{i,10}, S{i,11}, S{i,12}, S{i,13}, S{i,14}, S{i,15}, S{i,16}, S{i,17}, S{i,18}, ...
            S{i,19}, S{i,20}, S{i,21}, S{i,22}, S{i,23}, S{i,24}};
    end
    Ts = cell2table(S2, 'VariableNames', summ_cols);
    writetable(Ts, fullfile(OUT_DIR, sprintf('t53_summary_%s.csv', ts)));
end
if ~isempty(GATE_ROWS)
    Tg = cell2table(GATE_ROWS, 'VariableNames', gate_cols);
    writetable(Tg, fullfile(OUT_DIR, sprintf('t53_gates_%s.csv', ts)));
end
if ~isempty(REPRO_ROWS)
    Tr = cell2table(REPRO_ROWS, 'VariableNames', repro_cols);
    writetable(Tr, fullfile(OUT_DIR, sprintf('t53_repro_%s.csv', ts)));
end

fprintf('\n=============================================================\n');
fprintf('  GATE SUMMARY\n');
fprintf('=============================================================\n');
if ~isempty(GATE_ROWS)
    for i = 1:size(GATE_ROWS,1)
        fprintf('  %-12s leg=%-10s cell=%-6s verdict=%s\n', ...
            GATE_ROWS{i,1}, num2str(GATE_ROWS{i,2}), num2str(GATE_ROWS{i,3}), GATE_ROWS{i,6});
    end
end
fprintf('\n  elapsed %.1f s\n', toc(t_start));
end
