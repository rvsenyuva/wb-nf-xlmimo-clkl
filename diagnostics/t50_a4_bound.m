%% t50_a4_bound.m -- A4: T-50 geometry-matched CRB and T-57 decision
%
%  Spec    : PaperC_B4_T50_T57_BoundSpec.md (Lane B session B4, 2026-09-17).
%  Purpose : bound deliverable only. Post-processing of the A3 per-trial
%            geometries through the CRB code. Runs NO estimator and NO
%            channel generator. Forms NO estimator statistic of any kind.
%  Requires: - PATCHED wb_crb_compressed.m (spec Sec. 3, edits E1-E6)
%            - wb_crb_compressed_prepatch.m = the fd301c4 file with ONLY its
%              line-1 function name changed to wb_crb_compressed_prepatch
%            - live setup_production_P_v4.m, nf_usw_steer.m,
%              wb_nf_fresnel_steer.m on the MATLAB path
%            - the A3 archive CSVs in A3_DIR
%  Run     : author-side, MATLAB R2025b, 6-worker pool. Before running,
%            paste `git rev-parse --short HEAD` output into the diary.
%  Author of script: Claude Opus 5 (not executed by Claude; gates G0-G3
%            are designed to catch interface errors before any T-50 number
%            is read).
%
%  Encoding: 7-bit ASCII.

%% ---- 0. Configuration -------------------------------------------------
DRY_RUN = false;    % set false for the real run
A3_DIR   = 'C:\Users\senyu\Downloads\paperC_A3_20260907\paperC_A3_20260907';
OUT_DIR  = fullfile('results', 'paperC_A4_t50');
CTL_TS   = '20260908_145855';                          % seed base 1000
FRESH_TS = {'20260908_161501', '20260908_164830', ...  % 3000000, 3010000
            '20260908_172101', '20260908_175319'};     % 3020000, 3030000
T52B_TS  = '20260908_160303';                          % pinned gate leg
RUN_USW  = true;   % leg U (reported only); auto-disabled if PC-3 fails

if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end
ts = datestr(now, 'yyyymmdd_HHMMSS');
diary(fullfile(OUT_DIR, sprintf('A4_console_%s.txt', ts)));
fprintf('A4 (T-50/T-57) started %s\n', datestr(now));

%% ---- 1. Parameter structs, exactly as A3_console.txt L5, L37-40, L146, L179
P = setup_production_P_v4('snr', 'full');
if ~isfield(P, 'r_RD'), P.r_RD = 2*((P.M-1)*P.d_ant)^2 / P.lambda_c; end
P.r_hi_fac = 0.20;                                     % F-100 call-site override
P.u_min = 1/(P.r_hi_fac*P.r_RD*P.u_margin);
P.u_max = 1/(P.r_lo_fac*P.r_RD/P.u_margin);
Pc = setup_production_P_v4('convergence', 'restricted');
if ~isfield(Pc, 'r_RD'), Pc.r_RD = 2*((Pc.M-1)*Pc.d_ant)^2 / Pc.lambda_c; end
fprintf('CHECK P : r_lo=%.4f r_hi=%.4f N_seed=%d rng_seed_W=%d N_RF=%d K_s=%d\n', ...
    P.r_lo_fac*P.r_RD, P.r_hi_fac*P.r_RD, P.N_seed, P.rng_seed_W, P.N_RF, P.K_s);
fprintf('CHECK Pc: r_hi_fac=%.2f conv_theta=%.2f deg conv_r=%.3f m\n', ...
    Pc.r_hi_fac, Pc.conv_theta*180/pi, Pc.conv_r);

%% ---- 2. PC-1: provenance (line counts) ---------------------------------
files = {'wb_crb_compressed.m', 'wb_crb_compressed_prepatch.m', ...
         'run_monte_carlo_paperC.m', 'setup_production_P_v4.m', ...
         'wb_channel_gen_ofdm_nf.m', 'nf_usw_steer.m', 'wb_nf_fresnel_steer.m'};
for i = 1:numel(files)
    fp = which(files{i});
    if isempty(fp)
        fprintf('PC-1 %-30s NOT ON PATH\n', files{i});
    else
        txt = fileread(fp);
        fprintf('PC-1 %-30s lines=%d  %s\n', files{i}, sum(txt == newline), fp);
    end
end
% expect: prepatch 374, run_monte_carlo_paperC 870, setup 341, generator 259

%% ---- 3. PC-3 / PC-4: steering-vector forms (leg U precondition) --------
c1 = 2*pi*P.d_ant/P.lambda_c;                 % omega = c1*cos(theta)
c0 = pi*P.d_ant^2/P.lambda_c;                 % kappa = c0*sin^2(theta)/r
mb = ((0:P.M-1).' - (P.M-1)/2);
forms = struct('name', {'centred,-', 'centred,+', 'first,-', 'first,+'}, ...
               'delta', {mb*P.d_ant, mb*P.d_ant, (0:P.M-1).'*P.d_ant, (0:P.M-1).'*P.d_ant}, ...
               'sgn', {-1, +1, -1, +1});
err_f = zeros(1, numel(forms)); err_fres = [0 0]; nrm_dev = 0;
pc_ok = true;
try
tv = [20 40 60]*pi/180; rv = [1.0631 2.6578125 4.2525];
av = [P.alpha_k_vec(1) 1 P.alpha_k_vec(end)];
for it = 1:3
  for ir = 1:3
    for ia = 1:3
      th = tv(it); r = rv(ir); al = av(ia); kk = 2*pi*al/P.lambda_c;
      Pl = P; Pl.lambda = P.lambda_c / al;               % generator L250
      a_code = nf_usw_steer(th, r, Pl); a_code = a_code(:);
      nrm_dev = max(nrm_dev, abs(norm(a_code)^2 - P.M));
      for f = 1:numel(forms)
        dl = forms(f).delta;
        dist = sqrt(r^2 + dl.^2 - 2*r*dl*cos(th));
        a_f = exp(forms(f).sgn*1j*kk*(dist - r));
        err_f(f) = max(err_f(f), max(abs(a_code - a_f)));
      end
      a_fr = wb_nf_fresnel_steer(th, 1/r, al, P) * sqrt(P.M); a_fr = a_fr(:);
      om = c1*cos(th); ka = c0*sin(th)^2/r;
      err_fres(1) = max(err_fres(1), max(abs(a_fr - exp(1j*al*(om*mb - ka*mb.^2)))));
      err_fres(2) = max(err_fres(2), max(abs(a_fr - exp(1j*al*(om*mb + ka*mb.^2)))));
    end
  end
end
catch ME
    pc_ok = false;
    fprintf('PC-3/PC-4 ERROR (leg U will be VOID): %s\n', ME.message);
end
for f = 1:numel(forms)
    fprintf('PC-3 nf_usw_steer vs exact-distance form [%s]: max abs diff = %.3e\n', forms(f).name, err_f(f));
end
fprintf('PC-3 nf_usw_steer max |norm(a)^2 - M| = %.3e\n', nrm_dev);
fprintf('PC-4 wb_nf_fresnel_steer*sqrt(M) vs CRB form exp(j a(w m - k m^2)): %.3e ; vs (+k): %.3e\n', ...
    err_fres(1), err_fres(2));
hit = find(err_f <= 1e-12);
usw.valid = RUN_USW && pc_ok && numel(hit) == 1;
if usw.valid
    usw.delta = forms(hit).delta; usw.sgn = forms(hit).sgn;
    fprintf('PC-3 RESULT: leg U ENABLED with form [%s]\n', forms(hit).name);
else
    usw.delta = []; usw.sgn = NaN;
    fprintf('PC-3 RESULT: leg U VOID (%d forms matched)\n', numel(hit));
end

%% ---- 4. G1 / G2 / G3 at the nominal geometry ---------------------------
legs(1) = struct('name', 'CTL',  'Q', P,  'csv', fullfile(A3_DIR, ['mc_snr_sweep_' CTL_TS '.csv']));
legs(2) = struct('name', 'T52B', 'Q', Pc, 'csv', fullfile(A3_DIR, ['mc_convergence_sweep_' T52B_TS '.csv']));
rep_rows = {};
n_cells = 0; n_G1a = 0; n_G1b = 0; n_G2 = 0; n_G3C = 0; n_G3F = 0;
for L = 1:numel(legs)
    Q = legs(L).Q;
    [snr_s, pth_s, pr_s] = loc_pub_crb_strings(legs(L).csv);
    [th0, r0] = loc_nominal(Q);
    Wl = loc_W_fixed(Q);
    fprintf('G1 leg %s: nominal theta=%.6f deg r=%.7f m, %d SNR rows\n', ...
        legs(L).name, th0*180/pi, r0, numel(snr_s));
    for j = 1:numel(snr_s)
        snr = str2double(snr_s(j)); N0 = 10^(-snr/10); p0 = ones(Q.d,1)/Q.d;
        [cr_d, ct_d, in_d] = wb_crb_compressed(repmat(th0,Q.d,1), repmat(r0,Q.d,1), p0, N0, Q);
        [cr_o, ct_o, in_o] = wb_crb_compressed_prepatch(repmat(th0,Q.d,1), repmat(r0,Q.d,1), p0, N0, Q);
        o = struct(); o.W_list = Wl;
        [cr_w, ct_w, in_w] = wb_crb_compressed(repmat(th0,Q.d,1), repmat(r0,Q.d,1), p0, N0, Q, o);
        o = struct(); o.W_list = {eye(Q.M)};
        [cr_f, ct_f] = wb_crb_compressed(repmat(th0,Q.d,1), repmat(r0,Q.d,1), p0, N0, Q, o);
        s_r = sprintf('%.6g', mean(cr_d)); s_t = sprintf('%.6g', mean(ct_d));
        f_r = sprintf('%.6g', mean(cr_f)); f_t = sprintf('%.6g', mean(ct_f));
        g1a = [strcmp(s_r, pr_s(j)), strcmp(s_t, pth_s(j))];
        g1b = isequal(cr_d, cr_o) && isequal(ct_d, ct_o) && ...
              isequal(in_d.crb_r_all, in_o.crb_r_all) && isequal(in_d.crb_theta_all, in_o.crb_theta_all);
        g2  = isequal(cr_d, cr_w) && isequal(ct_d, ct_w) && ...
              isequal(in_d.crb_r_all, in_w.crb_r_all) && isequal(in_d.crb_theta_all, in_w.crb_theta_all);
        g3F = [strcmp(f_r, pr_s(j)), strcmp(f_t, pth_s(j))];
        n_cells = n_cells + 2; n_G1a = n_G1a + sum(g1a); n_G1b = n_G1b + 2*g1b;
        n_G2 = n_G2 + 2*g2; n_G3C = n_G3C + sum(g1a); n_G3F = n_G3F + sum(g3F);
        rep_rows(end+1, :) = {legs(L).name, char(snr_s(j)), char(pr_s(j)), s_r, char(pth_s(j)), s_t, ...
                               f_r, f_t, g1a(1), g1a(2), double(g1b), double(g2), g3F(1), g3F(2)}; %#ok<SAGROW>
    end
end
fprintf('G1a (patched default reproduces published, %%.6g string-exact): %d / %d\n', n_G1a, n_cells);
fprintf('G1b (patched default bit-identical to pre-patch):              %d / %d\n', n_G1b, n_cells);
fprintf('G2  (W_list seeds 0..49 bit-identical to default):             %d / %d\n', n_G2, n_cells);
fprintf('G3  reading C (compressed, 50-seed) matches: %d / %d ; reading F (W = I) matches: %d / %d\n', ...
    n_G3C, n_cells, n_G3F, n_cells);
if n_G3C == n_cells && n_G3F == 0
    fprintf('G3 DECISION: COMPRESSED\n');
elseif n_G3F == n_cells && n_G3C == 0
    fprintf('G3 DECISION: FULL-ARRAY\n');
else
    fprintf('G3 DECISION: UNRESOLVED\n');
end
Trep = cell2table(rep_rows, 'VariableNames', {'leg','SNR_dB','pub_crb_r_m','rep_crb_r_m', ...
    'pub_crb_theta_deg','rep_crb_theta_deg','fullarray_crb_r_m','fullarray_crb_theta_deg', ...
    'G1a_r','G1a_theta','G1b','G2','G3F_r','G3F_theta'});
writetable(Trep, fullfile(OUT_DIR, sprintf('t50_a4_g0_reproduction_%s.csv', ts)));
if ~(n_G1a == n_cells && n_G1b == n_cells && n_G2 == n_cells && n_G3C == n_cells && n_G3F == 0)
    fprintf('HALT: G1/G2/G3 not all PASS. No T-50 quantity is computed (spec Sec. 6.1).\n');
    diary off; return
end

%% ---- 5. Randomised blocks: per-trial bounds ----------------------------
blk_ts = [{CTL_TS}, FRESH_TS];
SNR_VEC = -5:2.5:25;  S = numel(SNR_VEC);  S_READ = 13;
if DRY_RUN, blk_ts = blk_ts(1); end
Wl = loc_W_fixed(P);
[th0, r0] = loc_nominal(P);
colnames = loc_colnames();
NC = numel(colnames);
summ = {};  gate_rows = {};
Dspread = nan(numel(blk_ts), S);  FAfac = nan(numel(blk_ts), 1);
all_ok_G4 = true; all_ok_G5 = true; max_cond = 0; recon_total = 0; n_total = 0;
for b = 1:numel(blk_ts)
    tcsv = fullfile(A3_DIR, ['mc_snr_trials_' blk_ts{b} '.csv']);
    [seeds, base, th_str, r_str, g0b] = loc_read_trials(tcsv, S_READ);
    n = numel(seeds);
    if DRY_RUN, n = 60; seeds = seeds(1:n); th_str = th_str(1:n); r_str = r_str(1:n); end
    fprintf('\nBLOCK %s  mc_seed_base=%d  trials=%d  G0b(geometry SNR-invariant, ordered)=%d\n', ...
        blk_ts{b}, base, n, g0b);
    res = nan(n, S, NC); recon = false(n,1); th_i = nan(n,1); r_i = nan(n,1);
    tic;
    parfor i = 1:n
        [th, r, Wi] = loc_recon(seeds(i), P);
        recon(i) = strcmp(sprintf('%.10g', th*180/pi), th_str(i)) && strcmp(sprintf('%.10g', r), r_str(i));
        th_i(i) = th; r_i(i) = r;
        row = nan(S, NC);
        for s = 1:S
            N0 = 10^(-SNR_VEC(s)/10);
            row(s, :) = loc_trial_row(th, r, th0, r0, N0, P, Wl, Wi, usw, c1, c0);
        end
        res(i, :, :) = row;
    end
    fprintf('BLOCK %s wall time %.1f s\n', blk_ts{b}, toc);
    recon_total = recon_total + sum(recon); n_total = n_total + n;
    fprintf('G0a (geometry reconstructed from mc_seed, %%.10g string-exact): %d / %d\n', sum(recon), n);
    loc_write_trials(fullfile(OUT_DIR, sprintf('t50_a4_trials_%s_from_%s.csv', ts, blk_ts{b})), ...
        base, seeds, th_i, r_i, recon, SNR_VEC, res, colnames);

    % ---- nominal-geometry references (same conventions) ----
    ref = nan(S, 6);   % [c50_eq_diag c50_eq_fj c50_eq_th fa_eq_diag fa_eq_fj fa_eq_th] variances
    ref_pub = nan(S, 2); % [c50 pinv meansqrtW diag, c50 eq meansqrtW diag] in m
    for s = 1:S
        N0 = 10^(-SNR_VEC(s)/10);
        v50 = loc_vars(th0, r0, N0, P, Wl, c1, c0);
        vfa = loc_vars(th0, r0, N0, P, {eye(P.M)}, c1, c0);
        ref(s, :) = [v50(5:7), vfa(5:7)];
        ref_pub(s, :) = [v50(4), v50(8)];
    end

    % ---- closed-form full-array geometry factors (spec Sec. 5, F-133) ----
    Fg  = r_i.^2 ./ sin(th_i).^2;  F0 = r0^2 / sin(th0)^2;
    cf_r  = sqrt(mean(Fg.^2)) / F0;
    cf_th = sqrt(mean(1 ./ sin(th_i).^2)) * sin(th0);
    FAfac(b) = cf_r;
    fprintf('CLOSED FORM full-array: r factor %.6g (%.4f dB), theta factor %.6g\n', ...
        cf_r, 20*log10(cf_r), cf_th);

    c = @(nm) find(strcmp(colnames, nm));
    for s = 1:S
        snr = SNR_VEC(s);
        V = squeeze(res(:, s, :));
        max_cond = max(max_cond, max(V(:, [c('c50_cond_eq_max'), c('cm_cond_eq'), c('fa_cond_eq')]), [], 'all'));
        % T-50, mean-sqrt, nominal: primary conventions (eq inverse)
        fam = {'c50','cm','fa'};
        for fi = 1:3
          for fm = {'r_diag','r_fj','theta'}
            for iv = {'eq','pinv'}
              nm = sprintf('%s_%s_v%s', fam{fi}, iv{1}, fm{1});
              if isempty(c(nm)), continue; end
              v = V(:, c(nm));
              T50 = sqrt(mean(v)); MSG = mean(sqrt(max(v,0)));
              if strcmp(iv{1}, 'eq')
                  switch fam{fi}
                    case 'c50', k = find(strcmp({'r_diag','r_fj','theta'}, fm{1}));
                                g0 = sqrt(ref(s, k));
                    case 'fa',  k = find(strcmp({'r_diag','r_fj','theta'}, fm{1}));
                                g0 = sqrt(ref(s, 3 + k));
                    case 'cm',  g0 = sqrt(mean(V(:, c(sprintf('cmg0_eq_v%s', fm{1})))));
                  end
              else
                  g0 = NaN;
              end
              unit = 'm'; if strcmp(fm{1}, 'theta'), unit = 'deg'; end
              summ(end+1, :) = {'snr', blk_ts{b}, base, snr, fam{fi}, iv{1}, fm{1}, 'T50',        T50, unit}; %#ok<SAGROW>
              summ(end+1, :) = {'snr', blk_ts{b}, base, snr, fam{fi}, iv{1}, fm{1}, 'meansqrt_g', MSG, unit}; %#ok<SAGROW>
              summ(end+1, :) = {'snr', blk_ts{b}, base, snr, fam{fi}, iv{1}, fm{1}, 'g0',         g0,  unit}; %#ok<SAGROW>
            end
          end
        end
        summ(end+1, :) = {'snr', blk_ts{b}, base, snr, 'c50', 'pinv', 'r_diag', 'g0_meansqrtW_published', ref_pub(s,1), 'm'}; %#ok<SAGROW>
        summ(end+1, :) = {'snr', blk_ts{b}, base, snr, 'c50', 'eq',   'r_diag', 'g0_meansqrtW',           ref_pub(s,2), 'm'}; %#ok<SAGROW>

        % ---- gates computed per cell ----
        fa_T50 = sqrt(mean(V(:, c('fa_eq_vr_diag')))); fa_g0 = sqrt(ref(s, 4));
        fa_th  = sqrt(mean(V(:, c('fa_eq_vtheta'))));  fa_g0t = sqrt(ref(s, 6));
        rel_r  = abs((fa_T50/fa_g0)/cf_r - 1);  rel_t = abs((fa_th/fa_g0t)/cf_th - 1);
        G4 = rel_r <= 1e-6 && rel_t <= 1e-6;
        c50_fac = sqrt(mean(V(:, c('c50_eq_vr_diag')))) / sqrt(ref(s, 1));
        c50_fj  = sqrt(mean(V(:, c('c50_eq_vr_fj'))))   / sqrt(ref(s, 2));
        cm_fac  = sqrt(mean(V(:, c('cm_eq_vr_diag'))))  / sqrt(mean(V(:, c('cmg0_eq_vr_diag'))));
        G5 = c50_fac > 1 && c50_fj > 1 && cm_fac > 1;
        D  = 20*log10(c50_fac) - 20*log10(fa_T50/fa_g0);
        Dspread(b, s) = D;
        X  = V(:, c('cm_eq_vr_fj')) - V(:, c('c50_eq_vr_fj'));
        z  = mean(X) / (std(X)/sqrt(numel(X)));
        tr_dB = 20*log10(sqrt(mean(V(:, c('c50_pinv_vr_diag')))) / sqrt(mean(V(:, c('c50_eq_vr_diag')))));
        all_ok_G4 = all_ok_G4 && G4; all_ok_G5 = all_ok_G5 && G5;
                uf_dB = 20*log10(sqrt(mean(V(:, c('ufa_eq_vr'))))   / sqrt(mean(V(:, c('fa_eq_vr_fj')))));
        uc_dB = 20*log10(sqrt(mean(V(:, c('ucm_eq_vr'))))   / sqrt(mean(V(:, c('cm_eq_vr_fj')))));
        summ(end+1, :) = {'snr', blk_ts{b}, base, snr, 'ufa', 'eq', 'r', 'T50', sqrt(mean(V(:, c('ufa_eq_vr')))), 'm'};
        summ(end+1, :) = {'snr', blk_ts{b}, base, snr, 'ucm', 'eq', 'r', 'T50', sqrt(mean(V(:, c('ucm_eq_vr')))), 'm'};
        fprintf('           legU (USW minus Fresnel, FJ): full array %+.4f dB | matched W %+.4f dB\n', uf_dB, uc_dB);
        gate_rows(end+1, :) = {blk_ts{b}, base, snr, rel_r, rel_t, double(G4), ...
            20*log10(c50_fac), 20*log10(c50_fj), 20*log10(cm_fac), double(G5), D, z, tr_dB}; %#ok<SAGROW>
        fprintf(['  SNR %+5.1f | G4 relerr r %.2e th %.2e %s | c50 fac %.4f dB (fj %.4f) cm fac %.4f dB %s' ...
                 ' | D %.4f dB | z(cm-c50) %+.2f | trunc(pinv-eq) %.4f dB\n'], ...
            snr, rel_r, rel_t, loc_pf(G4), 20*log10(c50_fac), 20*log10(c50_fj), 20*log10(cm_fac), loc_pf(G5), D, z, tr_dB);
    end
end
fprintf('\nG0a TOTAL: %d / %d geometries reconstructed\n', recon_total, n_total);
fprintf('G4 (full-array closed-form identity, all 65 cells): %s\n', loc_pf(all_ok_G4));
fprintf('G5 (compressed geometry-matched > nominal, all 65 cells, c50 diag, c50 fj, cm diag): %s\n', loc_pf(all_ok_G5));
fa_spread = 20*log10(max(FAfac)) - 20*log10(min(FAfac));
D_spread  = max(Dspread, [], 1) - min(Dspread, [], 1);
fprintf('G7 full-array closed-form block spread = %.4f dB; max over SNR of D spread = %.4f dB -> %s\n', ...
    fa_spread, max(D_spread), loc_pf(max(D_spread) <= fa_spread));
fprintf('G9 max equilibrated cond over all compressed and full-array evaluations = %.4g -> %s\n', ...
    max_cond, loc_pf(max_cond <= 4.5e9));

%% ---- 6. Pinned leg (T-52-B): bound at the pinned scene -----------------
Tp = loc_read_pinned_seeds(fullfile(A3_DIR, ['mc_convergence_trials_' T52B_TS '.csv']));
[th0c, r0c] = loc_nominal(Pc);
thp = Pc.conv_theta; rp = Pc.conv_r;
Wlc = loc_W_fixed(Pc);
for snr = [0 5 10 15]
    N0 = 10^(-snr/10);
    seeds = Tp.seeds(Tp.snr == snr);
    vcm = nan(numel(seeds), 3); vcm0 = nan(numel(seeds), 3);
    parfor i = 1:numel(seeds)
        [~, ~, Wi] = loc_recon(seeds(i), Pc);   % W from the same stream positions (driver L459, L609)
        a1 = loc_vars(thp,  rp,  N0, Pc, {Wi}, c1, c0);
        a2 = loc_vars(th0c, r0c, N0, Pc, {Wi}, c1, c0);
        vcm(i, :) = a1(5:7); vcm0(i, :) = a2(5:7);
    end
    v50p = loc_vars(thp, rp, N0, Pc, Wlc, c1, c0);  v500 = loc_vars(th0c, r0c, N0, Pc, Wlc, c1, c0);
    vfap = loc_vars(thp, rp, N0, Pc, {eye(Pc.M)}, c1, c0); vfa0 = loc_vars(th0c, r0c, N0, Pc, {eye(Pc.M)}, c1, c0);
    nm = {'r_diag','r_fj','theta'};
    for k = 1:3
        unit = 'm'; if k == 3, unit = 'deg'; end
        summ(end+1,:) = {'pinned', T52B_TS, 1000, snr, 'c50', 'eq', nm{k}, 'g_pin', sqrt(v50p(4+k)), unit}; %#ok<SAGROW>
        summ(end+1,:) = {'pinned', T52B_TS, 1000, snr, 'c50', 'eq', nm{k}, 'g0',    sqrt(v500(4+k)), unit}; %#ok<SAGROW>
        summ(end+1,:) = {'pinned', T52B_TS, 1000, snr, 'cm',  'eq', nm{k}, 'g_pin', sqrt(mean(vcm(:,k))),  unit}; %#ok<SAGROW>
        summ(end+1,:) = {'pinned', T52B_TS, 1000, snr, 'cm',  'eq', nm{k}, 'g0',    sqrt(mean(vcm0(:,k))), unit}; %#ok<SAGROW>
        summ(end+1,:) = {'pinned', T52B_TS, 1000, snr, 'fa',  'eq', nm{k}, 'g_pin', sqrt(vfap(4+k)), unit}; %#ok<SAGROW>
        summ(end+1,:) = {'pinned', T52B_TS, 1000, snr, 'fa',  'eq', nm{k}, 'g0',    sqrt(vfa0(4+k)), unit}; %#ok<SAGROW>
    end
    fprintf(['PINNED SNR %+3d | fa factor %.6g (closed form 0.330598) | c50 factor %.6g | ' ...
             'cm factor %.6g | c50 pinv meansqrtW at g_pin %.6g m\n'], snr, ...
        sqrt(vfap(5))/sqrt(vfa0(5)), sqrt(v50p(5))/sqrt(v500(5)), ...
        sqrt(mean(vcm(:,1)))/sqrt(mean(vcm0(:,1))), v50p(4));
end

%% ---- 7. Write summary and gate CSVs ------------------------------------
Ts = cell2table(summ, 'VariableNames', {'leg','source_ts','mc_seed_base','SNR_dB', ...
    'family','inverse','form','agg','value','unit'});
writetable(Ts, fullfile(OUT_DIR, sprintf('t50_a4_summary_%s.csv', ts)));
Tg = cell2table(gate_rows, 'VariableNames', {'source_ts','mc_seed_base','SNR_dB', ...
    'G4_relerr_r','G4_relerr_theta','G4_pass','c50_fac_diag_dB','c50_fac_fj_dB','cm_fac_diag_dB', ...
    'G5_pass','D_dB','z_cm_minus_c50_fj','trunc_pinv_minus_eq_dB'});
writetable(Tg, fullfile(OUT_DIR, sprintf('t50_a4_gates_%s.csv', ts)));
fprintf('A4 COMPLETE %s\n', datestr(now));
diary off;


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================
function [th0, r0] = loc_nominal(Q)
% run_monte_carlo_paperC.m L217-219 (live)
th0 = (Q.theta_lo + Q.theta_hi) / 2;
r0  = Q.r_lo_fac * Q.r_RD + (Q.r_hi_fac - Q.r_lo_fac) * Q.r_RD / 2;
end

function Wl = loc_W_fixed(Q)
% wb_crb_compressed.m L173-174 (live), seeds rng_seed_W .. rng_seed_W+N_seed-1
Wl = cell(1, Q.N_seed); s0 = rng;
for s = 1:Q.N_seed
    rng(Q.rng_seed_W + s - 1, 'twister');
    Wl{s} = (1/sqrt(Q.M)) * exp(1j * 2*pi * rand(Q.M, Q.N_RF));
end
rng(s0);
end

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

function v = loc_vars(theta, r, N0, Q, Wl, c1, c0)
% Returns [pinv: vr_diag vr_fj vtheta sr_diag_meansqrtW, eq: vr_diag vr_fj vtheta sr_diag_meansqrtW,
%          n_sing_mean, cond_eq_max]; variances averaged over the combiners in Wl.
o = struct(); o.W_list = Wl; o.return_raw = true;
[~, ~, info] = wb_crb_compressed(theta, r, 1, N0, Q, o);
kap  = c0 * sin(theta)^2 / r;
g    = [-2*r*cos(theta) / (c1*sin(theta)^2); -r/kap];   % d r / d(omega, kappa)
dkdr = c0 * sin(theta)^2 / r^2;                          % |d kappa / d r|, L155
dwdt = c1 * sin(theta);                                  % |d omega / d theta|
R = {info.raw.pinv, info.raw.eq};
v = nan(1, 10);
for e = 1:2
    Coo = R{e}.C_oo; Ckk = R{e}.C_kk; Cok = R{e}.C_ok;
    vr_d = max(Ckk, 0) / dkdr^2;                         % L203-206 convention
    vr_f = g(1)^2*Coo + 2*g(1)*g(2)*Cok + g(2)^2*Ckk;    % full-Jacobian
    vt   = max(Coo, 0) / dwdt^2 * (180/pi)^2;            % L202, L209 convention
    v(4*(e-1) + (1:4)) = [mean(vr_d), mean(vr_f), mean(vt), mean(sqrt(vr_d))];
end
v(9)  = mean(info.n_singular_all);
v(10) = max(info.raw.cond_eq_all);
end

function row = loc_trial_row(th, r, th0, r0, N0, Q, Wl, Wi, usw, c1, c0)
a50 = loc_vars(th,  r,  N0, Q, Wl,          c1, c0);
am  = loc_vars(th,  r,  N0, Q, {Wi},        c1, c0);
am0 = loc_vars(th0, r0, N0, Q, {Wi},        c1, c0);
af  = loc_vars(th,  r,  N0, Q, {eye(Q.M)},  c1, c0);
u = nan(1, 4);
if usw.valid
    [u(1), u(2)] = loc_crb_usw(th, r, N0, Q, eye(Q.M), usw.delta, usw.sgn);
    [u(3), u(4)] = loc_crb_usw(th, r, N0, Q, Wi,       usw.delta, usw.sgn);
end
row = [a50([1 2 3 4 5 6 7 8 9 10]), ...
       am([1 2 3 5 6 7 9 10]), ...
       am0([5 6 7]), ...
       af([1 5 6 7 9 10]), u];
end

function nm = loc_colnames()
nm = {'c50_pinv_vr_diag','c50_pinv_vr_fj','c50_pinv_vtheta','c50_pinv_sr_diag_meansqrtW', ...
      'c50_eq_vr_diag','c50_eq_vr_fj','c50_eq_vtheta','c50_eq_sr_diag_meansqrtW', ...
      'c50_nsing_mean','c50_cond_eq_max', ...
      'cm_pinv_vr_diag','cm_pinv_vr_fj','cm_pinv_vtheta','cm_eq_vr_diag','cm_eq_vr_fj','cm_eq_vtheta', ...
      'cm_nsing','cm_cond_eq', ...
      'cmg0_eq_vr_diag','cmg0_eq_vr_fj','cmg0_eq_vtheta', ...
      'fa_pinv_vr_diag','fa_eq_vr_diag','fa_eq_vr_fj','fa_eq_vtheta','fa_nsing','fa_cond_eq', ...
      'ufa_eq_vr','ufa_eq_vtheta','ucm_eq_vr','ucm_eq_vtheta'};
end

function [Vr, Vth] = loc_crb_usw(theta, r, N0, Q, W, delta, sgn)
% Slepian-Bangs FIM in (theta, r, p, N0) under the exact-distance phase-only
% steering matched to nf_usw_steer by PC-3; p = 1; same noise model and
% regularisation as wb_crb_compressed.m L290, L307; equilibrated exact inverse.
J = zeros(4);
WtW = W' * W;
for ks = 1:Q.K_s
    kk   = 2*pi*Q.alpha_k_vec(ks) / Q.lambda_c;
    dist = sqrt(r^2 + delta.^2 - 2*r*delta*cos(theta));
    a    = exp(sgn*1j*kk*(dist - r));
    da_t = sgn*1j*kk*(r*delta*sin(theta) ./ dist) .* a;
    da_r = sgn*1j*kk*((r - delta*cos(theta)) ./ dist - 1) .* a;
    dv = W' * a;  dt = W' * da_t;  dr = W' * da_r;
    Ry = N0*WtW + dv*dv';
    Ry = (Ry + Ry')/2 + 1e-12*eye(size(W, 2));
    Ri = inv(Ry); %#ok<MINV>
    dR = {dt*dv' + dv*dt', dr*dv' + dv*dr', dv*dv', WtW};
    for i1 = 1:4
        A = Ri * dR{i1};
        for i2 = i1:4
            val = Q.N * real(trace(A * Ri * dR{i2}));
            J(i1, i2) = J(i1, i2) + val;
            if i2 ~= i1, J(i2, i1) = J(i2, i1) + val; end
        end
    end
end
dj = sqrt(diag(J)); Dinv = diag(1 ./ dj);
Jn = Dinv * J * Dinv; Jn = (Jn + Jn.')/2;
C  = Dinv * (Jn \ eye(4)) * Dinv;
Vth = C(1,1) * (180/pi)^2;
Vr  = C(2,2);
end

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

function Tp = loc_read_pinned_seeds(csv_path)
L = splitlines(string(fileread(csv_path))); L = L(strlength(L) > 0);
hdr = split(L(1), ','); C = split(L(2:end), ',');
Tp.snr   = str2double(C(:, hdr == "SNR_dB"));
Tp.seeds = str2double(C(:, hdr == "mc_seed"));
end

function loc_write_trials(fn, base, seeds, th, r, recon, SNR_VEC, res, colnames)
fid = fopen(fn, 'w');
fprintf(fid, 'mc_seed_base,mc_idx,mc_seed,theta_true_deg,r_true_m,recon_ok,SNR_dB,%s\n', strjoin(colnames, ','));
fmt = ['%d,%d,%d,%.10g,%.10g,%d,%.6g', repmat(',%.10g', 1, numel(colnames)), '\n'];
for s = 1:numel(SNR_VEC)
    for i = 1:numel(seeds)
        fprintf(fid, fmt, base, i, seeds(i), th(i)*180/pi, r(i), recon(i), SNR_VEC(s), squeeze(res(i, s, :)));
    end
end
fclose(fid);
end

function s = loc_pf(tf)
if tf, s = 'PASS'; else, s = 'FAIL'; end
end
