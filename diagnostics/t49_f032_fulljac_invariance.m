% =========================================================================
% t49_f032_fulljac_invariance.m
% =========================================================================
% T-49 -- closure run for F-032 (Addendum B.3).
%
% QUESTION. F-028's full-array wideband range-gain invariance was measured in
% CURVATURE coordinates (the Eq. (19) convention). Does the same statement
% hold when the range bound is formed through the FULL Jacobian
%
%     V_r = g_omega^2 C_ww + 2 g_omega g_kappa C_wk + g_kappa^2 C_kk ,
%     g_omega = 2 r cot(theta) / (d omega / d theta),   g_kappa = -r / kappa,
%
% rather than through the published marginal propagation g_kappa^2 C_kk?
%
% STATUS ON ENTRY (see spec Sec. 2.1). PARTIALLY ANSWERED ALREADY. The D1b
% block of stageD_ctg_o2_o3_s0c.m ran exactly this factorial at W = I and
% B = 400 MHz; stageD_console.txt lines 50-55 record
%     G10_FJ = +2.174020e-07 dB, G01_FJ = +7.387041e-05 dB,
%     G11_FJ = +7.408781e-05 dB, |T_cross|/T_kappa = 2.660e-17.
% This run CONFIRMS that measurement and extends it over bandwidth, range,
% angle, and the Paper C array configuration, which the D1b block did not do.
%
% METHOD. Factorial over two switches applied to the subcarrier scaling:
%     bw_k = 1 + u*(alpha_k - 1)   (beam squint in omega)
%     bk_k = 1 + v*(alpha_k - 1)   (frequency-dependent curvature in kappa)
% (u,v) = (0,0) is the narrowband replicate: K_s identical manifolds, so its
% bound already contains the full 10log10(K_s) data-diversity factor. Every
% reported gain is therefore the RESIDUAL after data diversity, by
% construction, and is directly comparable to the closed form.
%
% CLOSED FORM. residual_dB = 10 log10( 1 + mean_k( (alpha_k - 1)^2 ) ).
% For an exactly uniform band this reduces to 10 log10(1 + (1/12)(B/f_c)^2).
% The grid's own second moment is used, not the idealised one: the two differ
% in the fourth significant figure whenever the grid is a subsample.
%
% PASS CRITERION (pre-registered, three-way, per configuration point):
%     SHARP-PASS  |G11_FJ - residual_dB| <= 1e-8 dB
%     PASS        |G11_FJ - residual_dB| <= 1e-6 dB
%     MARGINAL    1e-6 dB < |.| <= 3e-4 dB
%     FAIL        |.| > 3e-4 dB   (breaks the published claim tolerance)
% Additional required conditions at every point:
%     C1  |T_cross| / T_kappa <= 1e-12          (cross term vanishes at W = I)
%     C2  |G11_FJ - (G10_FJ + G01_FJ)| <= 1e-9 dB   (factorial additivity)
%     C3  loc_Cinv_raw's internal C*J = I assertion passes (built in)
%
% LOCAL HELPERS loc_subcarrier_grid, loc_fim_wb, loc_fim_single,
% loc_Cinv_raw and tern below are VERBATIM copies of the functions of the
% same names in stageD_ctg_o2_o3_s0c.m (lines 854-868, 887-925, 992-1006,
% 1313-1315). Do not edit them. loc_Cinv_raw is the CORRECTED inverse form;
% it is NOT the defective loc_Cinv_entry of stageA_s0_audit.m (Addendum B.6).
%
% Author-side execution. MATLAB R2025b, base MATLAB only. Provenance P-RUN.
% Encoding: 7-bit ASCII only.
% =========================================================================

clear; close all; clc;
stamp = datestr(now, 'yyyymmdd_HHMMSS');
diary(sprintf('t49_f032_console_%s.txt', stamp)); diary on;
t_start = tic;

fprintf('=============================================================\n');
fprintf('  T-49  F-032 CLOSURE -- full-Jacobian invariance at W = I\n');
fprintf('=============================================================\n\n');

%% ---- RUN CONTROL -------------------------------------------------------
NPAR = 4;                       % set to numcores-2 to use parfor in loc_fim_wb
TOL_SHARP = 1e-8;               % dB
TOL_PASS  = 1e-6;               % dB
TOL_FAIL  = 3e-4;               % dB, the published claim tolerance
TOL_CROSS = 1e-12;              % |T_cross| / T_kappa at W = I
TOL_ADD   = 1e-9;               % dB, factorial additivity

c0   = 3e8;                     % NOT 299792458 -- matches committed CSVs
fc   = 28e9;
lam  = c0 / fc;
dant = lam / 2;

% Two configurations. LEG 1 reproduces and extends the GLOBECOM companion's
% setting (the object F-032 is actually about). LEG 2 asks the same question
% at Paper C's own array and subcarrier grid.
LEG(1).name    = 'GLOBECOM';  LEG(1).M = 256; LEG(1).Delta_f = 120e3;
LEG(2).name    = 'PaperC';    LEG(2).M =  64; LEG(2).Delta_f =  25e6;

B_LIST     = [100 200 400 600 800] * 1e6;
R_LIST     = [1.50 2.13 3.00 5.00 10.00];       % m; 5 m is the D1b anchor
THETA_LIST = [20 40 60];                         % deg
SNR_dB_ref = 10;
N_snap     = 64;
p_true     = 1;
N0_ref     = p_true / (10^(SNR_dB_ref/10));

ROWS = {};
GATES = {};

for iL = 1:numel(LEG)
    M       = LEG(iL).M;
    Delta_f = LEG(iL).Delta_f;
    m_bar   = ((0:M-1) - (M-1)/2).';
    W_eye   = eye(M);
    q_scale = 2*pi*dant / lam;
    s0      = pi*dant^2 / lam;

    fprintf('### LEG %d -- %s : M = %d, Delta_f = %.3f MHz ###\n\n', ...
            iL, LEG(iL).name, M, Delta_f/1e6);

    for iB = 1:numel(B_LIST)
        B_val = B_LIST(iB);
        [al_sub, Ks] = loc_subcarrier_grid(B_val, fc, Delta_f);
        dl   = al_sub(:) - 1;
        md2  = mean(dl.^2);
        resid_dB = 10*log10(1 + md2);
        ideal_dB = 10*log10(1 + (1/12)*(B_val/fc)^2);

        fprintf('  B = %4.0f MHz | K_s = %4d | mean(delta^2) = %.9e\n', ...
                B_val/1e6, Ks, md2);
        fprintf('    closed form (grid)  = %+.9e dB\n', resid_dB);
        fprintf('    closed form (ideal) = %+.9e dB\n', ideal_dB);

        for iT = 1:numel(THETA_LIST)
            th_deg = THETA_LIST(iT);
            th     = th_deg * pi/180;
            om     = q_scale * cos(th);
            c_ch   = s0 * sin(th)^2;
            dom_dth = -q_scale * sin(th);

            for iR = 1:numel(R_LIST)
                rr = R_LIST(iR);
                ka = c_ch / rr;
                g_om = 2*rr*cot(th) / dom_dth;
                g_ka = -rr / ka;

                uv = [0 0; 1 0; 0 1; 1 1];
                Vc = zeros(4,1); Tc = zeros(4,1); Tk = zeros(4,1);
                To = zeros(4,1); cJ = zeros(4,1);
                for iu = 1:4
                    J = loc_fim_wb(om, ka, p_true, N0_ref, al_sub, ...
                                   uv(iu,1), uv(iu,2), m_bar, W_eye, ...
                                   N_snap, NPAR);
                    Dp    = diag([q_scale, ka, p_true, N0_ref]);
                    Craw  = loc_Cinv_raw(J, Dp);
                    To(iu) = g_om^2 * Craw(1,1);
                    Tk(iu) = g_ka^2 * Craw(2,2);
                    Tc(iu) = 2*g_om*g_ka * Craw(1,2);
                    Vc(iu) = To(iu) + Tk(iu) + Tc(iu);
                    cJ(iu) = cond(J);
                end

                G10 = 10*log10(Vc(1)/Vc(2));
                G01 = 10*log10(Vc(1)/Vc(3));
                G11 = 10*log10(Vc(1)/Vc(4));
                dev     = abs(G11 - resid_dB);
                cross_r = abs(Tc(4)) / max(Tk(4), realmin);
                add_err = abs(G11 - (G10 + G01));

                if     dev <= TOL_SHARP; vd = 'SHARP-PASS';
                elseif dev <= TOL_PASS;  vd = 'PASS';
                elseif dev <= TOL_FAIL;  vd = 'MARGINAL';
                else;                    vd = 'FAIL';
                end
                c1 = cross_r <= TOL_CROSS;
                c2 = add_err <= TOL_ADD;

                ROWS(end+1,:) = {LEG(iL).name, M, B_val, Ks, md2, ...
                    th_deg, rr, g_om, g_ka, To(4), Tk(4), Tc(4), Vc(4), ...
                    sqrt(Vc(4)), G10, G01, G11, resid_dB, ideal_dB, dev, ...
                    cross_r, add_err, max(cJ), vd, ...
                    double(c1), double(c2)}; %#ok<SAGROW>

                if rr == 5 || rr == 2.13
                    fprintf('    th=%2d r=%5.2f | G11=%+.9e | dev=%.3e | %s\n', ...
                            th_deg, rr, G11, dev, vd);
                end
            end
        end
        fprintf('\n');
    end
end

T = cell2table(ROWS, 'VariableNames', {'leg','M','B_hz','K_s','mean_delta2', ...
    'theta_deg','r_m','g_omega','g_kappa','T_omega','T_kappa','T_cross', ...
    'V_fullJac','sqrtV_m','G10_FJ_dB','G01_FJ_dB','G11_FJ_dB', ...
    'closed_grid_dB','closed_ideal_dB','dev_dB','cross_ratio', ...
    'additivity_err_dB','cond_J_max','verdict','C1_pass','C2_pass'});
writetable(T, sprintf('t49_f032_fulljac_%s.csv', stamp));

%% ---- GATE ADJUDICATION -------------------------------------------------
fprintf('=============================================================\n');
fprintf('  GATE SUMMARY\n');
fprintf('=============================================================\n');
n_fail = sum(strcmp(T.verdict, 'FAIL'));
n_marg = sum(strcmp(T.verdict, 'MARGINAL'));
n_pass = sum(strcmp(T.verdict, 'PASS') | strcmp(T.verdict, 'SHARP-PASS'));
fprintf('  points: %d | SHARP/PASS %d | MARGINAL %d | FAIL %d\n', ...
        height(T), n_pass, n_marg, n_fail);
fprintf('  max dev  = %.6e dB at B = %.0f MHz, r = %g m, theta = %g deg\n', ...
        max(T.dev_dB), T.B_hz(argmaxv(T.dev_dB))/1e6, ...
        T.r_m(argmaxv(T.dev_dB)), T.theta_deg(argmaxv(T.dev_dB)));
fprintf('  C1 (cross term) all pass : %s\n', tern(all(T.C1_pass==1),'YES','NO'));
fprintf('  C2 (additivity) all pass : %s\n', tern(all(T.C2_pass==1),'YES','NO'));

% Reproduction check against stageD_console.txt lines 50-55.
ix = strcmp(T.leg,'GLOBECOM') & T.B_hz==400e6 & T.r_m==5 & T.theta_deg==40;
if any(ix)
    fprintf('\n  REPRODUCTION of stageD D1b (GLOBECOM, 400 MHz, r=5 m, 40 deg):\n');
    fprintf('    G10_FJ now = %+.6e   stageD = +2.174020e-07\n', T.G10_FJ_dB(ix));
    fprintf('    G01_FJ now = %+.6e   stageD = +7.387041e-05\n', T.G01_FJ_dB(ix));
    fprintf('    G11_FJ now = %+.6e   stageD = +7.408781e-05\n', T.G11_FJ_dB(ix));
    rep = abs(T.G11_FJ_dB(ix) - 7.408781e-05) < 1e-10;
    fprintf('    GATE REPRO : %s\n', tern(rep, 'PASS', 'FAIL -- STOP'));
    GATES(end+1,:) = {'REPRO_stageD_D1b', tern(rep,'PASS','FAIL')}; %#ok<SAGROW>
end

if n_fail > 0
    verdict = 'FAIL';
elseif n_marg > 0
    verdict = 'MARGINAL';
else
    verdict = 'PASS';
end
GATES(end+1,:) = {'F032_overall', verdict};
fprintf('\n  F-032 OVERALL : %s\n', verdict);
fprintf('  Interpretation is in the T-45 spec Sec. 2.5; do not improvise it.\n');
writetable(cell2table(GATES, 'VariableNames', {'gate','verdict'}), ...
           sprintf('t49_f032_gates_%s.csv', stamp));

fprintf('\n  elapsed %.1f s\n', toc(t_start));
diary off;


%% ========================================================================
%  LOCAL HELPERS -- VERBATIM from stageD_ctg_o2_o3_s0c.m. Do not edit.
% =========================================================================
function i = argmaxv(v)
[~, i] = max(v);
end

function [al_sub, Ks] = loc_subcarrier_grid(B_val, fc, Delta_f)
% VERBATIM from stageA_s0_audit.m / wb_crb_globecom2026.m.
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

function J = loc_fim_wb(omega, kappa, p, N0, al_sub, u, v, m_bar, W, N_snap, npar)
% VERBATIM pattern from stageA_s0_audit.m.
    al_sub = al_sub(:);
    dl = al_sub - 1;
    bw = 1 + u*dl;
    bk = 1 + v*dl;
    Ks = numel(al_sub);
    J  = zeros(4);
    parfor (k = 1:Ks, npar)
        J = J + loc_fim_single(omega, kappa, p, N0, ...
                               bw(k), bk(k), m_bar, W, N_snap);
    end
    J = (J + J.')/2;
end

function J = loc_fim_single(omega, kappa, p, N0, bw, bk, m_bar, W, N_snap)
% VERBATIM from stageA_s0_audit.m.
    WtW = W' * W;
    a   = exp(1j*bw*omega*m_bar - 1j*bk*kappa*m_bar.^2);
    d   = W' * a;
    Ry  = N0 * WtW + p * (d * d');
    Ry_inv = inv(Ry); %#ok<MINV>
    dot_w = W' * ( 1j * bw * (m_bar    .* a) );
    dot_k = W' * (-1j * bk * (m_bar.^2 .* a) );
    dR = cell(4,1);
    dR{1} = p * (dot_w * d' + d * dot_w');
    dR{2} = p * (dot_k * d' + d * dot_k');
    dR{3} = d * d';
    dR{4} = WtW;
    J = zeros(4);
    for ii = 1:4
        Ri_dRi = Ry_inv * dR{ii};
        for jj = ii:4
            val = N_snap * real(trace(Ri_dRi * Ry_inv * dR{jj}));
            J(ii,jj) = val;  J(jj,ii) = val;
        end
    end
end

function Craw = loc_Cinv_raw(J, Dp)
% Full inverse in raw coordinates via the scaled solve.
%   Ju = Dp' J Dp  ==>  Ju^{-1} = Dp^{-1} J^{-1} Dp'^{-1}
%   ==>  J^{-1} = Dp * Ju^{-1} * Dp'
% NOTE: stageA_s0_audit.m loc_Cinv_entry used "Dp \ Ci / Dp" here, which
% is off by Dp^{-2} on each side. It is a latent defect in that script,
% affecting only the diagnostic C_wk column of stageA_A1_nrf_sweep.csv.
% Correlations and signs are invariant under the diagonal rescaling, so
% no Stage A conclusion is touched. Do not copy the old form back.
    Ju   = Dp' * ((J + J')/2) * Dp;
    Craw = Dp * (Ju \ eye(4)) * Dp';
    Craw = (Craw + Craw')/2;
    assert(norm(Craw*((J+J')/2) - eye(size(J,1)), 'fro') < 1e-6, ...
           'loc_Cinv_raw: C*J is not the identity');
end

function s = tern(c, a, b)
    if c; s = a; else; s = b; end
end
