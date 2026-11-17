function [theta_s, r_s, info_s] = new5pl_music_start(Y_full, W_comb, P)
%NEW5PL_MUSIC_START  The MUSIC-type comparator's estimate for one trial.
%
%   [THETA_S, R_S, INFO_S] = NEW5PL_MUSIC_START(Y_FULL, W_COMB, P)
%
%   Pilot of O-3 (paperC-lab sessions/2026-11-13_NEW5PL, SPEC_NEW5PL.md).
%   Forms, for ONE trial of run_monte_carlo_paperC, the estimate of the
%   comparator diagnostics/ctA/wb_nf_music_est.m from the compressed
%   snapshots of that trial, as diagnostics/ctA/ctA_sweep_paperC.m forms it
%   (its loc_one_seed, L221-261 at commit 7067ed9), with two differences:
%     (1) one trial per call, where the sweep stacks the 13 SNR values of a
%         seed in one call (its L227, L243, L251);
%     (2) OPT.eta_true is not set: no value of the trial's truth is passed
%         (wb_nf_music_est.m reads it only for its post-hoc capture flag).
%   The estimator and its helpers are called as they are; nothing of them
%   is copied here but the two local functions below, which another file
%   cannot reach.
%
%   DRAFT. This file was written as text. It has not been parsed or run.
%
%   INPUTS
%     Y_full   N_RF x N x K_s compressed snapshots of the trial
%              (wb_channel_gen_ofdm_nf, second output).
%     W_comb   M x N_RF combiner of the trial (its eighth output).
%     P        production struct of the run.
%   OUTPUTS
%     theta_s  angle [rad]: acos of omega_hat/q, clamped to [-1, 1], the map
%              of wb_nf_music_est.m L235.
%     r_s      range [m]: info.r_hat of the estimator (its L234, L267).
%     info_s   struct: eta (1 x 2), theta_hat_deg, n_iter, stop_reason,
%              verdict of the estimator, for a caller that wants them.
%
%   The constants K_COARSE = 8, N_CAND = 6 and maxit = 400 are those of the
%   CT-A run of record (RUN_SHEET_CTA.md L263-265 of paperC-lab
%   sessions/2026-10-30_CTA-predata; ctA_sweep_paperC.m L74).
%
%   No random number is drawn here. No bound is computed here.
%   Requires diagnostics/ctA on the path. Base MATLAB only. 7-bit ASCII.

narginchk(3, 3);
assert(P.d == 1, 'new5pl_music_start: the estimator takes one path');
assert(ndims(Y_full) == 3 && size(Y_full, 1) == P.N_RF && size(Y_full, 3) == P.K_s, ...
    'new5pl_music_start: Y_full must be N_RF x N x K_s');

L = ctA1_locals();

% The constants of ctA_sweep_paperC.m L67-69.
q_scale = 2*pi*P.d_ant/P.lambda_c;
s0      = pi*P.d_ant^2/P.lambda_c;
m_bar   = ((0:P.M-1) - (P.M-1)/2).';

cfg  = struct('K_COARSE', 8, 'N_CAND', 6);
GRID = loc_grid(P, s0, m_bar, cfg);

% The options of ctA_sweep_paperC.m L76-79, with batch 1 (one trial).
OPT = struct('s0', s0, 'q', q_scale, 'npar', 0, ...
             'batch', 1, 'assert_internal', true, ...
             'tangent', true, ...
             'tol_dec', 2e-12, 'tol_ml', 1e-9, 'maxit', 400);

[~, Tw] = L.whiten(W_comb);
U       = loc_eigvecs(Y_full, Tw, P.N);

[eta, info] = wb_nf_music_est(U, W_comb, Tw, P.alpha_k_vec, m_bar, GRID, OPT);

theta_s = acos(max(-1, min(1, eta(1)/q_scale)));
r_s     = info.r_hat(1);

info_s = struct();
info_s.eta           = eta;
info_s.theta_hat_deg = info.theta_hat_deg(1);
info_s.n_iter        = info.n_iter(1);
info_s.stop_reason   = info.stop_reason{1};
info_s.verdict       = info.verdict{1};

end

%% ========================================================================
%  LOCAL FUNCTIONS -- copied by script, byte for byte, from
%  diagnostics/ctA/ctA_sweep_paperC.m at commit 7067ed9 (L263-278, L280-313).
%  DO NOT EDIT: patch_new5pl.pl (mode copy) of the session compares them.
% =========================================================================

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
