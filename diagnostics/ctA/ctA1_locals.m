% PORT-PROVENANCE BEGIN
% Copied without change from repository rvsenyuva/wb-nf-crb-globecom26
% path matlab/ctA1/ctA1_locals.m
% commit ba1862cb39aa638c0dbc700f8f310e3e7dd12dcf
% blob 250084ae39548cb6fc83038d500d0be2e99f7d01
% CTAQ-D1 (a): port for the CT-A run at Paper C's geometry.
% No line below this block differs from the blob.
% PORT-PROVENANCE END
function L = ctA1_locals()
%CTA1_LOCALS  Verbatim Stage D local functions, exposed as function handles.
%
%   L = CTA1_LOCALS() returns a struct whose fields are handles to the
%   local functions that blocks E3 and E4 of task T-02 (CT-A1) need but
%   that live inside stageD_ctg_o2_o3_s0c.m, where MATLAB scoping makes
%   them unreachable from another file.
%
%   Usage:
%       L = ctA1_locals();
%       J = L.fim_wb(omega, kappa, p, N0, al_ES, 1, 1, m_bar, W, N_snap, 0);
%       C = L.Cinv_raw(J, diag([q, kappa, p, N0]));
%
% =========================================================================
%  DO NOT EDIT THE FUNCTION BODIES IN THIS FILE.
% =========================================================================
%   Every body below is a VERBATIM copy from stageD_ctg_o2_o3_s0c.m. This
%   file exists so there is exactly ONE diff target when checking for
%   drift, rather than one per extracted function. If Stage D changes,
%   re-copy; do not patch here and do not "improve" anything.
%
%   Provenance map (stageD_ctg_o2_o3_s0c.m, 1316 lines):
%
%     field              local function              source lines
%     -----------------  --------------------------  ------------
%     subcarrier_grid_Ks loc_subcarrier_grid_Ks       869-879
%     whiten             loc_whiten                  1008-1014
%     bvec               loc_bvec                    1016-1025
%     Qbar               loc_Qbar                    1027-1032
%     hessian            loc_hessian                 1127-1167
%     fim_wb             loc_fim_wb                   887-900
%     fim_single         loc_fim_single               902-924
%     Cinv_raw           loc_Cinv_raw                 992-1006
%     Vr_full            loc_Vr_full                  984-990
%     Vr_exact           loc_Vr_exact                 976-982
%     prctile            loc_prctile                 1303-1311
%
%   Internal call graph, all satisfied within this file:
%     fim_wb  -> fim_single
%     hessian -> bvec, Qbar
%
%   NOTES CARRIED FROM STAGE D, still binding:
%     - loc_Cinv_raw carries the correction to stageA_s0_audit.m's
%       loc_Cinv_entry, which used "Dp \ Ci / Dp" and was off by Dp^{-2} on
%       each side. Do not copy the old form back.
%     - loc_hessian is the POPULATION Hessian at the truth. It is not
%       superseded by wb_nf_music_grad; it is the independent anchor for
%       G-E2 (spec v2 Sec. 6.4).
%     - loc_Vr_exact is Eq. (19) of the paper, NOT the full-Jacobian bound.
%       The E4 comparator is Vr_full (spec v2 Sec. 3 and Sec. 11).
%
%   EXPECTED checkcode OUTPUT -- these three messages are CORRECT and must
%   NOT be silenced. Silencing them would require editing the bodies, which
%   would destroy the byte-identical relationship with Stage D that is this
%   file's entire purpose.
%
%     L 176  "message was once suppressed here, but is no longer generated"
%            -- the %#ok<MINV> on Ry_inv = inv(Ry). R2026's analyser moved
%               the message to the point of use.
%     L 186  INV(A)*b  ) both inside loc_fim_single's trace loop, where
%     L 188  INV(A)*b  ) Ry_inv is reused across ten products.
%
%   On the substance of the INV messages: Ry = N0*W^H W + p d d^H is well
%   conditioned in this configuration, the explicit inverse is reused ten
%   times per subcarrier so a solve-based form is not obviously cheaper,
%   and every committed P-RUN value in the Ledger was produced by exactly
%   this arithmetic. Whether the choice costs any accuracy is not a matter
%   of opinion here: G-E8 validates this FIM against an independent
%   finite-difference Slepian-Bangs construction at relative Frobenius
%   1e-6, which is three orders tighter than anything inv() could cost at
%   this conditioning. If G-E8 passes, the question is closed empirically.
%
%   Base MATLAB only. Encoding: 7-bit ASCII.
%
%   See also WB_NF_MUSIC_EST, WB_NF_MUSIC_GRAD, WB_NF_GEN_EIGVECS.

L = struct( ...
    'subcarrier_grid_Ks', @loc_subcarrier_grid_Ks, ...
    'whiten',             @loc_whiten, ...
    'bvec',               @loc_bvec, ...
    'Qbar',               @loc_Qbar, ...
    'hessian',            @loc_hessian, ...
    'fim_wb',             @loc_fim_wb, ...
    'fim_single',         @loc_fim_single, ...
    'Cinv_raw',           @loc_Cinv_raw, ...
    'Vr_full',            @loc_Vr_full, ...
    'Vr_exact',           @loc_Vr_exact, ...
    'prctile',            @loc_prctile);

end

% =========================================================================
%  VERBATIM BODIES -- stageD_ctg_o2_o3_s0c.m
% =========================================================================

function [al_sub, Ks] = loc_subcarrier_grid_Ks(B_val, fc, Delta_f, Ks_target)
% Estimator grid: same physical band, Ks_target retained subcarriers,
% endpoints included. Register decision D2 fixes Ks_target = 64.
    K_val = max(1, round(B_val / Delta_f));
    Ks    = min(Ks_target, K_val);
    k_v   = (0:K_val-1).';
    f_v   = fc + (k_v - (K_val-1)/2) * Delta_f;
    al_v  = f_v / fc;
    idx_s = round(linspace(1, K_val, Ks));
    al_sub = al_v(idx_s);
end

function [Gw, Tw] = loc_whiten(W)
% G_W = W^H W = L L^H,  T_W = L^{-1} by triangular solve (never inv()).
    Gw = W'*W;
    Gw = (Gw + Gw')/2;
    L  = chol(Gw, 'lower');
    Tw = L \ eye(size(L,1));
end

function b = loc_bvec(W, Tw, al_sub, m_bar, mb2, omega, kappa)
% Whitened compressed steering vectors, N_RF x K_s.
    al_sub = al_sub(:);
    Ks = numel(al_sub);
    b  = zeros(size(W,2), Ks);
    for k = 1:Ks
        a = exp(1j*al_sub(k)*omega*m_bar - 1j*al_sub(k)*kappa*mb2);
        b(:,k) = Tw * (W' * a);
    end
end

function Q = loc_Qbar(b, b0)
% Normalised population MUSIC cost, mean over subcarriers of 1 - |c_k|^2.
    num = sum(conj(b).*b0, 1);
    Q = mean(1 - abs(num).^2 ./ (sum(abs(b).^2,1) .* sum(abs(b0).^2,1)));
    Q = max(real(Q), 0);
end

function [H, Hnum] = loc_hessian(W, Tw, al_sub, m_bar, mb2, om0, ka0)
% Hessian of Qbar at the truth, analytic and by central differences.
% Analytic:  H_ij = 2 * mean_k Re{ bdot_i' P_n bdot_j } / ||b0||^2.
    al_sub = al_sub(:);
    Ks = numel(al_sub);
    H  = zeros(2);
    for k = 1:Ks
        al = al_sub(k);
        a  = exp(1j*al*om0*m_bar - 1j*al*ka0*mb2);
        b  = Tw*(W'*a);
        e  = real(b'*b);
        u  = b/sqrt(e);
        d1 = Tw*(W'*( 1j*al*(m_bar.*a)));
        d2 = Tw*(W'*(-1j*al*(mb2  .*a)));
        D  = [d1, d2];
        Pn = eye(numel(b)) - u*u';
        H  = H + 2*real(D'*Pn*D)/e;
    end
    H = H/Ks;
    H = (H + H')/2;

    b0 = loc_bvec(W, Tw, al_sub, m_bar, mb2, om0, ka0);
    h  = [1e-5*abs(om0), 1e-5*abs(ka0)];
    Hnum = zeros(2);
    for ii = 1:2
        for jj = ii:2
            s = 0;
            for si = [-1 1]
                for sj = [-1 1]
                    dv = [0 0]; dv(ii) = dv(ii)+si*h(ii);
                    dv(jj) = dv(jj)+sj*h(jj);
                    bb = loc_bvec(W, Tw, al_sub, m_bar, mb2, ...
                                  om0+dv(1), ka0+dv(2));
                    s = s + si*sj*loc_Qbar(bb, b0);
                end
            end
            val = s/(4*h(ii)*h(jj));
            Hnum(ii,jj) = val; Hnum(jj,ii) = val;
        end
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

function V = loc_Vr_full(J, g_omega, g_kappa, q, kappa, p, N0)
% VERBATIM from stageA_s0_audit.m.
    Dp = diag([q, kappa, p, N0]);
    Ju = Dp' * ((J + J')/2) * Dp;
    gu = Dp' * [g_omega; g_kappa; 0; 0];
    V  = gu' * (Ju \ gu);
end

function V = loc_Vr_exact(J, r, kappa, c_chirp, q, p, N0)
% VERBATIM from stageA_s0_audit.m. Eq. (19) of the paper.
    Dp = diag([q, kappa, p, N0]);
    Ju = Dp' * ((J + J')/2) * Dp;
    gu = Dp' * [0;1;0;0];
    V  = (r^2/c_chirp)^2 * (gu' * (Ju \ gu));
end

function y = loc_prctile(v, pct)
% VERBATIM from stageA_s0_audit.m.
    v = sort(v(:)); n = numel(v);
    if n == 1; y = v; return; end
    pos = pct/100 * n + 0.5;
    if pos <= 1;  y = v(1);  return; end
    if pos >= n;  y = v(n);  return; end
    lo = floor(pos); y = v(lo) + (pos-lo)*(v(lo+1)-v(lo));
end
