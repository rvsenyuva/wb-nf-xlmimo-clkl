% PORT-PROVENANCE BEGIN
% Copied without change from repository rvsenyuva/wb-nf-crb-globecom26
% path matlab/ctA1/wb_nf_music_grad.m
% commit ba1862cb39aa638c0dbc700f8f310e3e7dd12dcf
% blob 66484c0901f70ff484ff0b0c6996f03c67edb234
% CTAQ-D1 (a): port for the CT-A run at Paper C's geometry.
% No line below this block differs from the blob.
% PORT-PROVENANCE END
function [g, H, C] = wb_nf_music_grad(U, W, Tw, alpha_model, m_bar, ...
                                      omega, kappa, OPTG)
%WB_NF_MUSIC_GRAD  Sample-domain cost, gradient and Gauss-Newton Hessian of
%   the incoherent compressed-domain MUSIC null residual (T-02, block E2).
%
%   [G, H, C] = WB_NF_MUSIC_GRAD(U, W, TW, ALPHA_MODEL, M_BAR, OMEGA, KAPPA)
%   [G, H, C] = WB_NF_MUSIC_GRAD(..., OPTG)
%
%   Objective (spec v2 Sec. 3):
%
%     Chat(eta) = (1/K_s) sum_k [ 1 - |uhat_k^H u_k(eta)|^2 ],   eta = (omega,kappa)
%     b_k(eta)  = T_W W^H a_k(eta),   n_k = ||b_k||,   u_k = b_k/n_k
%     a_k(eta)  = exp(+1j*al_k*omega*m_bar - 1j*al_k*kappa*m_bar.^2)
%
%   Derivatives (spec v2 Sec. 6.3), with P_k = I - uhat_k uhat_k^H and
%   t_{k,i} = (1/n_k)(I - u_k u_k^H) db_k/d eta_i the TANGENT derivative:
%
%     g_i  = (2/K_s) sum_k Re{ t_{k,i}^H P_k u_k }
%     H_ij = (2/K_s) sum_k Re{ t_{k,i}^H P_k t_{k,j} }      (Gauss-Newton)
%
%     db_k/domega = T_W W^H ( +1j*al_k*diag(m_bar)      * a_k )
%     db_k/dkappa = T_W W^H ( -1j*al_k*diag(m_bar.^2)   * a_k )
%
%   INPUTS
%     U            N_RF x K_s unit-norm sample principal eigenvectors.
%     W            M x N_RF analog combiner.
%     Tw           N_RF x N_RF whitener, T_W G_W T_W^H = I.
%     alpha_model  K_s x 1 subcarrier frequency ratios USED BY THE MODEL.
%                  A scalar is expanded. This is deliberately distinct from
%                  the alpha used to build the data (spec v2 Sec. 9.2, item
%                  1): the R4-8 mismatched run sets alpha_model = 1.
%     m_bar        M x 1 centred element index.
%     omega,kappa  scalar evaluation point.
%     OPTG         optional struct:
%                    .tangent   (default true)  use tangent derivatives.
%                                false selects the norm-variation-omitted
%                                variant, for G-E2b ONLY. Never production.
%                    .cost_only (default false) skip the derivative blocks
%                                and return g = [], H = []. Used by Armijo
%                                backtracking, by the G-E2 central
%                                differences and by the G-E6 terminal
%                                certificate. Saves about two thirds of the
%                                per-k work.
%
%   OUTPUTS
%     g   2 x 1 gradient, or [] when OPTG.cost_only.
%     H   2 x 2 symmetric PSD Gauss-Newton Hessian, or [].
%     C   scalar Chat(eta). Always returned.
%
%   NOTE ON THE cost_only FLAG. The natural nargout branch is unavailable
%   because the signature frozen in spec v2 Sec. 9.1 puts C last, so
%   nargout == 1 returns g rather than C. An explicit flag is used instead.
%   No other deviation from the frozen signature.
%
%   RELATION TO loc_hessian (stageD_ctg_o2_o3_s0c.m line 1127). At the
%   population truth uhat_k equals u_k up to a phase, so P_k = I - u_k u_k^H,
%   the projector is idempotent, and t_i^H P_k t_j collapses to
%   d_i^H P_k d_j / n_k^2. This function then reproduces loc_hessian line for
%   line (spec v2 Sec. 6.4). loc_hessian is NOT superseded: it remains the
%   population anchor for G-E2.
%
%   Provenance class P-RUN. Base MATLAB only. 7-bit ASCII.
%
%   See also WB_NF_GEN_EIGVECS, WB_NF_MUSIC_EST.

% -------------------------------------------------------------------------
% Options and input conditioning
% -------------------------------------------------------------------------
if nargin < 8 || isempty(OPTG)
    OPTG = struct();
end
if ~isfield(OPTG, 'tangent');   OPTG.tangent   = true;  end
if ~isfield(OPTG, 'cost_only'); OPTG.cost_only = false; end

m_bar = m_bar(:);
mb2   = m_bar.^2;
Ks    = size(U, 2);

if isscalar(alpha_model)
    alpha_model = repmat(alpha_model, Ks, 1);
end
alpha_model = alpha_model(:);
assert(numel(alpha_model) == Ks, ...
    'wb_nf_music_grad: alpha_model must be scalar or K_s long');
assert(size(U,1) == size(W,2), ...
    'wb_nf_music_grad: U must be N_RF x K_s');
assert(isscalar(omega) && isscalar(kappa), ...
    'wb_nf_music_grad: omega and kappa must be scalars');

% -------------------------------------------------------------------------
% Accumulation over subcarriers
% -------------------------------------------------------------------------
Csum = 0;
g    = zeros(2,1);
H    = zeros(2);

for k = 1:Ks
    al = alpha_model(k);

    a  = exp(1j*al*omega*m_bar - 1j*al*kappa*mb2);
    b  = Tw*(W'*a);
    n2 = real(b'*b);
    n  = sqrt(n2);
    u  = b/n;

    uh = U(:,k);
    ip = uh'*u;

    Csum = Csum + (1 - abs(ip)^2);

    if OPTG.cost_only
        continue;
    end

    d1 = Tw*(W'*( 1j*al*(m_bar .* a)));       % db/domega
    d2 = Tw*(W'*(-1j*al*(mb2   .* a)));       % db/dkappa
    D  = [d1, d2];

    if OPTG.tangent
        T = (D - u*(u'*D))/n;                 % tangent, spec 6.3
    else
        T = D/n;                              % G-E2b variant ONLY
    end

    Pu = u - uh*(uh'*u);                      % P_k u_k
    PT = T - uh*(uh'*T);                      % P_k t_k

    g = g + 2*real(T'*Pu);
    H = H + 2*real(T'*PT);
end

C = Csum/Ks;

if OPTG.cost_only
    g = [];
    H = [];
    return;
end

g = g/Ks;
H = (H + H')/(2*Ks);

end
