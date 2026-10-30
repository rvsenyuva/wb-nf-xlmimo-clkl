% PORT-PROVENANCE BEGIN
% Copied without change from repository rvsenyuva/wb-nf-crb-globecom26
% path matlab/ctA1/wb_nf_gen_eigvecs.m
% commit ba1862cb39aa638c0dbc700f8f310e3e7dd12dcf
% blob 2aaa5e3643a6df31509381425389bfed31a8f1d8
% CTAQ-D1 (a): port for the CT-A run at Paper C's geometry.
% No line below this block differs from the blob.
% PORT-PROVENANCE END
function [U, meta] = wb_nf_gen_eigvecs(b_data, p, N0, N_snap, stream)
%WB_NF_GEN_EIGVECS  Sample principal eigenvectors for the CT-A1 validation
%   estimator (GLOBECOM 2026 camera-ready, task T-02, block E2).
%
%   [U, META] = WB_NF_GEN_EIGVECS(B_DATA, P, N0, N_SNAP, STREAM)
%
%   Draws one Monte Carlo realisation of the whitened compressed snapshot
%   matrix at every subcarrier and returns the sample principal eigenvector
%   of each per-subcarrier sample covariance.
%
%   INPUTS
%     b_data  N_RF x K_s whitened compressed steering matrix at the TRUE
%             parameters, built on the DATA subcarrier grid, i.e.
%             b_data = loc_bvec(W, Tw, al_DATA, m_bar, mb2, omega_0, kappa_0).
%             The data grid is deliberately not an argument here: the
%             estimator-side alpha_model is separate (spec v2 Sec. 9.2).
%     p       scalar path power.
%     N0      scalar noise power. N0 == 0 selects the noise-free E3 path.
%     N_snap  number of snapshots.
%     stream  RandStream object. Required when N0 > 0; ignored when N0 == 0.
%             An explicit stream (not the global one) is mandatory because
%             controls C-1, C-3 and C-4 and the CT-A2 matched/mismatched
%             pair all require the same realisation under a changed
%             processor.
%
%   OUTPUTS
%     U       N_RF x K_s unit-norm sample principal eigenvectors. The phase
%             gauge of each column is whatever EIG returns; the cost
%             1 - |uhat^H u|^2 is invariant to it.
%     meta    struct of convention diagnostics; see below.
%
%   PINNED CONVENTIONS (spec v2 Sec. 9.4) -- do not "fix" these:
%     (2) Rhat = (Y*Y')/N_snap with NO sample-mean removal. Mean removal
%         would give 63 effective degrees of freedom, a 10log10(64/63) =
%         0.0684 dB shift, comparable to the entire predicted high-SNR gap.
%     (3) s  = (randn + 1j*randn)/sqrt(2)
%         nz = (randn + 1j*randn)*sqrt(N0/2)
%         A missing sqrt(2) is a 3 dB SNR error that survives most
%         structural checks. meta.noise_power_rel_err and meta.gamma exist
%         to catch it empirically (expect gamma = 183.1 at r = 5 m,
%         SNR = 10 dB, combiner seed 42).
%
%   The draw order is s first, then nz, at each k in ascending k. This
%   matches loc_margin_mc in stageD_ctg_o2_o3_s0c.m exactly, so a given
%   stream state reproduces the Stage D realisations.
%
%   Provenance class P-RUN. Base MATLAB only. 7-bit ASCII.
%
%   See also WB_NF_MUSIC_GRAD, WB_NF_MUSIC_EST.

% -------------------------------------------------------------------------
% Input validation
% -------------------------------------------------------------------------
if nargin < 5
    stream = [];
end
assert(ismatrix(b_data) && ~isempty(b_data), ...
       'wb_nf_gen_eigvecs: b_data must be a non-empty N_RF x K_s matrix');
assert(isscalar(p) && isreal(p) && p > 0, ...
       'wb_nf_gen_eigvecs: p must be a positive real scalar');
assert(isscalar(N0) && isreal(N0) && N0 >= 0, ...
       'wb_nf_gen_eigvecs: N0 must be a nonnegative real scalar');
assert(isscalar(N_snap) && N_snap == round(N_snap) && N_snap >= 1, ...
       'wb_nf_gen_eigvecs: N_snap must be a positive integer');

[Nr, Ks] = size(b_data);
noise_free = (N0 == 0);

if noise_free
    % Deterministic unit-modulus row, shared across all k (spec v2 Sec. 10).
    % s*s' = N_snap exactly, so Rhat = p*b*b' and the principal eigenvector
    % is b/||b|| up to a phase, to machine precision.
    s_fixed = exp(1j*2*pi*(0:N_snap-1)/N_snap);
else
    assert(isa(stream, 'RandStream'), ...
        'wb_nf_gen_eigvecs: a RandStream is required when N0 > 0');
end

% -------------------------------------------------------------------------
% Per-subcarrier generation and eigendecomposition
% -------------------------------------------------------------------------
U   = zeros(Nr, Ks);
ev1 = zeros(1, Ks);
ev2 = zeros(1, Ks);

acc_noise = 0;  n_noise = 0;
acc_sig   = 0;  n_sig   = 0;

for k = 1:Ks
    if noise_free
        s = s_fixed;
        Y = sqrt(p)*b_data(:,k)*s;
    else
        s  = (randn(stream, 1,  N_snap) + 1j*randn(stream, 1,  N_snap))/sqrt(2);
        nz = (randn(stream, Nr, N_snap) + 1j*randn(stream, Nr, N_snap))*sqrt(N0/2);
        Y  = sqrt(p)*b_data(:,k)*s + nz;
        acc_noise = acc_noise + sum(abs(nz(:)).^2);
        n_noise   = n_noise   + numel(nz);
    end
    acc_sig = acc_sig + sum(abs(s).^2);
    n_sig   = n_sig   + numel(s);

    R = (Y*Y')/N_snap;                  % NO sample-mean removal (9.4.2)
    [Ve, De] = eig((R+R')/2);
    dv = real(diag(De));
    [dsrt, ord] = sort(dv, 'descend');

    uh       = Ve(:, ord(1));
    U(:,k)   = uh / norm(uh);
    ev1(k)   = dsrt(1);
    if Nr >= 2
        ev2(k) = dsrt(2);
    else
        ev2(k) = NaN;
    end
end

% -------------------------------------------------------------------------
% Convention diagnostics
% -------------------------------------------------------------------------
bn2 = sum(abs(b_data).^2, 1);

meta = struct();
meta.noise_free   = noise_free;
meta.N_RF         = Nr;
meta.K_s          = Ks;
meta.N_snap       = N_snap;
meta.p            = p;
meta.N0           = N0;

meta.bnorm2       = bn2;                    % ||b_k||^2, expect ~N_RF
meta.bnorm2_mean  = mean(bn2);              % Stage D reports 18.31 vs 16

meta.signal_power_meas = acc_sig/n_sig;     % assertion 9.4(3): expect 1

if noise_free
    meta.gamma                = Inf;
    meta.gamma_k              = inf(1, Ks);
    meta.noise_power_meas     = 0;
    meta.noise_power_rel_err  = 0;
else
    meta.gamma                = p*mean(bn2)/N0;   % expect 183.1 at the anchor
    meta.gamma_k              = p*bn2/N0;
    meta.noise_power_meas     = acc_noise/n_noise;
    meta.noise_power_rel_err  = meta.noise_power_meas/N0 - 1;
end

meta.eig1              = ev1;
meta.eig2              = ev2;
meta.eiggap_abs        = ev1 - ev2;                  % Ledger D-06
meta.eiggap_ratio      = ev1 ./ max(ev2, realmin);
meta.unit_norm_max_err = max(abs(sqrt(sum(abs(U).^2,1)) - 1));

end
