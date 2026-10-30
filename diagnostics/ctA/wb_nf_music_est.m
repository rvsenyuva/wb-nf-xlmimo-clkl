% PORT-PROVENANCE BEGIN
% Copied without change from repository rvsenyuva/wb-nf-crb-globecom26
% path matlab/ctA1/wb_nf_music_est.m
% commit ba1862cb39aa638c0dbc700f8f310e3e7dd12dcf
% blob da1c9ca20991bb8cdad349b0d6dc913d84f00829
% CTAQ-D1 (a): port for the CT-A run at Paper C's geometry.
% No line below this block differs from the blob.
% PORT-PROVENANCE END
function [eta_hat, info] = wb_nf_music_est(U, W, Tw, alpha_model, m_bar, ...
                                           GRID, OPT)
%WB_NF_MUSIC_EST  Incoherent compressed-domain MUSIC estimator for the
%   CT-A1 bound-attainability check (T-02, block E2).
%
%   [ETA_HAT, INFO] = WB_NF_MUSIC_EST(U, W, TW, ALPHA_MODEL, M_BAR, GRID, OPT)
%
%   Three stages (spec v2 Sec. 5):
%     1  detection   Chat on the full coarse grid using K_COARSE subcarriers;
%                    8-neighbourhood local minima; retain the N_CAND lowest.
%     2  selection   re-evaluate those N_CAND at all K_s; single global
%                    argmin -> eta_coarse.
%     3  refinement  exactly ONE Gauss-Newton refinement started at
%                    eta_coarse. No multi-start, no warm start, no post-hoc
%                    candidate comparison (register Sec. 10 firewall).
%
%   INPUTS
%     U            N_RF x K_s x n_trials unit-norm sample principal
%                  eigenvectors from WB_NF_GEN_EIGVECS. A 2-D array is the
%                  n_trials = 1 case.
%                  [Deviation from the spec v2 Sec. 9.1 signature, approved
%                   2026-08-06: the batch dimension is required because the
%                   Sec. 9.3 loop amortises the coarse manifold across
%                   trials. 2-D input behaves exactly as before.]
%     W            M x N_RF analog combiner.
%     Tw           N_RF x N_RF whitener.
%     alpha_model  K_s x 1 model subcarrier ratios, or scalar. Distinct from
%                  the data alpha (spec v2 Sec. 9.2 item 1).
%     m_bar        M x 1 centred element index.
%     GRID  .om_grid  n_om x 1 omega grid (endpoints already dropped)
%           .ka_grid  n_ka x 1 kappa grid
%           .idx_c    1 x K_COARSE detection subcarrier indices
%           .N_CAND   candidates carried to stage 2
%     OPT   .tol_dec   2e-12   scaled Newton decrement tolerance
%           .tol_ml    1e-9    terminal step, in mainlobes
%           .maxit     50
%           .armijo_c1 1e-4    .armijo_fac 0.5    .armijo_max 8
%           .tr_frac   0.5     trust radius, in mainlobes
%           .cond_max  1e8     ill-conditioning guard on cond(H_s)
%           .cert_frac 1e-3    terminal certificate offset, in mainlobes
%           .batch     1000    trials per coarse chunk. CHUNKING IS A
%                              MEMORY GUARD, NOT AN OPTIMISATION: the
%                              dominant coarse cost is trial-independent
%                              and is REPEATED once per chunk (measured
%                              2026-08-06 at 1.34 s per chunk for the S0-C
%                              grid). Keep batch >= n_trials whenever the
%                              working set allows; it is about 80 kB per
%                              trial, so 1000 trials is roughly 80 MB.
%           .npar      0       parfor workers for stages 2-3
%           .s0, .q            required for the range map
%           .eta_true  []      OPTIONAL [omega_0 kappa_0]. Post-hoc CAPTURED
%                              and RANK diagnostics only. Structurally it is
%                              read only AFTER the candidate set and the
%                              estimate are fixed; it cannot reach the
%                              estimate path.
%           .tangent   true    use tangent derivatives in the refinement.
%                              FALSE selects the norm-variation-omitted
%                              variant and exists ONLY for gate G-E2b. It
%                              must be true for every production run and
%                              for every number that reaches a figure.
%           .assert_internal true  run the free stage-2/stage-3 cost
%                              consistency check
%
%   OUTPUTS
%     eta_hat  n_trials x 2, columns [omega_hat, kappa_hat].
%     info     struct of per-trial arrays (Sec. 7.1 record plus diagnostics);
%              see the assembly block at the end of this file.
%
%   RANGE MAP (spec v2 Sec. 3), algebraic, no acos round trip:
%     r_hat = s0*(1 - (omega_hat/q)^2)/kappa_hat
%   Plug-in in BOTH coordinates, which is what makes the full-Jacobian bound
%   g^T C g the correct comparator. Holding theta at truth would silently
%   switch to the known-angle bound (Ledger A-29), 0.0117 dB away at r = 5 m.
%
%   PINNED CONVENTION (spec v2 Sec. 9.4 item 1). The coarse slice uses
%   B = X * Tw.' -- a PLAIN transpose -- because X holds (W^H a)^T as rows.
%   Replacing .' with ' silently breaks the whitening while leaving every
%   structural check intact. DO NOT "FIX" IT.
%
%   NO CRB VALUE ENTERS THIS FUNCTION AT ANY POINT (spec v2 Sec. 7.1).
%
%   Provenance class P-RUN. Base MATLAB only. 7-bit ASCII.
%
%   See also WB_NF_GEN_EIGVECS, WB_NF_MUSIC_GRAD.

% =========================================================================
% 0. Options, validation, derived constants
% =========================================================================
if nargin < 7 || isempty(OPT); OPT = struct(); end
OPT = loc_defaults(OPT);

m_bar = m_bar(:);
mb2   = m_bar.^2;
M     = numel(m_bar);

if ndims(U) == 2 %#ok<ISMAT>
    U = reshape(U, size(U,1), size(U,2), 1);
end
[Nr, Ks, n_tr] = size(U);

assert(size(W,1) == M && size(W,2) == Nr, ...
    'wb_nf_music_est: W must be M x N_RF consistent with m_bar and U');
assert(isequal(size(Tw), [Nr Nr]), 'wb_nf_music_est: Tw must be N_RF x N_RF');

if isscalar(alpha_model); alpha_model = repmat(alpha_model, Ks, 1); end
alpha_model = alpha_model(:);
assert(numel(alpha_model) == Ks, ...
    'wb_nf_music_est: alpha_model must be scalar or K_s long');

om_grid = GRID.om_grid(:);
ka_grid = GRID.ka_grid(:);
idx_c   = GRID.idx_c(:).';
N_CAND  = GRID.N_CAND;
n_om    = numel(om_grid);
n_ka    = numel(ka_grid);
KC      = numel(idx_c);

assert(n_om >= 3 && n_ka >= 1, 'wb_nf_music_est: grid too small');
assert(all(idx_c >= 1 & idx_c <= Ks), 'wb_nf_music_est: idx_c out of range');
assert(N_CAND >= 1 && N_CAND <= n_om*n_ka, 'wb_nf_music_est: bad N_CAND');
assert(~isempty(OPT.s0) && ~isempty(OPT.q), ...
    'wb_nf_music_est: OPT.s0 and OPT.q are required for the range map');

CONST = struct();
CONST.om_lo = om_grid(1);    CONST.om_hi = om_grid(end);
CONST.ka_lo = ka_grid(1);    CONST.ka_hi = ka_grid(end);
CONST.ml_om = 2*pi/M;                 % mainlobe width in omega
CONST.ml_ka = pi/max(mb2);            % mainlobe width in kappa

% =========================================================================
% 1. Stage 1 -- coarse detection, batched over trials
% =========================================================================
% The coarse manifold depends only on (W, alpha_model, grid): never on the
% trial, the SNR, or the true parameters. It is therefore built once per
% (kappa, subcarrier) and reduced against the whole trial batch in one
% matmul (spec v2 Sec. 9.3). Naive per-trial rebuilding costs about 113x
% more and puts CT-A2 out of reach.
%
% Under the mismatched model (alpha_model = 1) only the MANIFOLD BUILD
% collapses to a single subcarrier; the reduction still runs KC terms in
% stage 1 and K_s in stage 2, against distinct uhat_k.

al_c = alpha_model(idx_c);
[al_u, ~, map_c] = unique(al_c);
Ecache = cell(numel(al_u), 1);
for iu = 1:numel(al_u)
    Ecache{iu} = exp(1j*al_u(iu)*(om_grid*m_bar.'));      % n_om x M
end
Wc  = conj(W);
TwT = Tw.';                    % PINNED plain transpose, spec 9.4 item 1

candQ = inf(N_CAND, n_tr);
candI = zeros(N_CAND, n_tr);
candJ = zeros(N_CAND, n_tr);

n_chunk = ceil(n_tr/OPT.batch);
for ic = 1:n_chunk
    tr = ((ic-1)*OPT.batch + 1) : min(ic*OPT.batch, n_tr);
    [candQ(:,tr), candI(:,tr), candJ(:,tr)] = loc_stage1( ...
        U(:,idx_c,tr), Ecache, map_c, al_c, ka_grid, Wc, TwT, mb2, ...
        n_om, n_ka, KC, N_CAND);
end

% If the coarse surface yields fewer than N_CAND local minima, candI holds
% zeros and stage 2 would fail with an opaque index error. This cannot occur
% on the S0-C grid (roughly 9,188 effectively independent cells against
% N_CAND = 200) but is guarded so a grid change fails legibly.
assert(all(candI(:) > 0), ...
    ['wb_nf_music_est: fewer than N_CAND local minima found on the ' ...
     'coarse grid. Lower GRID.N_CAND or widen the grid.']);

% =========================================================================
% 2-3. Stage 2 selection and stage 3 refinement, per trial
% =========================================================================
eta_hat = zeros(n_tr, 2);
etac  = zeros(n_tr, 2);
Cco   = zeros(n_tr,1);   Cfi   = zeros(n_tr,1);
rha   = zeros(n_tr,1);   tha   = zeros(n_tr,1);
nit   = zeros(n_tr,1);   lam2  = zeros(n_tr,1);
ginf  = zeros(n_tr,1);   cmin  = zeros(n_tr,1);   cmax = zeros(n_tr,1);
bnd   = false(n_tr,1);   ill   = false(n_tr,1);   nst  = false(n_tr,1);
accf  = zeros(n_tr,1);   stpml = zeros(n_tr,1);
certd = zeros(n_tr,4);   certp = false(n_tr,1);
cmm   = zeros(n_tr,1);
srsn  = cell(n_tr,1);    vdct  = cell(n_tr,1);

s0v = OPT.s0;  qv = OPT.q;
npar = OPT.npar;

% Materialise the candidate coordinates so they enter the loop as SLICED
% variables. Indexing om_grid/ka_grid inside the loop would broadcast both
% grids to every worker (checkcode, 2026-08-06).
om_cand = om_grid(candI);
ka_cand = ka_grid(candJ);

parfor (it = 1:n_tr, npar)
    Uk = U(:,:,it);

    % ---- Stage 2: full-K_s re-evaluation of the retained candidates ----
    om_c = om_cand(:,it);
    ka_c = ka_cand(:,it);
    Qf   = loc_stage2(Uk, W, Tw, alpha_model, m_bar, mb2, om_c, ka_c);
    [Qc0, ib] = min(Qf);
    eta0 = [om_c(ib); ka_c(ib)];

    % ---- Stage 3: exactly one Gauss-Newton refinement ------------------
    out = loc_refine(eta0, Uk, W, Tw, alpha_model, m_bar, CONST, OPT);

    eta_hat(it,:) = out.eta.';
    etac(it,:)    = eta0.';
    Cco(it)   = Qc0;          Cfi(it)   = out.C;
    nit(it)   = out.n_iter;   lam2(it)  = out.lambda2;
    ginf(it)  = out.g_inf;    cmin(it)  = out.cond_min;
    cmax(it)  = out.cond_max; bnd(it)   = out.boundary;
    ill(it)   = out.illcond;  nst(it)   = out.nonstationary;
    accf(it)  = out.acc_frac; stpml(it) = out.step_ml;
    certd(it,:) = out.cert;   certp(it) = out.cert_pass;
    srsn{it}  = out.stop_reason;
    vdct{it}  = out.verdict;

    % Free internal check: stage 2's cost at eta_coarse must equal the
    % refinement's cost at its own starting point.
    cmm(it) = abs(Qc0 - out.C_start);

    % Range map, plug-in in BOTH coordinates
    rha(it) = s0v*(1 - (out.eta(1)/qv)^2)/out.eta(2);
    tha(it) = acos(max(-1, min(1, out.eta(1)/qv)));
end

if OPT.assert_internal
    assert(max(cmm) < 1e-10, ...
        'wb_nf_music_est: stage-2 and stage-3 costs disagree at eta_coarse');
end

% =========================================================================
% 4. Post-hoc CAPTURED and RANK (truth is read only here)
% =========================================================================
captured = false(n_tr,1);
rankv    = nan(n_tr,1);
if ~isempty(OPT.eta_true)
    [~, i0] = min(abs(om_grid - OPT.eta_true(1)));
    [~, j0] = min(abs(ka_grid - OPT.eta_true(2)));
    for it = 1:n_tr
        hit = (abs(candI(:,it) - i0) <= 3) & (abs(candJ(:,it) - j0) <= 3);
        if any(hit)
            captured(it) = true;
            rankv(it)    = find(hit, 1, 'first');   % candQ is sorted ascending
        end
    end
end

% =========================================================================
% 5. Assemble info
% =========================================================================
info = struct();
info.eta_coarse      = etac;      % stage-2 argmin, i.e. the refinement start
info.Chat_coarse     = Cco;
info.Chat_final      = Cfi;
info.r_hat           = rha;
info.theta_hat_deg   = tha*180/pi;
info.n_iter          = nit;
info.lambda2         = lam2;
info.g_inf           = ginf;
info.cond_Hs_min     = cmin;
info.cond_Hs_max     = cmax;
info.boundary        = bnd;
info.illcond         = ill;
info.nonstationary   = nst;
info.acc_frac        = accf;
info.step_ml         = stpml;
info.cert_delta      = certd;
info.cert_pass       = certp;
info.stop_reason     = srsn;
info.verdict         = vdct;
info.captured        = captured;
info.rank            = rankv;
info.coarse_cost_mismatch = cmm;
info.cand_Q          = candQ;
info.cand_iom        = candI;
info.cand_ika        = candJ;
info.CONST           = CONST;
info.OPT             = OPT;

end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function OPT = loc_defaults(OPT)
d = { 'tol_dec',2e-12; 'tol_ml',1e-9; 'maxit',50; ...
      'armijo_c1',1e-4; 'armijo_fac',0.5; 'armijo_max',8; ...
      'tr_frac',0.5; 'cond_max',1e8; 'cert_frac',1e-3; ...
      'batch',1000; 'npar',0; 's0',[]; 'q',[]; 'eta_true',[]; ...
      'tangent',true; ...
      'assert_internal',true };
for i = 1:size(d,1)
    if ~isfield(OPT, d{i,1}); OPT.(d{i,1}) = d{i,2}; end
end
end

% -------------------------------------------------------------------------
function [bQ, bI, bJ] = loc_stage1(Uc, Ecache, map_c, al_c, ka_grid, ...
                                   Wc, TwT, mb2, n_om, n_ka, KC, N_CAND)
% Coarse detection for one trial chunk.
%
% Memory strategy: the full 1534 x 288 surface is never materialised. An
% 8-neighbourhood local-minimum test needs only three consecutive kappa
% rows, and retaining the N_CAND lowest local minima seen so far is exact
% for any buffer of size >= N_CAND. Working set is three n_om x nt rows plus
% the buffer, roughly 4 MB at nt = 200, against 353 MB for the full surface.

nt = size(Uc,3);
Nr = size(Uc,1);

% conj(uhat) for each detection subcarrier, N_RF x nt
Ucj = cell(KC,1);
for j = 1:KC
    Ucj{j} = conj(reshape(Uc(:,j,:), Nr, nt));
end

bQ = inf(N_CAND, nt);
bI = zeros(N_CAND, nt);
bJ = zeros(N_CAND, nt);

rowIdx = repmat((1:n_om).', 1, nt);

Qm = [];  Qc = [];
for ik = 1:n_ka+1
    if ik <= n_ka
        Qn = loc_qrow(ik);
    else
        Qn = [];
    end

    if ~isempty(Qc)
        jk  = ik - 1;                     % the row now fully bracketed
        Qup = Qm;  if isempty(Qup); Qup = Qc; end    % replicate at edge
        Qdn = Qn;  if isempty(Qdn); Qdn = Qc; end

        isMin = loc_localmin(Qc, Qup, Qdn);
        newQ  = Qc;
        newQ(~isMin) = Inf;

        cQ = [bQ; newQ];
        cI = [bI; rowIdx];
        cJ = [bJ; jk*ones(n_om, nt)];
        [sQ, ord] = sort(cQ, 1);
        keep = ord(1:N_CAND, :);
        lin  = keep + (0:nt-1)*size(cQ,1);
        bQ = sQ(1:N_CAND, :);
        bI = cI(lin);
        bJ = cJ(lin);
    end

    Qm = Qc;  Qc = Qn;
end

    function Q = loc_qrow(ikk)
        Q = zeros(n_om, nt);
        for jj = 1:KC
            ch = exp(-1j*al_c(jj)*ka_grid(ikk)*mb2);      % M x 1
            X  = Ecache{map_c(jj)} * (Wc .* ch);          % n_om x N_RF
            B  = X * TwT;                                 % PINNED .'
            nb = sum(abs(B).^2, 2);
            CC = B * Ucj{jj};                             % n_om x nt
            Q  = Q + (1 - abs(CC).^2 ./ nb);
        end
        Q = Q/KC;
    end
end

% -------------------------------------------------------------------------
function isMin = loc_localmin(Qc, Qup, Qdn)
% 8-neighbourhood minimum test on the centre row, replicate padding in
% omega. Matches loc_landscape / padarray_local in stageD_ctg_o2_o3_s0c.m,
% including the <= convention on plateaus.
up = [Qc(1,:);     Qc(1:end-1,:)];
dn = [Qc(2:end,:); Qc(end,:)];
Au = [Qup(1,:);     Qup(1:end-1,:)];
Ad = [Qup(2:end,:); Qup(end,:)];
Bu = [Qdn(1,:);     Qdn(1:end-1,:)];
Bd = [Qdn(2:end,:); Qdn(end,:)];
isMin = (Qc <= up) & (Qc <= dn) & ...
        (Qc <= Qup) & (Qc <= Au) & (Qc <= Ad) & ...
        (Qc <= Qdn) & (Qc <= Bu) & (Qc <= Bd);
end

% -------------------------------------------------------------------------
function Qf = loc_stage2(Uk, W, Tw, alpha_model, m_bar, mb2, om_c, ka_c)
% Full-K_s re-evaluation of the retained candidates, vectorised over the
% candidate set. Cost N_CAND x K_s x M x N_RF MACs per trial.
Ks  = size(Uk,2);
Nc  = numel(om_c);
Phi = m_bar*om_c(:).' - mb2*ka_c(:).';        % M x Nc, real
Qf  = zeros(Nc,1);
for k = 1:Ks
    A   = exp(1j*alpha_model(k)*Phi);         % M x Nc
    Bc  = Tw*(W'*A);                          % N_RF x Nc
    ipv = Uk(:,k)'*Bc;                        % 1 x Nc
    nb  = sum(abs(Bc).^2, 1);
    Qf  = Qf + (1 - (abs(ipv).^2 ./ nb)).';
end
Qf = Qf/Ks;
end

% -------------------------------------------------------------------------
function out = loc_refine(eta0, Uk, W, Tw, alpha_model, m_bar, CONST, OPT)
% Exactly ONE Gauss-Newton refinement, spec v2 Sections 6-8.

OPTC = struct('tangent', true, 'cost_only', true);   % cost is tangent-free
OPTD = struct('tangent', OPT.tangent, 'cost_only', false);
Dv   = [CONST.ml_om; CONST.ml_ka];       % D_solve = D_stop, mainlobe metric

eta = eta0(:);
[~,~,C_start] = wb_nf_music_grad(Uk, W, Tw, alpha_model, m_bar, ...
                                 eta(1), eta(2), OPTC);
Ccur = C_start;

n_acc = 0;  n_try = 0;
illcond = false;  stall = false;  converged = false;
cond_min = Inf;   cond_max = 0;
lambda2 = NaN;    g_inf = NaN;   step_ml = Inf;
stop_reason = 'MAXIT';
n_iter = 0;

for it = 1:OPT.maxit+1
    [g, H, Ccur] = wb_nf_music_grad(Uk, W, Tw, alpha_model, m_bar, ...
                                    eta(1), eta(2), OPTD);
    Hs = (Dv*Dv.').*H;     Hs = (Hs + Hs.')/2;      % D H D, D diagonal
    gs = Dv.*g;
    g_inf = max(abs(gs));

    [Vh, Lh] = eig(Hs);
    lam = real(diag(Lh));
    if min(lam) <= 0
        cnd = Inf;
    else
        cnd = max(lam)/min(lam);
    end
    cond_min = min(cond_min, cnd);
    cond_max = max(cond_max, cnd);

    if ~isfinite(cnd) || cnd > OPT.cond_max
        % Step along the well-conditioned eigendirection only. This is a
        % 2x2 OPTIMISER Hessian; no covariance is loaded anywhere, so R3-1
        % item 5 is untouched. Because this branch can produce a small
        % NON-stationary step, the stopping rule is written on lambda2.
        illcond = true;
        [~, imx] = max(lam);
        v  = Vh(:,imx);
        gv = v'*gs;
        ps = -(gv/lam(imx))*v;
        lambda2 = real(gv'*gv)/lam(imx);
    else
        ps = -(Hs\gs);
        lambda2 = -real(gs.'*ps);
    end

    % ---- Stopping test, evaluated BEFORE stepping (spec 7.1) ----------
    % step_ml is the terminal ACCEPTED step in mainlobes from the previous
    % iteration, which is the quantity C-2 partitions.
    if lambda2 <= OPT.tol_dec || (step_ml <= OPT.tol_ml && lambda2 <= 1e-8)
        converged   = true;
        stop_reason = 'CONVERGED';
        break;
    end
    if it > OPT.maxit
        stop_reason = 'MAXIT';
        break;
    end

    delta = Dv.*ps;

    % ---- Trust region: UNIFORM scaling, never per-coordinate clipping --
    % Component-wise clipping can flip the sign of g'*delta and destroy
    % descent.
    st = min([1, ...
              OPT.tr_frac*CONST.ml_om/max(abs(delta(1)), realmin), ...
              OPT.tr_frac*CONST.ml_ka/max(abs(delta(2)), realmin)]);
    delta = st*delta;

    % ---- Armijo backtracking, evaluated at the PROJECTED point ---------
    s = 1;  accepted = false;  cand = eta;  Ct = Ccur;
    for ih = 0:OPT.armijo_max
        cand = loc_project(eta + s*delta, CONST);
        [~,~,Ct] = wb_nf_music_grad(Uk, W, Tw, alpha_model, m_bar, ...
                                    cand(1), cand(2), OPTC);
        n_try = n_try + 1;
        if Ct <= Ccur + OPT.armijo_c1*(g.'*(cand - eta))
            accepted = true;
            break;
        end
        s = s*OPT.armijo_fac;
    end

    if ~accepted
        stall = true;
        stop_reason = 'STALL';
        break;
    end

    step    = cand - eta;
    step_ml = max(abs(step./Dv));
    eta     = cand;
    Ccur    = Ct;
    n_acc   = n_acc + 1;
    n_iter  = n_iter + 1;
end

% ---- Terminal certificate, model-free (spec 7.1) ----------------------
% Uses ONLY the objective, so unlike lambda2 it is immune to a common-mode
% error in the analytic gradient.
cert = zeros(1,4);
pert = OPT.cert_frac*[CONST.ml_om 0; -CONST.ml_om 0; ...
                      0 CONST.ml_ka; 0 -CONST.ml_ka];
for ip = 1:4
    e = eta + pert(ip,:).';
    [~,~,Cp] = wb_nf_music_grad(Uk, W, Tw, alpha_model, m_bar, ...
                                e(1), e(2), OPTC);
    cert(ip) = Cp - Ccur;
end
cert_pass = all(cert > 0);

% ---- Boundary status ---------------------------------------------------
tol_b = 1e-12*max(1, abs([CONST.om_hi CONST.ka_hi]));
boundary = (abs(eta(1)-CONST.om_lo) <= tol_b(1)) || ...
           (abs(eta(1)-CONST.om_hi) <= tol_b(1)) || ...
           (abs(eta(2)-CONST.ka_lo) <= tol_b(2)) || ...
           (abs(eta(2)-CONST.ka_hi) <= tol_b(2));

% ---- Verdict, precedence BOUNDARY > STALL > MAXIT > ILLCOND > OK ------
% NONSTATIONARY is carried as a separate boolean rather than as a competing
% label, because spec 8.2 defines it as "termination for any other reason",
% which STALL and MAXIT both satisfy. ILLCOND is likewise a flag that can
% co-occur with a converged solve.
nonstationary = ~converged;
if boundary
    verdict = 'BOUNDARY';
elseif stall
    verdict = 'STALL';
elseif nonstationary
    verdict = 'MAXIT';
elseif illcond
    verdict = 'ILLCOND';
else
    verdict = 'OK';
end

out = struct();
out.eta           = eta;
out.C             = Ccur;
out.C_start       = C_start;
out.n_iter        = n_iter;
out.lambda2       = lambda2;
out.g_inf         = g_inf;
out.cond_min      = cond_min;
out.cond_max      = cond_max;
out.boundary      = boundary;
out.illcond       = illcond;
out.nonstationary = nonstationary;
out.acc_frac      = n_acc/max(n_try,1);
out.step_ml       = step_ml;
out.cert          = cert;
out.cert_pass     = cert_pass;
out.stop_reason   = stop_reason;
out.verdict       = verdict;
end

% -------------------------------------------------------------------------
function e = loc_project(e, CONST)
e(1) = min(max(e(1), CONST.om_lo), CONST.om_hi);
e(2) = min(max(e(2), CONST.ka_lo), CONST.ka_hi);
end
