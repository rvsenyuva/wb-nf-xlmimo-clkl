function P = setup_production_P(varargin)
%SETUP_PRODUCTION_P  Production parameter struct for Paper C Phase 3.
%
%  P = setup_production_P()
%  P = setup_production_P('snr')
%  P = setup_production_P('bandwidth')
%  P = setup_production_P('convergence')
%
%  Each sweep_mode may carry a second range_mode argument:
%
%  P = setup_production_P('snr',        'full')       % r_hi_fac=1.00  [DEFAULT]
%  P = setup_production_P('snr',        'restricted') % r_hi_fac=0.50
%  P = setup_production_P('bandwidth',  'restricted') % r_hi_fac=0.50
%  P = setup_production_P('convergence','restricted') % r_hi_fac=0.50
%
%  RANGE MODES
%  -----------
%  'full'       r_hi_fac=1.00:  full Rayleigh range [1.06, 21.26] m.
%               Use for CRB scripts (Tasks 12.4-12.6) and failure-mode
%               documentation in Paper C Sec.V.
%               B4 performance vs-CRB comparison NOT valid here (near-
%               Rayleigh Fresnel curvature too small for BPD warm-start).
%
%  'restricted' r_hi_fac=0.50:  strong near-field regime [1.06, 10.63] m.
%               At r=0.5*r_RD, theta=40 deg: kappa*m_bar_max^2 ~ 0.81 rad,
%               sufficient for reliable BPD range detection and WB-CL-KL
%               convergence.  Use for Figs 8-9 (RMSE vs SNR, NMSE vs BW).
%
%  DESIGN DECISION (Phase 3, May 5 2026)
%  --------------------------------------
%  Production run 1 (r_hi_fac=1.00, G_r=64) showed B4 B4/CRB gap of
%  31-96x at SNR=10-20 dB with conv_pct=24-27%.  Root cause: kappa at
%  r_RD is 0.00016 rad/m, giving kappa*m_bar_max^2 = 0.16 rad at r_hi.
%  BPD power map is physically flat in range at this curvature -- a finer
%  grid cannot recover what physics does not provide.  'restricted' mode
%  limits r to the regime where Fresnel curvature enables discrimination.
%
%  BPD DIAGNOSTIC CSV FIX (v2)
%  ----------------------------
%  P.bpd_write_csv = false suppresses the per-realisation CSV write inside
%  bpd_baseline.m that produced ~200 console messages per sweep point in
%  production run 1.  Set P.bpd_write_csv = true after calling this
%  function only when debugging a single realisation interactively.
%
%  GRID CHANGES FROM v1
%  --------------------
%  G_r: 64 -> 128.  Halves range-grid quantisation error in the detectable
%  regime (r < 0.5*r_RD) at 2x BPD cost.  G_r > 128 gives no additional
%  gain because physics (not grid density) limits range discrimination.
%  G_theta: 128 (unchanged -- angular resolution already adequate).
%
%  All fields comply with:
%    - Paper C Phase 3 Starter Prompt Sec.4 (production P-struct spec)
%    - Paper_C_Roadmap_v13.md Sec.9 and Sec.13
%    - Lessons L29-L36 (Roadmap v13 Sec.15)
%    - Parameter alias rule: P.c0/P.c and P.lambda_c/P.lambda both set
%
%  USAGE EXAMPLES
%  --------------
%  % Fig 8 -- RMSE vs SNR (restricted regime, publication run):
%    P = setup_production_P('snr', 'restricted');
%    run_monte_carlo_paperC(P, 'snr', [-5 0 5 10 15 20], 200, ...
%        'results/snr_sweep_restricted');
%
%  % Fig 9 -- NMSE vs BW (restricted regime):
%    P = setup_production_P('bandwidth', 'restricted');
%    run_monte_carlo_paperC(P, 'bandwidth', [100 200 400 600 800]*1e6, 200, ...
%        'results/bw_sweep_restricted');
%
%  % Fig 10 -- Convergence (r=5.0 m fixed < 0.5*r_RD; either mode valid):
%    P = setup_production_P('convergence', 'restricted');
%    run_monte_carlo_paperC(P, 'convergence', [0 5 10 15], 50, ...
%        'results/convergence_restricted');
%
%  % CRB scripts (Tasks 12.4-12.6) -- always full range:
%    P = setup_production_P('snr', 'full');
%
%  Author   : R. V. Senyuva (Maltepe University)
%  Date     : May 2026  (v2: CSV fix, restricted mode, G_r=128)
%  Ref      : Paper C Phase 3, Tasks 12.1-12.6 (Roadmap v13 Sec.9).
%  Called by: run_monte_carlo_paperC.m, wb_crb_multipath_sweep.m,
%             wb_crb_geometric_sweep.m

% =========================================================================
%  INPUT HANDLING
% =========================================================================
if nargin == 0
    sweep_mode = 'snr';
    range_mode = 'full';
elseif nargin == 1
    sweep_mode = lower(char(varargin{1}));
    range_mode = 'full';
else
    sweep_mode = lower(char(varargin{1}));
    range_mode = lower(char(varargin{2}));
end

valid_sweep = {'snr', 'bandwidth', 'convergence'};
valid_range = {'full', 'restricted'};
assert(any(strcmp(sweep_mode, valid_sweep)), ...
    'setup_production_P: sweep_mode must be ''snr'', ''bandwidth'', or ''convergence''.');
assert(any(strcmp(range_mode, valid_range)), ...
    'setup_production_P: range_mode must be ''full'' or ''restricted''.');

% =========================================================================
%  SECTION 1 -- Physical constants and carrier frequency
% =========================================================================
P.c  = 3e8;          % speed of light [m/s]
P.c0 = P.c;          % alias (parameter alias rule)
P.fc = 28e9;         % carrier frequency [Hz]  (28 GHz mmWave)

% =========================================================================
%  SECTION 2 -- Wavelength and element spacing
% =========================================================================
P.lambda   = P.c / P.fc;       % carrier wavelength [m]  ~ 0.010714 m
P.lambda_c = P.lambda;         % alias (parameter alias rule)
P.d_ant    = P.lambda_c / 2;   % half-wavelength element spacing [m]

% =========================================================================
%  SECTION 3 -- Array dimensions and path count
% =========================================================================
P.M    = 64;   % number of ULA elements
P.N_RF = 8;    % number of RF chains
P.N    = 64;   % snapshots per subcarrier
P.d    = 1;    % paths (single-path default; override externally for Task 12.4)

% =========================================================================
%  SECTION 4 -- Rayleigh distance and range bounds
% =========================================================================
D_ap   = (P.M - 1) * P.d_ant;           % aperture length [m]  = 0.3375 m
P.r_RD = 2 * D_ap^2 / P.lambda_c;       % Rayleigh distance [m] ~ 21.26 m

P.r_lo_fac = 0.05;   % inner range = 0.05 * r_RD ~ 1.063 m (both modes)

switch range_mode
    case 'full'
        P.r_hi_fac = 1.00;   % r_hi ~ 21.26 m
    case 'restricted'
        P.r_hi_fac = 0.50;   % r_hi ~ 10.63 m
end

P.u_margin = 2.0;
P.u_min    = 1 / (P.r_hi_fac * P.r_RD * P.u_margin);
P.u_max    = 1 / (P.r_lo_fac * P.r_RD / P.u_margin);

% =========================================================================
%  SECTION 5 -- Angular bounds
% =========================================================================
P.theta_lo = 20 * pi/180;   % [rad]  20 deg
P.theta_hi = 60 * pi/180;   % [rad]  60 deg

% =========================================================================
%  SECTION 6 -- Wideband subcarrier parameters (B = 400 MHz default)
%  Lesson L29: Delta_f = 25 MHz (NOT 5G NR 120 kHz numerology).
% =========================================================================
P.B       = 400e6;                   % bandwidth [Hz]
P.K_s     = 16;                      % number of active subcarriers
P.Delta_f = P.B / P.K_s;            % subcarrier spacing [Hz] = 25 MHz
P.K       = P.K_s;                  % alias

k_idx         = (-(P.K_s/2) : (P.K_s/2 - 1)).';   % K_s x 1 integer indices
P.k_indices   = k_idx;
P.alpha_k_vec = 1 + k_idx * (P.Delta_f / P.fc);    % K_s x 1  f_k/f_c ratios
% Range: [1 - 7*25e6/28e9, 1 + 8*25e6/28e9] = [0.99286, 1.00625]

% =========================================================================
%  SECTION 7 -- Grid resolutions
%  G_r upgraded 64->128 in v2.  See header for rationale.
% =========================================================================
P.Q_theta    = 256;   % WB-P-SOMP polar angle grid (unchanged)
P.beta_delta = 1.2;   % range oversampling factor for u-domain grid

P.G_theta    = 128;   % BPD angle grid points  (unchanged from v1)
P.G_r        = 128;   % BPD range grid points  (v1: 64 -> v2: 128)
P.d_max      = P.d;   % BPD max path count (= d; update externally if d changes)

% =========================================================================
%  SECTION 8 -- WB-CL-KL optimiser parameters
% =========================================================================
P.lambda_reg = 1e-4;   % L1 regularisation weight
P.max_iter   = 200;    % maximum gradient iterations (Task 11.3 rev3 C1)
P.tol_clkl   = 1e-5;   % convergence tolerance (relative KL change)
P.alpha_p    = 0.5;    % Armijo line-search initial step size
P.ls_beta    = 0.5;    % line-search backtrack factor
P.ls_sigma   = 1e-4;   % Armijo sufficient-decrease coefficient
P.eps_reg    = 1e-3;   % SNR-adaptive regularisation coeff (Task 11.3 rev2)

% =========================================================================
%  SECTION 9 -- Sweep-mode-specific fields
% =========================================================================
switch sweep_mode

    case 'snr'
        % No extra fields.  SNR sweep vector passed externally.

    case 'bandwidth'
        P.snr_fixed = 10;   % [dB]; used when sweep_type = 'bandwidth'
        % K_s and alpha_k_vec overridden per sweep point inside the driver
        % using Delta_f_fixed = 25 MHz (Lesson L29).

    case 'convergence'
        % Fixed scene: theta=35 deg, r=3.0 m.
        % r=3.0 m is inside the locked operating box [1.063, 4.2525] m
        % (r_hi_fac=0.20, r_RD=21.2625 m).  Prior runs used r=5.0 m which
        % falls outside this box (5.0 > 4.2525 m) -- consistency fix (Sprint A).
        P.convergence_fixed_scene = true;
        P.conv_theta = 35 * pi/180;   % [rad]
        P.conv_r     = 3.0;           % [m]  -- was 5.0 (Sprint A fix)

end

% =========================================================================
%  CONVERGENCE-SCENE FRESNEL CHECK (P19: explicit post-override print + assert)
% =========================================================================
% This block runs AFTER the switch so P.conv_r reflects the authoritative value
% set above.  The main Fresnel check (Section 12) uses r_hi_fac and is NOT
% sufficient for the fixed convergence scene.  Print and assert here so that
% any future change to P.conv_r produces an immediate programmatic failure
% rather than a silent out-of-regime run.
if strcmp(sweep_mode, 'convergence')
    kappa_conv    = pi * P.d_ant^2 / P.lambda_c ...
                    * sin(P.conv_theta)^2 / P.conv_r;
    quad_phase_conv = kappa_conv * ((P.M - 1) / 2)^2;
    fprintf('  Conv scene    : theta=%.1f deg, r=%.2f m\n', ...
        P.conv_theta * 180/pi, P.conv_r);
    fprintf('  Fresnel (conv): kappa*m_bar_max^2 = %.4f rad', quad_phase_conv);
    assert(quad_phase_conv > 0.25, ...
        ['setup_production_P: conv_r too large -- insufficient Fresnel ' ...
         'curvature at fixed convergence scene (quad_phase=%.4f < 0.25).'], ...
        quad_phase_conv);
    fprintf('  [PASS]\n');
end

% =========================================================================
%  SECTION 10 -- CRB helper defaults (Tasks 12.4, 12.5, 12.6)
% =========================================================================
P.N_seed     = 50;    % random-W draws for CRB envelope
P.rng_seed_W = 0;     % base RNG seed for W draws inside wb_crb_compressed

% =========================================================================
%  SECTION 11 -- Diagnostic / CSV suppression  [v2 additions]
% =========================================================================
% wb_gen_write_csv = false: suppress per-realisation channel-gen CSV
%   (one file per MC realisation floods disk during N_MC=200 runs).
P.wb_gen_write_csv = false;

% bpd_write_csv = false: suppress per-realisation BPD diagnostic CSV
%   (bpd_baseline.m writes one row per realisation; with N_MC=200 and
%   6 parfor workers this produced ~200 console messages per sweep point
%   in production run 1, making the output unreadable).
%   To re-enable for single-realisation debugging:
%     P = setup_production_P(...); P.bpd_write_csv = true;
P.bpd_write_csv = false;

% use_riviello_snr_axis = false: when true, the SNR axis in the 'snr'
%   sweep is replaced by an empirical per-UT distribution drawn from the
%   3GPP UMi path-loss + shadow-fading model (Section VI-D robustness
%   experiment, Fig. 11).  Requires riviello_snr_axis.m on the MATLAB
%   path.  See riviello_snr_axis.m for the full parameter description.
%   Default: false (backward-compatible with all existing runs).
%   Session 4.9 (D1 execution).
P.use_riviello_snr_axis = false;

% =========================================================================
%  SECTION 12 -- Console summary and sanity assertions
% =========================================================================
fprintf('============================================================\n');
fprintf('  setup_production_P (v4):  Paper C Phase 3\n');
fprintf('  sweep_mode    = %s\n',  sweep_mode);
fprintf('  range_mode    = %s\n',  range_mode);
fprintf('------------------------------------------------------------\n');
fprintf('  fc            = %.1f GHz\n',   P.fc/1e9);
fprintf('  lambda_c      = %.6f m\n',     P.lambda_c);
fprintf('  d_ant         = %.6f m\n',     P.d_ant);
fprintf('  M             = %d\n',          P.M);
fprintf('  N_RF          = %d\n',          P.N_RF);
fprintf('  N             = %d\n',          P.N);
fprintf('  d             = %d\n',          P.d);
fprintf('  r_RD          = %.4f m\n',      P.r_RD);
fprintf('  r_lo_fac      = %.2f  ->  r_lo = %.4f m\n', ...
    P.r_lo_fac, P.r_lo_fac * P.r_RD);
fprintf('  r_hi_fac      = %.2f  ->  r_hi = %.4f m', ...
    P.r_hi_fac, P.r_hi_fac * P.r_RD);
if strcmp(range_mode, 'restricted')
    fprintf('  [RESTRICTED]\n');
else
    fprintf('  [FULL]\n');
end
fprintf('  B             = %.0f MHz\n',   P.B/1e6);
fprintf('  K_s           = %d\n',          P.K_s);
fprintf('  Delta_f       = %.1f MHz\n',   P.Delta_f/1e6);
fprintf('  alpha_k range = [%.6f, %.6f]\n', ...
    min(P.alpha_k_vec), max(P.alpha_k_vec));
fprintf('  Q_theta       = %d\n',          P.Q_theta);
fprintf('  G_theta       = %d\n',          P.G_theta);
fprintf('  G_r           = %d   (v1:64->v2:128)\n', P.G_r);
fprintf('  max_iter      = %d\n',          P.max_iter);
fprintf('  eps_reg       = %.4g  [v3: SNR-adaptive loading]\n', P.eps_reg);
  fprintf('  bpd_write_csv = false  [v2: console-spam fix]\n');
fprintf('------------------------------------------------------------\n');

% --- Alias assertions ---
assert(P.lambda_c == P.lambda, ...
    'setup_production_P: lambda_c/lambda alias mismatch.');
assert(P.c0 == P.c, ...
    'setup_production_P: c0/c alias mismatch.');
assert(numel(P.alpha_k_vec) == P.K_s, ...
    'setup_production_P: alpha_k_vec length != K_s.');
assert(P.K == P.K_s, ...
    'setup_production_P: K/K_s alias mismatch.');
assert(P.d_max == P.d, ...
    'setup_production_P: d_max/d mismatch (update externally if d changes).');
fprintf('  Alias checks  : PASS (lambda, c, K, d_max)\n');

% --- r_RD consistency ---
r_RD_check = 2 * ((P.M-1)*P.d_ant)^2 / P.lambda_c;
assert(abs(P.r_RD - r_RD_check) < 1e-9, ...
    'setup_production_P: r_RD computation inconsistency.');
fprintf('  r_RD sanity   : PASS (%.4f m)\n', P.r_RD);

% --- alpha_k range ---
assert(all(P.alpha_k_vec > 0.98) && all(P.alpha_k_vec < 1.02), ...
    'setup_production_P: alpha_k_vec out of expected range [0.98, 1.02].');
fprintf('  alpha_k sanity: PASS (range [%.5f, %.5f])\n', ...
    min(P.alpha_k_vec), max(P.alpha_k_vec));

% --- Fresnel curvature check at r_hi ---
kappa_r_hi = pi * P.d_ant^2 / P.lambda_c ...
             * sin(40*pi/180)^2 / (P.r_hi_fac * P.r_RD);
quad_phase = kappa_r_hi * ((P.M-1)/2)^2;
fprintf('  Fresnel check : kappa*m_bar_max^2 at r_hi (theta=40 deg) = %.4f rad', ...
    quad_phase);
if quad_phase >= 0.25
    fprintf('  [PASS: sufficient curvature for BPD]\n');
else
    fprintf('  [WARN: near-Rayleigh -- BPD warm-start unreliable]\n');
end

fprintf('============================================================\n');

end  % setup_production_P
