function [snr_per_ut, info] = riviello_snr_axis(P, N_UT)
%RIVIELLO_SNR_AXIS  Realistic 3GPP UMi per-UT SNR distribution for Paper C.
%
%  [snr_per_ut, info] = riviello_snr_axis(P, N_UT)
%
%  Generates an empirical per-user-terminal (UT) SNR distribution by
%  sampling UT positions in a hexagonal UMi cell and applying the 3GPP
%  TR 38.901 path-loss / shadow-fading model at 28 GHz.  Used as an
%  alternative SNR axis for the Section VI-D robustness experiment
%  (Fig. 11), replacing the controlled uniform SNR grid with a
%  deployment-realistic distribution.
%
%  Implements:
%    - UMi LOS probability  : Riviello 2022 Eq. (31) / TR 38.901 Table 7.4.2-1
%    - UMi path loss (LOS)  : Riviello 2022 Eqs. (32)-(33) / TR 38.901 Table 7.4.1-1
%    - UMi path loss (NLOS) : TR 38.901 Table 7.4.1-1 (Riviello omits explicit form)
%    - O2I penetration loss : TR 38.901 Table 7.4.3-1 high-loss concrete model
%    - Shadow fading        : Riviello 2022 Eq. (36) / TR 38.901 Table 7.5-6
%
%  INPUTS
%    P     : Paper C parameter struct (needs P.fc [Hz], P.M, P.B [Hz]).
%            The SNR budget uses P.M for array gain and P.B for
%            thermal noise floor.  Other fields unused here.
%    N_UT  : number of UTs to draw (scalar positive integer).
%            Set equal to N_MC for the Fig. 11 run.
%
%  OUTPUTS
%    snr_per_ut : N_UT x 1 vector of per-UT SNR values [dB]
%    info       : struct with fields
%                   median_snr_dB  -- median of snr_per_ut [dB]
%                   mean_snr_dB    -- mean of snr_per_ut [dB]
%                   p10_snr_dB     -- 10th percentile [dB]
%                   p90_snr_dB     -- 90th percentile [dB]
%                   pct_indoor     -- fraction of indoor UTs [0,100]
%                   pct_los        -- fraction of LOS UTs [0,100]
%                   N_UT           -- N_UT (echo for logging)
%
%  DESIGN NOTES
%    - Fully vectorized; no per-UT for-loops.
%    - Self-contained; no external toolbox or QuaDRiGa dependency.
%    - Does NOT draw the full 12-step TR 38.901 cluster model; only the
%      path-loss and single-LSP (sigma_SF) shadow-fading subset is used.
%    - Path-loss recipe is attribution-only from Riviello 2022 (MIT licence).
%      The function is an independent implementation of TR 38.901 Sec. 7.4.
%
%  USAGE (Fig. 11 run)
%    P = setup_production_P_v4('snr', 'restricted');
%    P.use_riviello_snr_axis = true;
%    P.N_MC = 200;
%    [snr_per_ut, snr_info] = riviello_snr_axis(P, P.N_MC);
%
%  Path-loss / shadowing recipe adapted from:
%    D. G. Riviello, F. Di Stasio, R. Tuninato, "Performance Analysis of
%    Multi-User MIMO Schemes under Realistic 3GPP 3-D Channel Model for
%    5G mmWave Cellular Networks," Electronics, 11(3):330, 2022.
%    MIT-licensed at gitlab.com/daniel.riviello/3gpp-channel-model-tr-38901
%
%  Author   : R. V. Senyuva (Maltepe University)
%  Date     : June 2026
%  Ref      : Paper C Phase 4, Session 4.9 (D1 robustness experiment).
%  Called by: run_monte_carlo_paperC.m (when P.use_riviello_snr_axis = true)
%             fig11_riviello_snr_axis.m (plotting script)

% =========================================================================
%  INPUT VALIDATION
% =========================================================================
narginchk(2, 2);
assert(isstruct(P),        'riviello_snr_axis: P must be a struct.');
assert(isfield(P, 'fc'),   'riviello_snr_axis: P must have field fc.');
assert(isfield(P, 'M'),    'riviello_snr_axis: P must have field M.');
assert(isfield(P, 'B'),    'riviello_snr_axis: P must have field B.');
assert(isscalar(N_UT) && N_UT >= 1, ...
    'riviello_snr_axis: N_UT must be a positive scalar integer.');
N_UT = round(N_UT);

% =========================================================================
%  SECTION 1 -- Physical constants and carrier
% =========================================================================
f_c_GHz = P.fc / 1e9;   % carrier frequency [GHz]   (28 GHz)
B_Hz    = P.B;           % bandwidth [Hz]             (400e6)
M_ant   = P.M;           % number of ULA elements     (64)

% =========================================================================
%  SECTION 2 -- Link budget constants
% =========================================================================
% TR 38.901 UMi cell geometry
h_BS  = 10.0;   % BS height [m]
h_UT_outdoor = 1.5;   % outdoor UT height [m] (TR 38.901 default)

% Transmit power
P_TX_dBm = 30.0;   % [dBm]  typical 5G UMi small cell (1 W)

% Array gain: M isotropic elements => G_array = 10*log10(M) dBi
G_array_dBi = 10 * log10(M_ant);   % [dBi]

% Thermal noise floor: N_thermal = -174 + 10*log10(B) [dBm]
N_thermal_dBm = -174 + 10 * log10(B_Hz);   % [dBm]

% Noise figure (typical mmWave LNA)
NF_dB = 7.0;   % [dB]

% Noise power
N_power_dBm = N_thermal_dBm + NF_dB;   % [dBm]

% O2I penetration: concrete high-loss model (TR 38.901 Table 7.4.3-1)
% Simplified: 20 dB fixed (conservative; dominates over glass 10 dB)
PL_O2I_dB = 20.0;   % [dB]

% Indoor shadow fading sigma (TR 38.901 Table 7.5-6 UMi NLOS):
sigma_SF_nlos = 7.82;   % [dB]
sigma_SF_los  = 4.00;   % [dB]
% (Outdoor NLOS UTs: use NLOS sigma; LOS UTs: use LOS sigma)

% LOS breakpoint distance (TR 38.901 Table 7.4.1-1):
% d_BP = 4 * h_BS_eff * h_UT_eff * f_c / c
% Effective heights: h_BS_eff = h_BS - 1.0; h_UT_eff = h_UT - 1.0
h_BS_eff = h_BS - 1.0;   % [m]
h_UT_eff = h_UT_outdoor - 1.0;   % [m]  (use outdoor height for d_BP)
h_UT_eff = max(h_UT_eff, 0.01);  % guard against non-positive
c0       = 3e8;   % speed of light [m/s]
d_BP     = 4 * h_BS_eff * h_UT_eff * P.fc / c0;   % [m]

% =========================================================================
%  SECTION 3 -- UT placement (vectorized over N_UT)
%  TR 38.901: d_2D uniform in [10 m, 100 m] (UMi ISD=200 m hex cell)
%  Azimuth uniform in [0, 2*pi)
% =========================================================================
d_2D_min = 10.0;    % [m]
d_2D_max = 100.0;   % [m]

d_2D    = d_2D_min + (d_2D_max - d_2D_min) * rand(N_UT, 1);   % N_UT x 1
% azimuth not needed for path loss (isotropic scenario)

% =========================================================================
%  SECTION 4 -- Indoor / outdoor assignment
%  TR 38.901 Table 7.4.3-1 UMi: 80% indoor, 20% outdoor
% =========================================================================
p_indoor   = 0.8;
is_indoor  = rand(N_UT, 1) < p_indoor;   % logical N_UT x 1

% Indoor UT height: ground-floor assumption (no floor-level elevation for
% conservative worst-case; floor number not needed when using fixed O2I)
h_UT = h_UT_outdoor * ones(N_UT, 1);   % all UTs: 1.5 m above grade

% 3D distance
d_3D = sqrt(d_2D.^2 + (h_BS - h_UT).^2);   % N_UT x 1

% =========================================================================
%  SECTION 5 -- LOS probability (Riviello Eq. 31 / TR 38.901 Table 7.4.2-1)
% =========================================================================
% P_LOS = 1                                          if d_2D <= 18 m
% P_LOS = 18/d_2D + exp(-d_2D/36) * (1 - 18/d_2D)  if d_2D >  18 m

p_los = (18 ./ d_2D) + exp(-d_2D / 36) .* (1 - 18 ./ d_2D);
p_los(d_2D <= 18) = 1.0;
p_los = min(max(p_los, 0), 1);   % clip to [0,1] for safety

is_los = rand(N_UT, 1) < p_los;  % logical N_UT x 1 (Bernoulli draw)

% =========================================================================
%  SECTION 6 -- LOS path loss (Riviello Eqs. 32-33 / TR 38.901 Table 7.4.1-1)
%
%  PL_LOS_1 = 32.4 + 21*log10(d_3D) + 20*log10(f_c_GHz)   d_2D <= d_BP
%  PL_LOS_2 = 32.4 + 40*log10(d_3D) + 20*log10(f_c_GHz)
%             - 9.5*log10(d_BP^2 + (h_BS - h_UT)^2)        d_2D >  d_BP
% =========================================================================
log10_fc = log10(f_c_GHz);
log10_d3 = log10(d_3D);

PL_LOS_1 = 32.4 + 21 * log10_d3 + 20 * log10_fc;
d_BP_sq  = d_BP^2 + (h_BS - h_UT_outdoor)^2;   % scalar breakpoint factor
PL_LOS_2 = 32.4 + 40 * log10_d3 + 20 * log10_fc - 9.5 * log10(d_BP_sq);

PL_LOS   = PL_LOS_1;
PL_LOS(d_2D > d_BP) = PL_LOS_2(d_2D > d_BP);

% =========================================================================
%  SECTION 7 -- NLOS path loss (TR 38.901 Table 7.4.1-1 UMi NLOS)
%
%  PL_NLOS = max(PL_LOS, PL_NLOS_3GPP) per standard recommendation.
%  PL_NLOS_3GPP = 35.3*log10(d_3D) + 22.4 + 21.3*log10(f_c_GHz) - 0.3*(h_UT-1.5)
%  Valid for 10 m <= d_2D <= 5 km (TR 38.901 Table 7.4.1-1 row 5)
% =========================================================================
PL_NLOS_3GPP = 35.3 * log10_d3 + 22.4 ...
             + 21.3 * log10_fc ...
             - 0.3 * (h_UT - 1.5);

PL_NLOS = max(PL_LOS, PL_NLOS_3GPP);

% =========================================================================
%  SECTION 8 -- Assign path loss per UT based on LOS/NLOS state
% =========================================================================
PL_basic = is_los .* PL_LOS + (~is_los) .* PL_NLOS;   % N_UT x 1 [dB]

% =========================================================================
%  SECTION 9 -- O2I penetration loss for indoor UTs
% =========================================================================
PL_total = PL_basic + is_indoor * PL_O2I_dB;   % N_UT x 1 [dB]

% =========================================================================
%  SECTION 10 -- Shadow fading (Riviello Eq. 36 / TR 38.901 Table 7.5-6)
%  Independent lognormal per UT.
%  LOS sigma = 4.0 dB; NLOS sigma = 7.82 dB.
%  Indoor UTs: NLOS sigma (conservative; most indoor UTs are NLOS).
% =========================================================================
sigma_SF = sigma_SF_los * is_los + sigma_SF_nlos * (~is_los);   % N_UT x 1
SF_dB    = sigma_SF .* randn(N_UT, 1);   % zero-mean lognormal shadow [dB]

% =========================================================================
%  SECTION 11 -- Per-UT SNR
%  SNR_dB = P_TX_dBm - PL_total - SF_dB + G_array_dBi - N_power_dBm
% =========================================================================
snr_per_ut = P_TX_dBm - PL_total - SF_dB + G_array_dBi - N_power_dBm;
% snr_per_ut is N_UT x 1, in dB

% =========================================================================
%  SECTION 12 -- Info struct for caption reporting
% =========================================================================
info.median_snr_dB = median(snr_per_ut);
info.mean_snr_dB   = mean(snr_per_ut);
info.p10_snr_dB    = prctile(snr_per_ut, 10);
info.p90_snr_dB    = prctile(snr_per_ut, 90);
info.pct_indoor    = 100 * mean(is_indoor);
info.pct_los       = 100 * mean(is_los);
info.N_UT          = N_UT;

% =========================================================================
%  SECTION 13 -- Console summary
% =========================================================================
fprintf('[riviello_snr_axis] N_UT=%d  fc=%.1f GHz  B=%.0f MHz  M=%d\n', ...
    N_UT, f_c_GHz, B_Hz/1e6, M_ant);
fprintf('  P_TX=%.1f dBm  G_array=%.2f dBi  N_thermal=%.2f dBm  NF=%.1f dB\n', ...
    P_TX_dBm, G_array_dBi, N_thermal_dBm, NF_dB);
fprintf('  O2I=%.1f dB  d_BP=%.1f m\n', PL_O2I_dB, d_BP);
fprintf('  SNR [dB]: median=%.2f  mean=%.2f  p10=%.2f  p90=%.2f\n', ...
    info.median_snr_dB, info.mean_snr_dB, info.p10_snr_dB, info.p90_snr_dB);
fprintf('  Indoor: %.1f%%   LOS: %.1f%%\n', info.pct_indoor, info.pct_los);

end  % riviello_snr_axis
