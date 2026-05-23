function a = wb_nf_fresnel_steer(theta, u, alpha_k_val, P)
%WB_NF_FRESNEL_STEER  Wideband Fresnel near-field steering vector.
%
%  a = wb_nf_fresnel_steer(theta, u, alpha_k_val, P)
%
%  Shared Paper C Phase 2 utility.  Extends the Paper B narrowband
%  Fresnel formula to subcarrier k by substituting an effective
%  wavelength lambda_eff = lambda_c / alpha_k_val.  This scales both
%  omega and kappa by alpha_k_val, implementing the locked Paper C
%  wideband steering vector (GLOBECOM eq. a_k):
%
%    [a_{l,k}]_m = exp( j*alpha_k*(m_bar*omega_l + m_bar^2*kappa_l) )
%
%  where
%    omega_l  = (2*pi*d_ant/lambda_c)*cos(theta)     [linear phase]
%    kappa_l  = (pi*d_ant^2/lambda_c)*sin(theta)^2*u [quadratic phase]
%    alpha_k  = f_k/fc = lambda_c/lambda_eff
%
%  At alpha_k_val = 1 (carrier frequency), this reduces exactly to the
%  Paper B narrowband formula in nf_fresnel_steer.m.
%
%  SIGN CONVENTION
%  ---------------
%  omega uses POSITIVE cosine: omega = +(2*pi*d_ant/lambda_eff)*cos(theta).
%  This matches Paper B's nf_fresnel_steer convention.  All Paper C
%  estimators and the channel generator wb_channel_gen_ofdm_nf (Task 11.2)
%  MUST use the same sign convention.
%
%  OUTPUT NORMALISATION
%  --------------------
%  a is normalised to unit L2 norm: norm(a) = 1.
%  To recover the unnormalised atom (||a||^2 = M), multiply by sqrt(M).
%
%  INPUTS
%  ------
%  theta      : scalar angle of arrival [rad],  0 < theta < pi/2
%  u          : scalar inverse range 1/r [1/m], u > 0
%  alpha_k_val: scalar frequency scaling f_k/fc, typically in [0.99, 1.01]
%               for 5G NR FR2 with K=2048, Delta_f=120 kHz at 28 GHz
%  P          : parameter struct; only three fields are accessed:
%                 .M        -- ULA element count
%                 .lambda_c -- carrier wavelength [m]  (= c0/fc)
%                 .d_ant    -- element spacing [m]     (= lambda_c/2)
%
%  OUTPUT
%  ------
%  a : M x 1 complex unit-norm steering vector
%
%  USAGE EXAMPLES
%  --------------
%  % Carrier-frequency atom (alpha_k = 1):
%    a_c = wb_nf_fresnel_steer(30*pi/180, 1/5.0, 1.0, P);
%
%  % Atom at subcarrier k with scaling alpha_k(k):
%    a_k = wb_nf_fresnel_steer(theta, 1/r, P.alpha_k_vec(k), P);
%
%  % Unnormalised atom for channel generation (norm = sqrt(M)):
%    a_unnorm = wb_nf_fresnel_steer(theta, 1/r, alpha_k, P) * sqrt(P.M);
%
%  DEPENDENCIES
%  ------------
%  None.  Pure MATLAB, no toolboxes required.
%
%  CALLED BY
%  ---------
%  bpd_baseline.m         -- baseline B4 dictionary and residual update
%  wb_channel_gen_ofdm_nf.m (Task 11.2) -- channel generator
%  (future) wb_cl_kl.m    -- proposed estimator
%
%  Author : R. V. Senyuva (Maltepe University)
%  Date   : May 2026

M          = P.M;
lambda_c   = P.lambda_c;
d_ant      = P.d_ant;

lambda_eff = lambda_c / alpha_k_val;          % effective wavelength at f_k

m_bar = ((0:M-1).' - (M-1)/2);               % M x 1, centred element indices

omega = (2*pi * d_ant / lambda_eff) * cos(theta);       % linear phase slope
c_i   = (pi  * d_ant^2 / lambda_eff) * sin(theta)^2;   % curvature scale
kappa = c_i * u;                                        % quadratic coefficient

a = exp(1j * omega * m_bar - 1j * kappa * m_bar.^2);   % M x 1
a = a / norm(a);                                        % unit L2 norm

end
