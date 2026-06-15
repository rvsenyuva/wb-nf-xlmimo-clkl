# wb-nf-xlmimo-clkl

**Wideband Near-Field Channel Estimation under Hybrid Compression:
Cross-Subcarrier KL Covariance Fitting with OFDM Fresnel Model**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20356437.svg)](https://doi.org/10.5281/zenodo.20356437)

MATLAB simulation code and production CSV results for Paper C,
submitted to IEEE Open Journal of the Communications Society (OJ-COMS), 2026.

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{senyuva2026wbnf,
  author  = {{\c{S}}enyuva, R. Volkan},
  title   = {Wideband Near-Field Channel Estimation under Hybrid Compression:
             Cross-Subcarrier {KL} Covariance Fitting with {OFDM} {F}resnel Model},
  journal = {IEEE Open Journal of the Communications Society},
  year    = {2026},
  note    = {Submitted. arXiv preprint arXiv:[to be assigned]}
}
```

Also cite the companion GLOBECOM 2026 paper for the CRB derivation:

```bibtex
@inproceedings{senyuva2026globecom,
  author    = {{\c{S}}enyuva, R. Volkan},
  title     = {Wideband Compressed-Domain {C}ram\'{e}r--{R}ao Bounds for
               Near-Field {XL-MIMO}: Data and Geometric Diversity Decomposition},
  booktitle = {Proc. IEEE Global Commun. Conf. (GLOBECOM)},
  year      = {2026},
  note      = {arXiv:2604.08531}
}
```

---

## Prerequisites

- MATLAB R2023b or later (tested on R2025b)
- Signal Processing Toolbox
- Parallel Computing Toolbox (required for `run_monte_carlo_paperC.m`; 6 workers recommended)

No additional toolboxes are required. All near-field steering vector computations
use closed-form Fresnel approximations implemented in `src/wb_nf_fresnel_steer.m`.

---

## Quick start

```matlab
% 1. Clone the repo and navigate to it
% 2. Add src/ to the MATLAB path (run from repo root)
addpath('src')

% 3. Verify all modules pass their unit tests
cd tests
wb_clkl_toy_test_v2           % Expected: T1-T5 all PASS
wb_clkl_driver_toy_test_v2    % Expected: T1-T6 all PASS
wb_crb_compressed_toy_test    % Expected: T1-T4 all PASS
wb_psomp_toy_test              % Expected: T1-T4 all PASS
wb_dl_omp_toy_test             % Expected: T1-T4 all PASS
wb_channel_gen_toy_test        % Expected: smoke test passes
bpd_toy_test                   % Expected: smoke test passes
mc_toy_test                    % Expected: smoke test passes (5 MC trials)
cd ..

% 4. Reproduce all paper figures from the committed CSVs
cd plotting
fig8a_rmse_r_snr       % Fig 1: Range RMSE vs SNR, all methods
fig8b_nmse_r_bw        % Fig 2: Range NMSE vs bandwidth, all methods
fig9a_rmse_theta_snr   % Fig 3: Angle RMSE vs SNR, all methods
fig9b_convpct_snr_b4   % Fig 4: WB-CL-KL convergence rate vs SNR
fig10_kl_convergence   % Fig 5: KL objective vs iteration
fig6_crb_d2_snr        % Fig 6: CRB compression gap vs SNR
fig7_crb_geodiv_bw     % Fig 7: Geometric diversity CRB vs bandwidth
fig11_riviello_snr_axis % Fig 8: Robustness under 3GPP UMi SNR distribution
cd ..
```

Each plotting script reads directly from the committed CSVs in `results/` and
exports a PDF to the `plotting/` folder. No re-running the Monte Carlo is needed
to reproduce the figures.

To run the full Monte Carlo from scratch (SNR sweep, ~8 min on 6 parfor workers):

```matlab
addpath('src')
% Run from repo root (not from inside src/)
run_monte_carlo_paperC
```

The Monte Carlo driver reads parameters from `setup_production_P_v4.m`.

---

## Repository structure

```
wb-nf-xlmimo-clkl/
|
+-- README.md
+-- LICENSE
+-- .gitignore
|
+-- src/                              12 canonical MATLAB source files
|   +-- wb_clkl_estimator.m           Proposed WB-CL-KL estimator (v2; SNR-adaptive loading)
|   +-- wb_clkl_driver.m              Driver with C3-patch (3-way Phase D argmin)
|   +-- setup_production_P_v4.m       Production parameter struct (r_hi_fac=0.20)
|   +-- run_monte_carlo_paperC.m      Master Monte Carlo driver (parfor, N_MC=600)
|   +-- wb_crb_compressed.m           Wideband compressed-domain CRB
|   +-- bpd_baseline.m                BPD baseline estimator (B1: WB-BPD full-array)
|   +-- wb_psomp.m                    WB-P-SOMP baseline (B2: compressed)
|   +-- wb_dl_omp.m                   WB-DL-OMP baseline (B5: full-array)
|   +-- wb_nf_fresnel_steer.m         Shared Fresnel/USW steering vector utility
|   +-- wb_channel_gen_ofdm_nf.m      Wideband near-field OFDM channel generator
|   +-- nf_usw_steer.m                USW steering vector (dependency of wb_channel_gen_ofdm_nf)
|   +-- riviello_snr_axis.m           3GPP UMi per-UT SNR distribution (Fig 8 robustness)
|
+-- tests/                            8 unit test / toy test scripts
|   +-- wb_clkl_toy_test_v2.m
|   +-- wb_clkl_driver_toy_test_v2.m
|   +-- wb_crb_compressed_toy_test.m
|   +-- wb_psomp_toy_test.m
|   +-- wb_dl_omp_toy_test.m
|   +-- wb_channel_gen_toy_test.m
|   +-- bpd_toy_test.m
|   +-- mc_toy_test.m
|
+-- diagnostics/                      1 forensic diagnostic script
|   +-- wb_clkl_forensic_single.m     Single-realisation KL trajectory diagnostic
|
+-- plotting/                         8 publication figure scripts
|   +-- fig8a_rmse_r_snr.m            Fig 1: RMSE_r vs SNR, all methods
|   +-- fig8b_nmse_r_bw.m             Fig 2: NMSE_r vs bandwidth, all methods
|   +-- fig9a_rmse_theta_snr.m        Fig 3: RMSE_theta vs SNR, all methods
|   +-- fig9b_convpct_snr_b4.m        Fig 4: WB-CL-KL convergence rate vs SNR
|   +-- fig10_kl_convergence.m        Fig 5: KL objective delta_L vs iteration
|   +-- fig6_crb_d2_snr.m             Fig 6: sqrt-CRB_r vs SNR, N_RF in {8,16,32,64}, d=2
|   +-- fig7_crb_geodiv_bw.m          Fig 7: sqrt-CRB_r vs bandwidth, M in {32,64,128}
|   +-- fig11_riviello_snr_axis.m     Fig 8: Robustness, 3GPP UMi SNR distribution
|
+-- results/                          Production CSV files
    +-- regime_probe_rhi020/
    |   +-- mc_snr_sweep_20260528_204613.csv             (65 rows x 24 cols; Figs 1,3,4)
    |
    +-- paperC_phase35_fig8b_bw_r020_N600/
    |   +-- mc_bandwidth_sweep_20260529_210446.csv       (40 rows x 24 cols; Fig 2)
    |
    +-- paperC_sprint_a_fig10_convergence_r3m/
    |   +-- mc_convergence_sweep_20260531_215757.csv     (20 rows x 24 cols; Fig 5)
    |   +-- mc_convergence_Lhist_SNR+0dB_20260531_215757.csv   (9431 rows x 4 cols)
    |   +-- mc_convergence_Lhist_SNR+5dB_20260531_215757.csv   (9642 rows x 4 cols)
    |   +-- mc_convergence_Lhist_SNR+10dB_20260531_215757.csv  (9621 rows x 4 cols)
    |   +-- mc_convergence_Lhist_SNR+15dB_20260531_215757.csv  (9663 rows x 4 cols)
    |
    +-- paperC_phase35_fig6_crb_d2_v2/
    |   +-- wb_crb_d2_sweep_20260526_225015_v3.csv       (52 rows x 10 cols; Fig 6)
    |
    +-- paperC_phase3_fig7_crb_geodiv_v1/
    |   +-- wb_crb_geodiv_sweep_20260520_184631.csv      (18 rows x 7 cols; Fig 7)
    |
    +-- riviello_snr/
        +-- mc_snr_sweep_20260610_203626.csv             (50 rows x 24 cols; Fig 8 top)
        +-- riviello_snr_distribution_2000.csv           (2000 rows x 1 col; Fig 8 bottom)
```

---

## Reproducing the paper figures

All figures are generated directly from the CSVs in `results/`. The table below
maps each paper figure to its plotting script and source CSV.

| Figure | Description | Script | Source CSV |
|--------|-------------|--------|------------|
| Fig 1 | Range RMSE vs SNR, all methods, SNR in [-5, 17.5] dB | `plotting/fig8a_rmse_r_snr.m` | `results/regime_probe_rhi020/mc_snr_sweep_20260528_204613.csv` |
| Fig 2 | Range NMSE vs bandwidth, all methods, B in {100,...,800} MHz | `plotting/fig8b_nmse_r_bw.m` | `results/paperC_phase35_fig8b_bw_r020_N600/mc_bandwidth_sweep_20260529_210446.csv` |
| Fig 3 | Angle RMSE vs SNR, all methods | `plotting/fig9a_rmse_theta_snr.m` | same as Fig 1 |
| Fig 4 | WB-CL-KL convergence rate vs SNR (strict tolerance) | `plotting/fig9b_convpct_snr_b4.m` | same as Fig 1 |
| Fig 5 | Median KL objective delta_L vs iteration, SNR in {0,5,10,15} dB | `plotting/fig10_kl_convergence.m` | `results/paperC_sprint_a_fig10_convergence_r3m/mc_convergence_Lhist_SNR+*dB_20260531_215757.csv` |
| Fig 6 | Compression gap: sqrt-CRB_r vs SNR, N_RF in {8,16,32,64}, d=2 | `plotting/fig6_crb_d2_snr.m` | `results/paperC_phase35_fig6_crb_d2_v2/wb_crb_d2_sweep_20260526_225015_v3.csv` |
| Fig 7 | Geometric diversity: sqrt-CRB_r vs bandwidth, M in {32,64,128} | `plotting/fig7_crb_geodiv_bw.m` | `results/paperC_phase3_fig7_crb_geodiv_v1/wb_crb_geodiv_sweep_20260520_184631.csv` |
| Fig 8 | Robustness under 3GPP UMi SNR distribution (two-panel) | `plotting/fig11_riviello_snr_axis.m` | `results/riviello_snr/mc_snr_sweep_20260610_203626.csv` + `riviello_snr_distribution_2000.csv` |

Each script prints anchor-verification lines (OK / ANCHOR MISMATCH) to the
command window before exporting the PDF.

---

## Key verified numerical claims

All values are CSV-sourced and anchor-verified by the plotting scripts.

**Estimator performance (r_hi_fac=0.20, strong near-field regime, N_MC=600):**
- B4 (WB-CL-KL) RMSE_r at SNR=10 dB: **0.01982 m** vs CRB 0.01991 m (ratio 0.996; CRB-efficient)
- B4 NMSE_r at B=400 MHz, SNR=10 dB: **-43.16 dB**
- B4 convergence rate (strict tolerance) at SNR=10 dB: **73.0%** (fail_rate 0.0%)
- B4 convergence rate at SNR=0 dB: **100%** (all trials converge at low SNR)
- B4/CRB at median deployment SNR (9.6 dB, 3GPP UMi): **0.959** (Fig 8 robustness gate)

**CRB decomposition (B=400 MHz, M=64, N_RF=8, d=1):**
- Total wideband CRB gain: **+27.793 dB**
  (data diversity: +27.093 dB; geometric diversity: +0.701 dB)
- Compression gap N_RF=8 vs N_RF=32 at SNR=10 dB (d=2): **8.12 dB**
- Compression gap N_RF=8 vs N_RF=64 at SNR=10 dB (d=2): **11.40 dB**

**CRB absolute values at B=400 MHz, SNR=10 dB, r=5 m:**
- sqrt-CRB_r (narrowband, B→0): **11.948 mm**
- sqrt-CRB_r (wideband, B=400 MHz): **487.12 µm**

**Geometric diversity (d=1, SNR=10 dB, theta=35 deg, r=5 m):**
- CRB gain from B=50 to B=400 MHz (M=64): **9.125 dB**
- M-doubling gap M=64 to M=128 at B=400 MHz: **12.629 dB**

---

## Method labels

| Label | Method | Architecture | Original source |
|-------|--------|--------------|-----------------|
| B1 | WB-BPD | Full-array (M RF chains) | Cui & Dai, Sci. China Inf. Sci. 2023 |
| B2 | WB-P-SOMP | Hybrid compressed (N_RF RF chains) | Cui & Dai, IEEE Trans. Commun. 2022 |
| B4 | WB-CL-KL (proposed) | Hybrid compressed (N_RF RF chains) | This work |
| B5 | WB-DL-OMP | Full-array (M RF chains) | Zhang et al., IEEE Trans. Commun. 2024 |

B3 (WB-BF-SOMP) is defined in the paper but not implemented; its entries in
the CSVs are NaN placeholders.

---

## System parameters (production run)

| Parameter | Value | Description |
|-----------|-------|-------------|
| f_c | 28 GHz | Carrier frequency |
| M | 64 | Number of ULA antenna elements |
| N_RF | 8 | Number of RF chains (compressed) |
| N | 64 | Number of pilot snapshots |
| d | 1 | Path count (single-path default) |
| r_hi_fac | 0.20 | Range window upper factor (r_hi = 0.20 × r_RD) |
| r_RD | 21.26 m | Rayleigh distance (M=64, f_c=28 GHz) |
| r range | [0.2126, 4.2525] m | Scene range window (strong near-field) |
| B | 400 MHz | Bandwidth (SNR sweep default) |
| K_s | 16 | Number of OFDM subcarriers (at B=400 MHz) |
| SNR range | -5 to 17.5 dB | 10-point window (2.5 dB step) |
| N_MC | 600 | Monte Carlo trials per point (Figs 1–5) |
| N_MC | 200 | Monte Carlo trials per point (Fig 8 robustness) |

---

## Operating regime note

All estimator results (Figs 1–5) are evaluated in the **strong near-field regime**
defined by r_hi_fac=0.20, placing targets well within the Rayleigh distance
(r ≤ 0.20 × r_RD = 4.25 m). This is the physically valid and scientifically
interesting operating regime for the WB-CL-KL estimator, where the Fresnel
quadratic phase is large enough for reliable range discrimination. CRB figures
(Figs 6–7) are evaluated at representative scene points. Fig 8 extends the
SNR sweep to [-20, +35] dB and overlays the 3GPP UMi empirical SNR distribution
(N_UT=2000, median 9.6 dB) for robustness verification.

---

## Related repositories

- **Paper B** (narrowband near-field, KL covariance fitting):
  https://github.com/rvsenyuva/nearfield-clkl

- **GLOBECOM 2026** (wideband compressed-domain CRB derivation):
  https://github.com/rvsenyuva/wb-nf-crb-globecom26
  DOI: 10.5281/zenodo.19487208

- **Paper A** (beamspace ESPRIT, mmWave sensor arrays):
  https://github.com/rvsenyuva/CovGuided-ESPRIT
  DOI: 10.5281/zenodo.19553631

---

## License

MIT License. See [LICENSE](LICENSE).

This code is for academic research only. If you use it, please cite the paper above.
