# wb-nf-xlmimo-clkl

**Wideband Near-Field Channel Estimation under Hybrid Compression:
Cross-Subcarrier KL Covariance Fitting with OFDM Fresnel Model**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20356437.svg)](https://doi.org/10.5281/zenodo.20356437)

MATLAB simulation code and production CSV results for Paper C.
Companion to arXiv preprint [to be assigned].

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{senyuva2026wbnf,
  author  = {Senyuva, R. Volkan},
  title   = {Wideband Near-Field Channel Estimation under Hybrid Compression:
             Cross-Subcarrier {KL} Covariance Fitting with {OFDM} {F}resnel Model},
  journal = {IEEE Open Journal of the Communications Society},
  year    = {2026},
  note    = {Submitted. Preprint: arXiv:[to be assigned]}
}
```

Also cite the companion GLOBECOM 2026 paper for the CRB derivation:

```bibtex
@inproceedings{senyuva2026globecom,
  author    = {Senyuva, R. Volkan},
  title     = {Wideband Compressed-Domain {C}ram\'{e}r--{R}ao Bounds for
               Near-Field {XL-MIMO}: Data and Geometric Diversity Decomposition},
  booktitle = {Proc. IEEE GLOBECOM},
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
% 2. Add src/ to the MATLAB path
addpath('src')

% 3. Verify all modules pass their unit tests
cd tests
wb_clkl_toy_test_v2           % Expected: T1-T5 all PASS
wb_clkl_driver_toy_test_v2    % Expected: T1-T6 all PASS
wb_crb_compressed_toy_test    % Expected: T1-T4 all PASS
wb_psomp_toy_test              % Expected: T1-T4 all PASS
wb_dl_omp_toy_test             % Expected: T1-T4 all PASS
mc_toy_test                    % Expected: smoke test passes (5 MC trials)
cd ..

% 4. Run the full Monte Carlo (SNR sweep, ~5 min on 6 parfor workers)
cd src
run_monte_carlo_paperC
```

The Monte Carlo driver reads parameters from `setup_production_P_v4.m`.
Default sweep: SNR in [-5, 0, 5, 10, 15] dB, N_MC=200, B=400 MHz, M=64, N_RF=8.

---

## Repository structure

```
wb-nf-xlmimo-clkl/
|
+-- README.md
+-- LICENSE
+-- .gitignore
|
+-- src/                          11 canonical MATLAB source files
|   +-- wb_clkl_estimator.m       Proposed WB-CL-KL estimator (v2; SNR-adaptive loading)
|   +-- wb_clkl_driver.m          Driver with C3-patch (3-way Phase D argmin)
|   +-- setup_production_P_v4.m   Production parameter struct
|   +-- run_monte_carlo_paperC.m  Master Monte Carlo driver (745 lines; parfor)
|   +-- wb_crb_compressed.m       Wideband compressed-domain CRB
|   +-- bpd_baseline.m            BPD baseline estimator (B1: WB-BPD full-array)
|   +-- wb_psomp.m                WB-P-SOMP baseline (B2: compressed)
|   +-- wb_dl_omp.m               WB-DL-OMP baseline (B5: full-array)
|   +-- wb_nf_fresnel_steer.m     Shared Fresnel/USW steering vector utility
|   +-- wb_channel_gen_ofdm_nf.m  Wideband near-field OFDM channel generator
|   +-- nf_usw_steer.m            Paper B USW steering vector (dependency of wb_channel_gen_ofdm_nf)
|
+-- tests/                        8 unit test / toy test scripts
|   +-- wb_clkl_toy_test_v2.m
|   +-- wb_clkl_driver_toy_test_v2.m
|   +-- wb_crb_compressed_toy_test.m
|   +-- wb_psomp_toy_test.m
|   +-- wb_dl_omp_toy_test.m
|   +-- wb_channel_gen_toy_test.m
|   +-- bpd_toy_test.m
|   +-- mc_toy_test.m
|
+-- diagnostics/                  1 forensic diagnostic script
|   +-- wb_clkl_forensic_single.m Single-realisation KL trajectory diagnostic
|
+-- results/                      9 production CSV files (Phase 3 complete)
    +-- mc_snr_sweep_20260511_193038.csv
    +-- mc_bandwidth_sweep_20260520_180317.csv
    +-- mc_convergence_sweep_20260520_181757.csv
    +-- mc_convergence_Lhist_SNR0dB_20260520_181757.csv
    +-- mc_convergence_Lhist_SNR5dB_20260520_181757.csv
    +-- mc_convergence_Lhist_SNR10dB_20260520_181757.csv
    +-- mc_convergence_Lhist_SNR15dB_20260520_181757.csv
    +-- wb_crb_d2_sweep_20260520_183239.csv
    +-- wb_crb_geodiv_sweep_20260520_184631.csv
```

---

## Reproducing the paper figures

All figures are generated from the CSVs in `results/`. The table below maps each
paper figure to the script that generates it and the CSV it reads.

| Figure | Description | Generating script | Source CSV(s) |
|--------|-------------|-------------------|---------------|
| Fig 5  | Compressed-domain CRB vs SNR (d=1,2) | `src/wb_crb_compressed.m` | `wb_crb_d2_sweep_20260520_183239.csv` |
| Fig 6  | Compression gap: CRB vs N_RF (d=2, multi-path) | `src/wb_crb_compressed.m` | `wb_crb_d2_sweep_20260520_183239.csv` |
| Fig 7  | Geometric diversity: CRB vs M and B | `src/wb_crb_compressed.m` | `wb_crb_geodiv_sweep_20260520_184631.csv` |
| Fig 8a | RMSE vs SNR (full range, all methods) | `src/run_monte_carlo_paperC.m` | `mc_snr_sweep_20260511_193038.csv` |
| Fig 8b | NMSE vs bandwidth (full range, all methods) | `src/run_monte_carlo_paperC.m` | `mc_bandwidth_sweep_20260520_180317.csv` |
| Fig 10 | Convergence: KL objective vs iteration, 4 SNR | `src/run_monte_carlo_paperC.m` | `mc_convergence_sweep_20260520_181757.csv`, `mc_convergence_Lhist_SNR*dB_20260520_181757.csv` |

**Key verified numerical claims (all CSV-sourced):**
- B4 (WB-CL-KL) RMSE_r at SNR=10 dB: 0.300 m vs CRB 0.344 m (ratio 0.87x)
- B4 NMSE_r at B=400 MHz: -32.74 dB
- B4 convergence rate (full range, SNR=10 dB): 82.5%
- Total wideband CRB gain at B=400 MHz: +27.793 dB
  (data diversity: +27.093 dB; geometric diversity: +0.701 dB)
- Compression gap N_RF=8 vs N_RF=32 at SNR=10 dB: 8.12 dB

---

## Method labels

| Label | Method | Architecture |
|-------|--------|--------------|
| B1 | WB-BPD | Full-array (M RF chains) |
| B2 | WB-P-SOMP | Hybrid compressed (N_RF RF chains) |
| B4 | WB-CL-KL (proposed) | Hybrid compressed (N_RF RF chains) |
| B5 | WB-DL-OMP | Full-array (M RF chains) |

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

## System parameters (default production run)

| Parameter | Value | Description |
|-----------|-------|-------------|
| f_c | 28 GHz | Carrier frequency |
| M | 64 | Number of ULA antenna elements |
| N_RF | 8 | Number of RF chains |
| N | 64 | Number of pilot snapshots |
| d | 1 | Element spacing (x lambda/2) |
| r_RD | 21.26 m | Rayleigh distance |
| B | 400 MHz | Bandwidth (default) |
| K_s | 64 | Number of OFDM subcarriers |
| N_MC | 200 | Monte Carlo trials per SNR point |
| SNR range | -5 to 15 dB | 5-point sweep |

---

## License

MIT License. See [LICENSE](LICENSE).

This code is for academic research only. If you use it, please cite the paper above.
