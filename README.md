# ESPRIT for pulses of unknown shape from multiple snapshots

MATLAB code for

> M. Kalra and K. Lee, "Stable Estimation of Pulses of Unknown Shape From Multiple Snapshots via ESPRIT," *IEEE Transactions on Signal Processing*, vol. 72, pp. 2637–2648, 2024. doi: [10.1109/TSP.2024.3403494](https://doi.org/10.1109/TSP.2024.3403494)

## Contents

- `esprit_pulse_estimation.m` — self-contained implementation with a short example: the ESPRIT estimator, the two sub-array designs (maximum overlapping and random doublets), the matching-distance error metric, and data generation. Run it directly.
- `paper_figures/` — scripts for the numerical results in the paper.

## Figures

Run the scripts from inside `paper_figures/`. Scripts with saved results have a flag `use_saved` at the top: `true` plots the results in `data/`, `false` reruns the Monte Carlo simulation.

| Figure | Script | Saved results |
|---|---|---|
| Fig. 1 (left) | `Fig1a_phase_transition_Dirac.m` | `data/Fig1a_phase_transition_Dirac.mat` |
| Fig. 1 (right) | `Fig1b_phase_transition_nonDirac.m` | `data/Fig1b_phase_transition_nonDirac.mat` |
| Fig. 2(a) | `Fig2a_error_vs_L.m` | `data/Fig2a_error_vs_L.mat` |
| Fig. 2(b) | `Fig2b_error_vs_M.m` | — |
| Fig. 2(c) | `Fig2c_error_vs_pulse_width.m` | — |
| Fig. 2(d) | `Fig2d_error_vs_kappa_RX.m` | — |
| Fig. 3(a) | `Fig3a_error_vs_Mtilde.m` | `data/Fig3a_error_vs_Mtilde_nonDirac.mat`, `data/Fig3a_error_vs_Mtilde_Dirac.mat` |
| Fig. 3(b) | `Fig3b_phase_transition_M_vs_S.m` | `data/Fig3b_phase_transition_M_vs_S.mat` |

Shared functions in `paper_figures/`:

- `Esprit.m` — ESPRIT estimator
- `design_subarray.m` — sub-array designs (mode 1: maximum overlapping, mode 2: random doublets)
- `data_gen_new_pulse.m` — pulse spectrum and Fourier matrix
- `location_error.m` — matching distance between true and estimated locations
- `tau_k_gen.m` — pulse location generator

## Requirements

MATLAB R2016b or newer. `sinc` is part of the Signal Processing Toolbox.

## Citation

```bibtex
@article{kalra2024stable,
  author  = {Kalra, Meghna and Lee, Kiryung},
  title   = {Stable Estimation of Pulses of Unknown Shape From Multiple Snapshots via {ESPRIT}},
  journal = {IEEE Transactions on Signal Processing},
  volume  = {72},
  pages   = {2637--2648},
  year    = {2024},
  doi     = {10.1109/TSP.2024.3403494}
}
```
