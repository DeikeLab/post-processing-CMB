# Interface statistics

Windowed statistics of the air–water interface in the wind-wave simulations: moments of the
elevation $\eta(x,z,t)$, moments of the surface slopes $\partial_x\eta$ and $\partial_z\eta$,
and the running PDFs of the normalised elevation
$\xi = (\eta - \langle\eta\rangle)/\sigma_\eta$.

Unlike the spectral notebooks, this one does **not** read from `/projects/DEIKE/...` at run
time. The statistics are accumulated on the fly by the solver into small text tables, and
`sync_case.sh` copies those tables off `/scratch` into `data/` here — a few MB per case
against ~450 GB of dumps and restarts left behind. So the notebook runs anywhere, including
off the cluster.

## Files

- **eta_stats_analysis.ipynb** — the analysis. Windowed $\langle\eta\rangle$, $\sigma_\eta$,
  skewness and kurtosis; the same moments for both slope components; a twin-axis cross-check
  of $\gamma_1(\eta)$ against $\kappa(\partial_x\eta)$; per-window elevation PDFs against a
  reference Gaussian, and an overlay of all windows on a log axis.
- **generate_eta_plots.py** — the same figures headless, called from the run script after
  every job so `statistics/plots/` stays current on the cluster:
  `python generate_eta_plots.py <case>/statistics`
- **sync_case.sh** — copy a case's statistics into `data/`:
  `./sync_case.sh /scratch/gpfs/DEIKE/cm6797/multiphase_cases/<case>`

## Cases in `data/`

| case | Re_τ | Bo | k_p H_s | u_*/c | Re_wave | level |
| --- | --- | --- | --- | --- | --- | --- |
| `re720_bo200_kpHs0p16_uoc0p50_reW2.5e4_L10_stats` | 720 | 200 | 0.161 | 0.50 | 2.56e4 | 10 |

Each case folder holds `case_params.txt` (the parameter header sliced out of the multi-GB
`out.log`, which carries `T0` = $T_p$), `statistics/` with the tables, the per-rank
accumulator checkpoints and the generated PNGs, `logs/`, and the SLURM script.

The notebook locates its data automatically: `data/<case>/` when run from this folder,
`./statistics/` when the notebook is copied into a case directory on the cluster. Set
`CASE_DIR` in section 1 to pick a specific case.

## Reading `eta_stats_window.out`

Worth knowing before quoting any number from it: the file is an append-only trace of the
*running* window accumulator, not one row per window. Each dump appends the moments
accumulated so far in the window currently being filled; `is_partial = 1` marks an
in-progress row and `is_partial = 0` the closing row of a window. The case above has 228
rows but only 9 closed windows, and restarts replay rows already written (73 exact
duplicates here). The notebook drops the duplicates on load and marks window closures;
`COMPLETE_ONLY = True` reduces it to one point per finished window. Partial rows can be
averages over as few as a couple of hundred snapshots, so they carry early-window noise
that is easy to mistake for a physical transient.

Columns:

| col | quantity | col | quantity |
| --- | --- | --- | --- |
| 0 | `t_end` | 11 | κ(ξ) kurtosis |
| 1 | `i` (step) | 12 | window weight |
| 2 | `t_rel` | 13–16 | ⟨s_x⟩, σ_sx, γ₁(s_x), κ(s_x) |
| 3 | window id | 17–20 | ⟨s_z⟩, σ_sz, γ₁(s_z), κ(s_z) |
| 4 | `is_partial` | | |
| 5 | `nsnap` | | with s_x = ∂η/∂x, s_z = ∂η/∂z |
| 6–7 | `t_start`, `t_last` | | |
| 8–9 | ⟨η⟩, σ_η | | |
| 10 | γ₁(ξ) skewness | | |

Runs predating the slope diagnostics write only the first 13 columns; the notebook detects
this and skips the slope sections. Slope moments are normalised over the interface area with
$|n_y| > 0.10$, so near-vertical cells are excluded from the slope statistics.

`eta_pdf_window_*.out` has columns
`0:eta_center 1:xi_center 2:pdf_xi 3:pdf_gauss 4:count 5:pdf_eta`; column 3 is the unit
normal in $\xi$, so the departure from it reads directly as non-Gaussianity of the surface.
