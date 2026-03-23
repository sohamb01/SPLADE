# SPLADE
This repository contains the official implementation for the paper: **Fast localization of anomalous patches in spatial data under dependence.**

**Authors:** S. Bonnerjee (University of Chicago), S. Karmakar (University of Florida), G. Michailidis (UCLA).

---

## Overview

SPLADE (**SP**atial patch **L**ocalization of **A**nomalies under **DE**pendence) is a two-stage, scalable algorithm for detecting and localizing multiple axis-aligned rectangular anomalous patches in 2D spatial grids under general spatial dependence. The method is provably accurate, achieves near-linear runtime O(n) in favorable signal regimes, and accommodates heavy-tailed noise under only mild moment assumptions — substantially generalizing the existing literature which largely assumes independence or restrictive m-dependence structures. The following implementation tackles the 2-dimensional case ($d=2$). 

---

## Repository Structure
```
.
├── multiple_cp_helper.R     # Core Algorithm 1 and Algorithm 2 implementations
├── multiple-final-RI.R      # Simulation driver: evaluation over rho/jump/gridsize
├── plot_helper.R            # Visualization utilities (ggplot2-based)
├── LRV_helper.R             # Long-run variance (LRV) estimation
└── quantile_calc_simu.R     # Monte Carlo threshold calibration for Algorithm 2
```

---

## File Descriptions

### `multiple_cp_helper.R`
The central implementation file. Contains:

- **`prefix_sum2d(M)`** — Computes the 2D prefix sum (integral image) of a matrix, enabling O(1) rectangle sum queries.
- **`rect_sum2d(cs, s, t)`** — Returns the sum over a rectangle with corners `s` and `t` using pre-computed prefix sums.
- **`label_components(mask, connectivity)`** — BFS-based connected-component labeling on a logical 2D mask, supporting 4-connectivity (rook) or 8-connectivity (queen).
- **`single_sp_changepoint(X, alpha, kappa, C, connectivity)`** — Implements **Algorithm 1**: single anomalous patch localization via intelligent sub-sampling. Takes a 2D data matrix and hyperparameters `alpha` (subsampling exponent) and `kappa` (band-width exponent), and returns the estimated bounding box of the patch.
- **`multiple_sp_changepoints(X, alpha, kappa, Q, c, connectivity, ...)`** — Implements **Algorithm 2** (SPLADE): multiple anomalous patch localization. Performs block-mean thresholding, connected-component extraction, and then applies Algorithm 1 within each candidate region. Returns the estimated number of patches `K_hat`, their bounding boxes `I_hat` and `I_bounds`, and the search-window bounds `D_bounds`.
- **`plot_matrix(M, ...)`** — Plots a matrix using base R `image()` with matrix-consistent row/column orientation.
- **`draw_rect_indices(i1, i2, j1, j2, nrow_total, ...)`** — Overlays a rectangle on an existing `plot_matrix` call using row/column index coordinates.

### `multiple-final-RI.R`
End-to-end simulation driver used to reproduce the paper's Table 1 (Config 1, SAR noise). Loops over grid sizes (`N ∈ {250, 500}`), jump sizes, and SAR spatial correlation levels (`rho ∈ {0, 0.25, 0.5, 0.75}`). For each configuration it runs `totaliter = 100` parallel replicates using `foreach`/`doSNOW` and records:
- `K_hat`: estimated number of patches
- `symdiff_max`: maximum normalized symmetric difference between true and estimated patches (when `K_hat` is correct)
- `RI`: Adjusted Rand Index between true and estimated pixel-level labeling

Also contains the data-generating helpers:
- **`simulate.sar(N, rho, spatial.var)`** — Generates an N×N spatial autoregressive (SAR) field with correlation parameter `rho` via iterative neighbor-averaging until convergence.
- **`simulate.fulldata(N, spatial.var)`** — Generates an i.i.d. Gaussian N×N field (baseline, rho=0).
- **`epi_mean(tau1, tau2, delta, N)`** — Constructs the N×N mean matrix for K patches, specified as fractional endpoint pairs `tau1`, `tau2` with jump sizes `delta`.
- **`rects_to_label(tau1, tau2, N)`** — Converts fractional patch endpoints to an N×N integer label matrix (0 = background, k = patch k), for Rand Index computation.
- **`rect_symdiff_area_bltr(bl1, tr1, bl2, tr2)`** — Computes the symmetric-difference area between two axis-aligned rectangles specified by fractional (bottom-left, top-right) corner pairs.

### `plot_helper.R`
ggplot2-based plotting utilities for visualizing detection results:

- **`plot_detection_summary(X, tau1, tau2, res, line_width)`** — Produces a publication-ready heatmap of the data matrix with overlaid true patch boundaries (blue) and estimated patch boundaries (red). Accepts a result object from `multiple_sp_changepoints`.
- **`indices_to_rect(i1, i2, j1, j2, nr, nc)`** — Converts integer matrix index bounds to fractional coordinates for ggplot2 `geom_rect`.
- **`matrix_to_df(M)`** — Reshapes a matrix into a long data frame suitable for `geom_raster`.

### `LRV_helper.R`
Long-run variance (LRV) estimation, required for calibrating the first-stage threshold Q in Algorithm 2:

- **`lrv_estimator_outer_shell(X, alpha, Bn_pow, kernel)`** — Estimates σ² using observations from the outer border band of the grid (which is guaranteed to be anomaly-free under Assumption 3), with an Epanechnikov kernel HAC-style estimator. This avoids bias from anomalous pixels.
- **`lrv_estimator_fft_new(X, alpha, Bn_pow, kernel)`** — An FFT-accelerated version of the same LRV estimator. Computes the kernel-weighted autocovariance sum via zero-padded 2D FFT, substantially faster for large grids.
- **`border_band(n1, n2, w_row, w_col)`** — Returns the pixel indices forming the outer border band of an n1×n2 grid (default thickness `ceil(sqrt(n_k))`), used by `lrv_estimator_outer_shell`.
- **`K_epanechnikov(u)`** — Evaluates the 1D Epanechnikov kernel K(u) = (1 − u²) · 1(|u| ≤ 1).

### `quantile_calc_simu.R`
Monte Carlo calibration of the first-stage detection threshold Q:

- **`simulate_Q_blockmeans(n, alpha, sigma2, nsim, seed, alpha_q, progress)`** — Simulates `nsim` i.i.d. Gaussian null fields of size n₁×n₂, computes the maximum absolute block mean over all Algorithm-2-style blocks, and returns the (1 − α)-quantile of this null distribution as the threshold Q. Supports heteroskedastic fields via a spatially varying `sigma2` matrix. In practice, this is called once with `sigma2 = 1` and the result is rescaled by the estimated `sqrt(lrv_est)` at runtime.
- **`simulate.sar(N, rho, spatial.var)`** — (also defined here for standalone use) SAR field generator.
- **`rect_sums_from_cs(cs, i1, i2, j1, j2)`** — Vectorized rectangle sum queries from prefix sums, used internally for fast block-mean computation during threshold simulation.

---

## Method Summary

SPLADE operates in two stages:

**Stage 1 — Block screening.** The grid is partitioned into blocks of side length ~n^α. Block means are computed and compared against a threshold Q (calibrated via `simulate_Q_blockmeans` and scaled by `sqrt(lrv_est)`). Blocks whose absolute mean exceeds Q are flagged; connected components of flagged blocks with sufficient size form candidate regions C₁, …, C_K̂. The number of components K̂ is the estimated patch count.

**Stage 2 — Refined localization.** Each candidate region Cⱼ is expanded into a search window Dⱼ. Algorithm 1 (intelligent sub-sampling) is applied independently within each Dⱼ: a coarse sub-sampled estimate identifies approximate patch boundaries, which are then refined via a local scan over a narrow band, yielding the final estimate Î_j.

Under mild dependence assumptions (finite p-th moments, long-run variance, uniform Gaussian approximation), SPLADE consistently estimates both K and each Î_j at the minimax-optimal rate, in approximately O(n) time.

---

## Installation and Dependencies

No package installation is required for the core algorithm. All computations in `multiple_cp_helper.R` use **base R** only.

The simulation driver and helpers require:
```r
install.packages(c("parallel", "foreach", "doSNOW", "mclust",
                   "ggplot2", "viridis", "tidyr", "dplyr"))
```

---

## Quick Start
```r
source("multiple_cp_helper.R")
source("LRV_helper.R")
source("quantile_calc_simu.R")

set.seed(42)
N <- 256
# Define two anomalous patches in fractional coordinates
tau1 <- rbind(c(0.20, 0.20), c(0.60, 0.60))
tau2 <- rbind(c(0.45, 0.70), c(0.85, 0.85))
delta <- c(1.0, -1.0)

# Simulate a SAR(rho=0.4) field with two anomalous patches
mu  <- epi_mean(tau1, tau2, delta, N)
Z   <- simulate.sar(N = N, rho = 0.4)
X   <- Z + mu

# Stage 0: calibrate threshold (run once per grid size / alpha)
Q_base <- simulate_Q_blockmeans(n = c(N, N), alpha = 0.5, alpha_q = 0.5,
                                sigma2 = 1, nsim = 1000, seed = 123)

# Stage 0: estimate long-run variance from border band
lrv_est <- lrv_estimator_outer_shell(X, alpha = 0.5, Bn_pow = 1/3)
Q       <- Q_base$Q * sqrt(lrv_est)

# Run SPLADE
res <- multiple_sp_changepoints(X, alpha = 0.5, kappa = 0.01, Q = Q,
                                connectivity = 8L,
                                plot_components = TRUE,
                                plot_estimates  = TRUE)

cat("Estimated number of patches:", res$K_hat, "\n")
print(res$I_bounds)   # rows are (i1, j1, i2, j2) for each estimated patch
```

### Key hyperparameters

| Parameter | Role | Recommended default |
|-----------|------|---------------------|
| `alpha` | Block size exponent: L_k = floor(n_k^alpha) | 0.5 |
| `kappa` | Band-width exponent in Algorithm 1 | 0.01 |
| `alpha_q` | Quantile level for threshold Q | 0.5 |
| `nsim` | Monte Carlo simulations for Q | 1000 |
| `Bn_pow` | LRV bandwidth exponent | 1/3 |
| `connectivity` | CC labeling: 4 (rook) or 8 (queen) | 8 |

---

## Reproducing Simulation Results

The full simulation (Table 1 of the paper) is run via `multiple-final-RI.R`. Edit the `num_core`, `totaliter`, grid sizes, jump sizes, and `rho` values at the top of the file as needed, then:
```r
source("multiple-final-RI.R")
```

Results are collected in the data frame `datfin`, which is printed at the end and contains columns `N`, `jump`, `cor`, `norm_max_sym_diff`, `RI`, and `counts_K_hat`.

---

## Citation
```bibtex
@article{bonnerjee2025splade,
  title   = {Fast and provably accurate detection of anomalous patches for spatial data},
  author  = {Bonnerjee, S. and Karmakar, S. and Michailidis, G.},
  journal = {Preprint},
  year    = {2025},
}
```

## 3D Extension: Fibre System Analysis

The following files extend SPLADE to **three-dimensional** spatial arrays (volumetric grids), targeting the fibre-reinforced polymer application described in Appendix G of the paper. The algorithms are direct 3D analogues of Algorithms 1 and 2, replacing 2D rectangles with axis-aligned rectangular prisms and using 3D prefix sums and 3D FFT-based LRV estimation.

### Additional Files
```
├── multiple_cp_helper_3d.R    # Core 3D Algorithm 1 and Algorithm 2
├── LRV_helper_3d.R            # 3D FFT-based long-run variance estimator
├── quantile_calc_simu_3d.R    # Monte Carlo threshold calibration (3D)
├── run_one_simulation_3d.R    # Self-contained 3D simulation example
└── fibre_anomaly_main.R       # Real-data driver for fibre direction datasets
```

### `multiple_cp_helper_3d.R`
3D counterpart of `multiple_cp_helper.R`. Contains:

- **`prefix_sum3d(A)`** — Builds the 3D summed-volume table (prefix sum) of an n1×n2×n3 array via three sequential cumsum passes, enabling O(1) box sum queries via inclusion-exclusion over 8 corners.
- **`box_sum3d(cs, s, t)`** — Returns the sum over the rectangular prism with corners `s = (i1,j1,k1)` and `t = (i2,j2,k2)` from the 3D prefix array.
- **`box_sums_from_cs3d(cs, i1, i2, j1, j2, k1, k2)`** — Vectorized version of `box_sum3d` for simultaneous computation over many boxes, used in threshold simulation.
- **`single_sp_changepoint_3d(X, alpha, kappa, ...)`** — Implements **Algorithm 1** in 3D: sub-samples the volume on a coarse grid, finds a preliminary prism estimate, then refines within local band neighborhoods around each face.
- **`multiple_sp_changepoints_3d(X, alpha, kappa, Q, connectivity, ...)`** — Implements **Algorithm 2** in 3D (SPLADE-3D): partitions the volume into blocks, flags blocks exceeding threshold Q, extracts 3D connected components (supporting 6-, 18-, or 26-connectivity), and applies `single_sp_changepoint_3d` within each candidate region.
- **`epi_mean_3d(tau1, tau2, delta, N)`** — Constructs the n1×n2×n3 mean array for K rectangular prism patches specified by fractional endpoints `tau1`, `tau2` (K×3 matrices) and jump sizes `delta`.
- **`simulate.sar.3d(n, rho, spatial.var)`** — Generates a 3D SAR field of size `n = c(n1, n2, n3)` with spatial correlation `rho` via iterative 6-neighbor averaging until convergence.

### `LRV_helper_3d.R`
Long-run variance estimation for 3D fields:

- **`lrv_estimator_fft_3d(X, alpha, Bn_pow, kernel)`** — Estimates σ² from a 3D array using zero-padded 3D FFT autocorrelation with a separable Epanechnikov product kernel. Block-mean centering is applied first (using the same block construction as Algorithm 2) to reduce bias from the anomalous region. The bandwidth vector is `Bn = (n1^Bn_pow, n2^Bn_pow, n3^Bn_pow)`.
- **`K_epanechnikov(u)`** — 1D Epanechnikov kernel (shared with 2D version).

### `quantile_calc_simu_3d.R`
Monte Carlo threshold calibration for 3D block means:

- **`simulate_Q_blockmeans_3d(n, alpha, sigma2, nsim, seed, alpha_q, progress)`** — Direct 3D analogue of `simulate_Q_blockmeans`. Simulates `nsim` Gaussian null volumes, computes the maximum absolute block mean over all Algorithm-2-style prism blocks (using `box_sums_from_cs3d`), and returns the (1 − α)-quantile as Q. Accepts either a scalar or a 3D array for `sigma2`.

### `run_one_simulation_3d.R`
Self-contained end-to-end script for a synthetic 3D experiment. Demonstrates the full pipeline:
1. Define K=3 rectangular prism patches via fractional endpoints `tau1`, `tau2`
2. Simulate a 3D SAR(ρ) noise field with `simulate.sar.3d`
3. Estimate the LRV with `lrv_estimator_fft_3d`
4. Calibrate Q with `simulate_Q_blockmeans_3d` (both estimated and oracle LRV)
5. Run `multiple_sp_changepoints_3d` and evaluate

### `fibre_anomaly_main.R`
Real-data driver for the glass fibre-reinforced polymer analysis (Appendix G). Reads MAVI-processed fibre direction CSV data and runs the full 3D SPLADE pipeline on a chosen direction attribute. Steps:
1. **Load & preprocess** — calls `load_fibre_data()` (from `preprocess_fibre_data.R`) to read the CSV and populate six direction arrays (`x_abs`, `y_abs`, `z_abs`, etc.) on the voxel grid
2. **Select attribute** — one of `"x"`, `"y"`, `"z"`, `"x_abs"`, `"y_abs"`, `"z_abs"`, or `"aniso"`
3. **LRV estimation** — `lrv_estimator_fft_3d`
4. **Threshold calibration** — `simulate_Q_blockmeans_3d` with the estimated LRV
5. **Detection** — `multiple_sp_changepoints_3d` with 26-connectivity
6. **Output** — prints estimated patch bounding boxes in both grid indices and fractional coordinates

### Quick Start (3D)
```r
source("multiple_cp_helper_3d.R")
source("LRV_helper_3d.R")
source("quantile_calc_simu_3d.R")

N   <- c(60, 60, 60)
rho <- 0.3
tau1 <- rbind(c(0.10, 0.10, 0.10), c(0.50, 0.10, 0.55))
tau2 <- rbind(c(0.30, 0.35, 0.35), c(0.75, 0.35, 0.85))
delta <- c(0.8, -0.8)

mu  <- epi_mean_3d(tau1, tau2, delta, N)
Z   <- simulate.sar.3d(n = N, rho = rho, spatial.var = 1)
X   <- Z + mu

lrv_est <- lrv_estimator_fft_3d(X, alpha = 0.5, Bn_pow = 1/3)
resQ    <- simulate_Q_blockmeans_3d(n = N, alpha = 0.5, sigma2 = lrv_est,
                                    nsim = 500, seed = 123)

res <- multiple_sp_changepoints_3d(X, alpha = 0.5, kappa = 0.01,
                                   Q = resQ$Q, connectivity = 26L)
cat("K_hat:", res$K_hat, "\n")
print(res$I_bounds)
```

### Fibre Data Analysis
```r
# Requires: preprocess_fibre_data.R, plot_helper_3d.R
source("fibre_anomaly_main.R")
# Edit 'attribute' at the top of fibre_anomaly_main.R to switch
# between fibre directions: "x_abs", "y_abs", "z_abs", etc.
