###############################################
# LRV estimation for 3D spatial random fields
###############################################

# --- Epanechnikov kernel (1D), support [-1, 1] ---
K_epanechnikov <- function(u) {
  out <- numeric(length(u))
  idx <- abs(u) <= 1
  out[idx] <- (1 - u[idx]^2)
  out
}

# --- Helper: group indices for block aggregation ---
.make_groups <- function(n, L) {
  m <- ceiling(n / L)
  lens <- rep(L, m)
  if (m > 1L) lens[m] <- n - sum(lens[1:(m - 1)]) else lens[m] <- n
  rep.int(seq_len(m), times = lens)
}

# --- 3D FFT-based LRV estimator ---
# Uses block-mean centering + zero-padded 3D FFT autocorrelation + product kernel.
lrv_estimator_fft_3d <- function(X, alpha, Bn_pow = 1/3, kernel = c("epanechnikov"),
                                 show_progress = TRUE) {
  kernel <- match.arg(kernel)
  stopifnot(length(dim(X)) == 3L, is.numeric(X))
  d <- dim(X)
  n1 <- d[1]; n2 <- d[2]; n3 <- d[3]
  n_tot <- n1 * n2 * n3
  Bn <- c(n1^Bn_pow, n2^Bn_pow, n3^Bn_pow)

  # 1) Block means & centering
  L1 <- max(1L, floor(n1^alpha)); L2 <- max(1L, floor(n2^alpha)); L3 <- max(1L, floor(n3^alpha))
  g1 <- .make_groups(n1, L1)
  g2 <- .make_groups(n2, L2)
  g3 <- .make_groups(n3, L3)
  m1 <- max(g1); m2 <- max(g2); m3 <- max(g3)

  # Block sums via sequential aggregation along each dimension
  # Reshape to matrix, aggregate dim 1
  Xmat <- matrix(X, nrow = n1, ncol = n2 * n3)
  agg1 <- rowsum(Xmat, g1)  # m1 x (n2*n3)

  # Aggregate dim 2: reshape so dim2 is rows
  agg1_arr <- array(agg1, dim = c(m1, n2, n3))
  agg12 <- array(0, dim = c(m1, m2, n3))
  for (k in 1:n3) {
    agg12[, , k] <- t(rowsum(t(agg1_arr[, , k]), g2))
  }

  # Aggregate dim 3: reshape so dim3 is rows
  agg123 <- array(0, dim = c(m1, m2, m3))
  for (i in 1:m1) {
    agg123[i, , ] <- t(rowsum(t(agg12[i, , ]), g3))
  }

  # Block sizes
  n1_per <- tabulate(g1, nbins = m1)
  n2_per <- tabulate(g2, nbins = m2)
  n3_per <- tabulate(g3, nbins = m3)
  block_sizes <- outer(outer(n1_per, n2_per, `*`), n3_per, `*`)  # m1 x m2 x m3

  block_mean <- agg123 / block_sizes

  # Broadcast block means back to full grid
  # R arrays are column-major: fastest-varying index is dim1, then dim2, then dim3
  # For an n1 x n2 x n3 array, linear index = i + (j-1)*n1 + (k-1)*n1*n2
  idx1 <- rep(g1, times = n2 * n3)                    # repeats g1 for each (j,k)
  idx2 <- rep(rep(g2, each = n1), times = n3)          # repeats g2-expanded for each k
  idx3 <- rep(g3, each = n1 * n2)                      # repeats g3 for each (i,j)
  Xbar <- array(block_mean[cbind(idx1, idx2, idx3)], dim = c(n1, n2, n3))
  Y <- X - Xbar

  # 2) Zero-padded 3D FFT autocorrelation
  nextpow2 <- function(x) 2^ceiling(log2(x))
  M1 <- nextpow2(2 * n1 - 1L)
  M2 <- nextpow2(2 * n2 - 1L)
  M3 <- nextpow2(2 * n3 - 1L)

  Ypad <- array(0.0, dim = c(M1, M2, M3))
  Ypad[1:n1, 1:n2, 1:n3] <- Y

  if (show_progress) message("Computing 3D FFT for LRV...")
  FY <- stats::fft(Ypad)
  R_full <- Re(stats::fft(FY * Conj(FY), inverse = TRUE)) / (M1 * M2 * M3)

  # 3) Kernel-weighted lag sum
  H1 <- min(floor(Bn[1]), n1 - 1L)
  H2 <- min(floor(Bn[2]), n2 - 1L)
  H3 <- min(floor(Bn[3]), n3 - 1L)
  h1_vals <- (-H1):H1
  h2_vals <- (-H2):H2
  h3_vals <- (-H3):H3

  k1 <- K_epanechnikov(h1_vals / Bn[1])
  k2 <- K_epanechnikov(h2_vals / Bn[2])
  k3 <- K_epanechnikov(h3_vals / Bn[3])

  # 3D product kernel
  K3D <- outer(outer(k1, k2, `*`), k3, `*`)

  # Map lags to FFT indices (modular wrap)
  rows <- ((h1_vals %% M1) + 1L)
  cols <- ((h2_vals %% M2) + 1L)
  deps <- ((h3_vals %% M3) + 1L)
  R_sub <- R_full[rows, cols, deps]

  Sigma_hat <- as.numeric(sum(K3D * R_sub) / n_tot)
  if (show_progress) message(paste0("LRV (FFT) estimate: ", round(Sigma_hat, 6)))
  Sigma_hat
}

# --- 3D outer-shell LRV estimator ---
# Uses the boundary shell of the 3D array (voxels near any face) for centering,
# then direct lag-sum computation.
border_band_3d <- function(n1, n2, n3,
                           w1 = ceiling(sqrt(n1)),
                           w2 = ceiling(sqrt(n2)),
                           w3 = ceiling(sqrt(n3))) {
  w1 <- max(0L, min(w1, floor(n1 / 2)))
  w2 <- max(0L, min(w2, floor(n2 / 2)))
  w3 <- max(0L, min(w3, floor(n3 / 2)))

  # Build full grid indices
  grid <- expand.grid(i = seq_len(n1), j = seq_len(n2), k = seq_len(n3))
  inner <- (grid$i > w1) & (grid$i <= n1 - w1) &
           (grid$j > w2) & (grid$j <= n2 - w2) &
           (grid$k > w3) & (grid$k <= n3 - w3)
  as.matrix(grid[!inner, , drop = FALSE])
}

lrv_estimator_outer_shell_3d <- function(X, alpha, Bn_pow = 1/3,
                                         kernel = c("epanechnikov"),
                                         show_progress = TRUE) {
  kernel <- match.arg(kernel)
  stopifnot(length(dim(X)) == 3L, is.numeric(X))
  d <- dim(X)
  n1 <- d[1]; n2 <- d[2]; n3 <- d[3]
  n_tot <- n1 * n2 * n3

  outer_shell <- border_band_3d(n1, n2, n3)
  nloc <- nrow(outer_shell)

  # Mean over boundary shell
  shell_idx <- cbind(outer_shell[, "i"], outer_shell[, "j"], outer_shell[, "k"])
  Xbar <- mean(X[shell_idx])

  Bn <- rep((nloc)^(Bn_pow / 3), 3)

  # Center: zero inside, centered on shell
  Y <- array(0, dim = d)
  Y[shell_idx] <- X[shell_idx] - Xbar

  # Lag-sum computation
  H1 <- min(floor(Bn[1]), n1 - 1L)
  H2 <- min(floor(Bn[2]), n2 - 1L)
  H3 <- min(floor(Bn[3]), n3 - 1L)
  h1_vals <- (-H1):H1
  h2_vals <- (-H2):H2
  h3_vals <- (-H3):H3

  k1 <- K_epanechnikov(h1_vals / Bn[1])
  k2 <- K_epanechnikov(h2_vals / Bn[2])
  k3 <- K_epanechnikov(h3_vals / Bn[3])

  acc <- 0.0
  total_steps <- length(h1_vals) * length(h2_vals) * length(h3_vals)
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = total_steps, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  step <- 0L

  for (a in seq_along(h1_vals)) {
    h1 <- h1_vals[a]
    if (h1 >= 0L) { rows <- 1:(n1 - h1) } else { rows <- (1 - h1):n1 }
    if (length(rows) == 0L) {
      step <- step + length(h2_vals) * length(h3_vals)
      if (show_progress) utils::setTxtProgressBar(pb, step)
      next
    }
    for (b in seq_along(h2_vals)) {
      h2 <- h2_vals[b]
      if (h2 >= 0L) { cols <- 1:(n2 - h2) } else { cols <- (1 - h2):n2 }
      if (length(cols) == 0L) {
        step <- step + length(h3_vals)
        if (show_progress) utils::setTxtProgressBar(pb, step)
        next
      }
      for (cc in seq_along(h3_vals)) {
        h3 <- h3_vals[cc]
        w <- k1[a] * k2[b] * k3[cc]
        step <- step + 1L
        if (w == 0) {
          if (show_progress) utils::setTxtProgressBar(pb, step)
          next
        }
        if (h3 >= 0L) { deps <- 1:(n3 - h3) } else { deps <- (1 - h3):n3 }
        if (length(deps) == 0L) {
          if (show_progress) utils::setTxtProgressBar(pb, step)
          next
        }

        A <- Y[rows, cols, deps, drop = FALSE]
        B <- Y[rows + h1, cols + h2, deps + h3, drop = FALSE]
        acc <- acc + w * sum(A * B)

        if (show_progress) utils::setTxtProgressBar(pb, step)
      }
    }
  }

  Sigma_hat <- as.numeric(acc / nloc)
  if (show_progress) message(paste0("\nLRV (outer shell) estimate: ", round(Sigma_hat, 6)))
  Sigma_hat
}
