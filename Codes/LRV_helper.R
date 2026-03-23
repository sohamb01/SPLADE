######## LRV estimation ###################

# --- Kernel: Epanechnikov (1D) with support [-1,1] ---

prefix_sum2d <- function(M) {
  cs <- apply(M, 2, cumsum)
  cs <- t(apply(cs, 1, cumsum))
  cs
}

# Fast rectangle sum from prefix sums; s=(i1,j1), t=(i2,j2), inclusive corners
rect_sum2d <- function(cs, s, t) {
  i1 <- s[1]; j1 <- s[2]; i2 <- t[1]; j2 <- t[2]
  if (i1 > i2) { tmp <- i1; i1 <- i2; i2 <- tmp }
  if (j1 > j2) { tmp <- j1; j1 <- j2; j2 <- tmp }
  tot   <- cs[i2, j2]
  left  <- if (j1 > 1) cs[i2, j1 - 1] else 0
  up    <- if (i1 > 1) cs[i1 - 1, j2] else 0
  upL   <- if (i1 > 1 && j1 > 1) cs[i1 - 1, j1 - 1] else 0
  tot - left - up + upL
}

# --- Epanechnikov kernel (1D), support [-1, 1] ---
K_epanechnikov <- function(u) {
  out <- numeric(length(u))
  idx <- abs(u) <= 1
  out[idx] <-  (1 - u[idx]^2)  # 3/4 * (1 - u^2), |u|<=1
  out
}

# --- 2D FFT helpers using mvfft (true 2D) ---
fft2  <- function(Z) { mvfft(t(mvfft(Z))) |> t() }
ifft2 <- function(Z) {
  Z <- mvfft(t(mvfft(Z, inverse = TRUE)), inverse = TRUE) |> t()
  Z / (nrow(Z) * ncol(Z))
}

.make_groups <- function(n, L) {
  m <- ceiling(n / L)
  lens <- rep(L, m)
  if (m > 1L) lens[m] <- n - sum(lens[1:(m-1)]) else lens[m] <- n
  rep.int(seq_len(m), times = lens)
}


lrv_estimator_fft_new <- function(X, alpha, Bn_pow = 1/3, kernel = c("epanechnikov"),
                                  show_progress = TRUE) {
  kernel <- match.arg(kernel)
  stopifnot(is.matrix(X), is.numeric(X))
  n1 <- nrow(X); n2 <- ncol(X)
  n_tot <- n1 * n2
  Bn <- c(n1^Bn_pow, n2^Bn_pow)  # anisotropic bandwidth (vector), same as naive
  
  # ------------------------------
  # 1) Block means & centering (Alg. steps 1–9) — identical to naive
  # ------------------------------
  L1 <- max(1L, floor(n1^alpha)); L2 <- max(1L, floor(n2^alpha))
  m1 <- ceiling(n1 / L1);         m2 <- ceiling(n2 / L2)
  
  I1 <- function(s) { a <- (s - 1L) * L1 + 1L; b <- min(s * L1, n1); a:b }
  I2 <- function(s) { a <- (s - 1L) * L2 + 1L; b <- min(s * L2, n2); a:b }
  
  g1 <- .make_groups(n1, L1)     # length n1, values in {1,...,m1}
  g2 <- .make_groups(n2, L2)     # length n2, values in {1,...,m2}
  
  # (A) block sums via two-stage grouping (C-level, fast)
  # first aggregate rows by g1  -> m1 x n2
  row_aggr   <- rowsum(X, g1)                    
  # then aggregate columns by g2 (work on transpose) -> m2 x m1, transpose back
  block_sums <- t(rowsum(t(row_aggr), g2))       # m1 x m2
  
  # (B) block sizes so we get exact means even for edge blocks
  nrows_per_block <- as.integer(tabulate(g1, nbins = m1))  # length m1
  ncols_per_block <- as.integer(tabulate(g2, nbins = m2))  # length m2
  block_sizes     <- outer(nrows_per_block, ncols_per_block, `*`)  # m1 x m2
  
  # (C) block means
  block_mean <- block_sums / block_sizes                      # m1 x m2
  
  # (D) broadcast block means back to full grid to form Xbar
  #     (exactly the same assignment Xbar[I1(s1), I2(s2)] <- block_mean[s1,s2])
  Xbar <- block_mean[cbind(rep(g1, times = n2), rep(g2, each = n1))]
  dim(Xbar) <- c(n1, n2)
  
  Y <- X - Xbar  # centered residual field
  
  # ------------------------------
  # 2) Linear autocorrelation via zero-padded 2D FFT
  #    (matches naive "valid overlap" sums)
  # ------------------------------
  nextpow2 <- function(x) 2^ceiling(log2(x))
  M1 <- nextpow2(2 * n1 - 1L)
  M2 <- nextpow2(2 * n2 - 1L)
  
  Ypad <- matrix(0.0, M1, M2)
  Ypad[1:n1, 1:n2] <- Y
  
  FY <- fft2(Ypad)
  R_full <- Re(ifft2(FY * Conj(FY)))  # zero-lag at [1,1]; no wrap thanks to padding
  
  # ------------------------------
  # 3) Kernel-weighted lag sum (same lag set & weights as naive)
  #    NOTE: the naive code defines H1/H2 using min(floor(Bn), n_k-1),
  #    where floor(Bn) is a *vector*. We replicate that behavior exactly.
  # ------------------------------
  H1 <- min(floor(Bn), n1 - 1L)
  H2 <- min(floor(Bn), n2 - 1L)
  h1_vals <- (-H1):H1
  h2_vals <- (-H2):H2
  
  if (kernel == "epanechnikov") {
    k1 <- K_epanechnikov(h1_vals / Bn[1])
    k2 <- K_epanechnikov(h2_vals / Bn[2])
  } else {
    stop("Only 'epanechnikov' is implemented.")
  }
  K2D <- outer(k1, k2, `*`)
  
  # Map integer lags to indices in R_full using modulo on the padded sizes
  rows <- ((h1_vals %% M1) + 1L)
  cols <- ((h2_vals %% M2) + 1L)
  R_sub <- R_full[rows, cols]
  
  Sigma_hat <- as.numeric(sum(K2D * R_sub) / n_tot)
  Sigma_hat
}


