#############################################
# Q_{1-α} via simulation of block-mean maxima
# (matching Algorithm 2 block construction)
#############################################

# 2D prefix sums (integral image)
prefix_sum2d <- function(M) {
  cs <- apply(M, 2, cumsum)
  cs <- t(apply(cs, 1, cumsum))
  cs
}

# Vectorized rectangle sums from prefix sums.
# i1,i2,j1,j2 are equal-length integer vectors (inclusive corners).
rect_sums_from_cs <- function(cs, i1, i2, j1, j2) {
  stopifnot(length(i1) == length(i2), length(j1) == length(j2))
  n <- length(i1)
  # ensure i1<=i2, j1<=j2
  ii1 <- pmin(i1, i2); ii2 <- pmax(i1, i2)
  jj1 <- pmin(j1, j2); jj2 <- pmax(j1, j2)
  
  tot  <- cs[cbind(ii2, jj2)]
  left <- numeric(n); idx <- jj1 > 1
  left[idx] <- cs[cbind(ii2[idx], jj1[idx] - 1)]
  
  up   <- numeric(n); idx <- ii1 > 1
  up[idx] <- cs[cbind(ii1[idx] - 1, jj2[idx])]
  
  upL  <- numeric(n); idx <- (ii1 > 1) & (jj1 > 1)
  upL[idx] <- cs[cbind(ii1[idx] - 1, jj1[idx] - 1)]
  
  tot - left - up + upL
}

# Main: simulate Q threshold for the max block mean
# n        : c(n1, n2) grid size
# alpha    : block-length exponent (L_k = floor(n_k^alpha));
#            also used as the quantile level -> Q_{1-alpha}.
# sigma2   : variance of the Gaussian field. Either a scalar or an n1 x n2 matrix.
# nsim     : number of Monte Carlo simulations
# seed     : optional random seed for reproducibility
# alpha_q  : (optional) use a different quantile level; default = alpha
# progress : show a simple progress bar
simulate_Q_blockmeans <- function(n, alpha, sigma2 = 1, nsim = 1000, seed = NULL,
                                  alpha_q = alpha, progress = TRUE) {
  stopifnot(length(n) == 2)
  n1 <- as.integer(n[1]); n2 <- as.integer(n[2])
  
  # --- Blocks per Algorithm 2 ---
  L1 <- max(1L, floor(n1^alpha)); L2 <- max(1L, floor(n2^alpha))
  m1 <- ceiling(n1 / L1);         m2 <- ceiling(n2 / L2)
  
  # Precompute all block rectangles (inclusive corners) and areas
  i1_vec <- integer(m1 * m2)
  i2_vec <- integer(m1 * m2)
  j1_vec <- integer(m1 * m2)
  j2_vec <- integer(m1 * m2)
  area   <- integer(m1 * m2)
  
  idx <- 0L
  for (s1 in 1:m1) for (s2 in 1:m2) {
    idx <- idx + 1L
    i1 <- (s1 - 1L) * L1 + 1L
    i2 <- min(s1 * L1, n1)
    j1 <- (s2 - 1L) * L2 + 1L
    j2 <- min(s2 * L2, n2)
    i1_vec[idx] <- i1; i2_vec[idx] <- i2
    j1_vec[idx] <- j1; j2_vec[idx] <- j2
    area[idx]   <- (i2 - i1 + 1L) * (j2 - j1 + 1L)
  }
  
  # --- Simulation loop ---
  if (!is.null(seed)) set.seed(seed)
  maxima <- numeric(nsim)
  if (progress) pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  
  use_scalar_var <- length(sigma2) == 1L
  if (!use_scalar_var) {
    stopifnot(is.matrix(sigma2), nrow(sigma2) == n1, ncol(sigma2) == n2)
    sd_mat <- sqrt(sigma2)
  } else {
    sd_scalar <- sqrt(sigma2)
  }
  message("estimating threshold:");
  for (b in 1:nsim) {
    Z <- matrix(stats::rnorm(n1 * n2), n1, n2)
    if (use_scalar_var) {
      Xnull <- Z * sd_scalar
    } else {
      Xnull <- Z * sd_mat
    }
    cs <- prefix_sum2d(Xnull)
    sums <- rect_sums_from_cs(cs, i1_vec, i2_vec, j1_vec, j2_vec)
    means <- sums / area
    maxima[b] <- max(abs(means))
    if (progress)  setTxtProgressBar(pb, b)
  }
  if (progress) close(pb)
  
  Q <- as.numeric(stats::quantile(maxima, probs = 1 - alpha_q, names = FALSE))
  list(
    Q = Q,
    maxima = maxima,
    alpha_block = alpha,
    alpha_quantile = alpha_q,
    L = c(L1, L2),
    m = c(m1, m2),
    blocks = data.frame(i1 = i1_vec, i2 = i2_vec, j1 = j1_vec, j2 = j2_vec, area = area)
  )
}


### simulation helper

simulate.fulldata <- function(N, spatial.var = 1) {
  matrix(rnorm(N * N, mean = 0, sd = sqrt(spatial.var)), nrow = N)
}

simulate.sar <- function(N, rho = 0.5, spatial.var = 1, tol = 1e-8, max_iter = 1000) {
  epsilon <- matrix(rnorm(N * N, mean = 0, sd = sqrt(spatial.var)), N, N)
  X <- epsilon
  
  # Create convolution kernel for averaging neighbors
  # This represents W operation without building the matrix
  
  for (iter in 1:max_iter) {
    X_old <- X
    
    # Compute neighbor averages using shifting
    neighbor_sum <- matrix(0, N, N)
    n_neighbors <- matrix(0, N, N)
    
    # Shift operations to get neighbor sums
    # Up neighbors
    neighbor_sum[2:N, ] <- neighbor_sum[2:N, ] + X[1:(N-1), ]
    n_neighbors[2:N, ] <- n_neighbors[2:N, ] + 1
    
    # Down neighbors  
    neighbor_sum[1:(N-1), ] <- neighbor_sum[1:(N-1), ] + X[2:N, ]
    n_neighbors[1:(N-1), ] <- n_neighbors[1:(N-1), ] + 1
    
    # Left neighbors
    neighbor_sum[, 2:N] <- neighbor_sum[, 2:N] + X[, 1:(N-1)]
    n_neighbors[, 2:N] <- n_neighbors[, 2:N] + 1
    
    # Right neighbors
    neighbor_sum[, 1:(N-1)] <- neighbor_sum[, 1:(N-1)] + X[, 2:N]
    n_neighbors[, 1:(N-1)] <- n_neighbors[, 1:(N-1)] + 1
    
    # Update: X = epsilon + rho * W * X
    n_neighbors[n_neighbors == 0] <- 1  # avoid division by zero
    X <- epsilon + rho * (neighbor_sum / n_neighbors)
    
    # Check convergence
    if (max(abs(X - X_old)) < tol) {
      break
    }
  }
  
  return(X)
}


# Multiple-patch version (backward compatible):
# - tau1, tau2: either length-2 numeric vectors (single patch) OR
#               two-column matrices/data.frames with one row per patch:
#               tau1[p,] = (row_start_frac, col_start_frac)
#               tau2[p,] = (row_end_frac,   col_end_frac)
# - delta: scalar or length = #patches (additive; overlapping patches sum)
epi_mean <- function(tau1, tau2, delta, N) {
  # normalize inputs to matrix form
  if (is.numeric(tau1) && length(tau1) == 2) {
    tau1 <- matrix(tau1, ncol = 2, byrow = TRUE)
  } else {
    tau1 <- as.matrix(tau1)
  }
  if (is.numeric(tau2) && length(tau2) == 2) {
    tau2 <- matrix(tau2, ncol = 2, byrow = TRUE)
  } else {
    tau2 <- as.matrix(tau2)
  }
  stopifnot(nrow(tau1) == nrow(tau2), ncol(tau1) == 2, ncol(tau2) == 2)
  
  P <- nrow(tau1)
  if (length(delta) == 1L) delta <- rep(delta, P)
  stopifnot(length(delta) == P)
  
  mean_mat <- matrix(0, nrow = N, ncol = N)
  
  for (p in seq_len(P)) {
    i1 <- max(1L, ceiling(N * tau1[p, 1]))
    i2 <- min(N,   floor(  N * tau2[p, 1]))
    j1 <- max(1L, ceiling(N * tau1[p, 2]))
    j2 <- min(N,   floor(  N * tau2[p, 2]))
    
    if (i1 <= i2 && j1 <= j2) {
      mean_mat[i1:i2, j1:j2] <- mean_mat[i1:i2, j1:j2] + delta[p]
    }
    # else: empty interval after rounding; skip silently
  }
  mean_mat
}
