###############################################
# Q threshold via simulation of block-mean maxima (3D)
# (matching Algorithm 2 block construction)
###############################################

source("multiple_cp_helper_3d.R")

simulate_Q_blockmeans_3d <- function(n, alpha, sigma2 = 1, nsim = 1000, seed = NULL,
                                     alpha_q = alpha, progress = TRUE) {
  stopifnot(length(n) == 3)
  n1 <- as.integer(n[1]); n2 <- as.integer(n[2]); n3 <- as.integer(n[3])

  L1 <- max(1L, floor(n1^alpha)); L2 <- max(1L, floor(n2^alpha)); L3 <- max(1L, floor(n3^alpha))
  m1 <- ceiling(n1 / L1); m2 <- ceiling(n2 / L2); m3 <- ceiling(n3 / L3)

  # Precompute all block rectangles and volumes
  nblocks <- m1 * m2 * m3
  i1_vec <- integer(nblocks); i2_vec <- integer(nblocks)
  j1_vec <- integer(nblocks); j2_vec <- integer(nblocks)
  k1_vec <- integer(nblocks); k2_vec <- integer(nblocks)
  vol    <- integer(nblocks)

  idx <- 0L
  for (s1 in 1:m1) for (s2 in 1:m2) for (s3 in 1:m3) {
    idx <- idx + 1L
    i1 <- (s1 - 1L) * L1 + 1L; i2 <- min(s1 * L1, n1)
    j1 <- (s2 - 1L) * L2 + 1L; j2 <- min(s2 * L2, n2)
    k1 <- (s3 - 1L) * L3 + 1L; k2 <- min(s3 * L3, n3)
    i1_vec[idx] <- i1; i2_vec[idx] <- i2
    j1_vec[idx] <- j1; j2_vec[idx] <- j2
    k1_vec[idx] <- k1; k2_vec[idx] <- k2
    vol[idx] <- (i2 - i1 + 1L) * (j2 - j1 + 1L) * (k2 - k1 + 1L)
  }

  if (!is.null(seed)) set.seed(seed)
  maxima <- numeric(nsim)
  if (progress) {
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  use_scalar_var <- (length(sigma2) == 1L)
  if (!use_scalar_var) {
    stopifnot(length(dim(sigma2)) == 3L,
              dim(sigma2)[1] == n1, dim(sigma2)[2] == n2, dim(sigma2)[3] == n3)
    sd_arr <- sqrt(sigma2)
  } else {
    sd_scalar <- sqrt(sigma2)
  }

  message("Estimating threshold Q (3D):")
  for (b in 1:nsim) {
    Z <- array(stats::rnorm(n1 * n2 * n3), dim = c(n1, n2, n3))
    if (use_scalar_var) {
      Xnull <- Z * sd_scalar
    } else {
      Xnull <- Z * sd_arr
    }
    cs <- prefix_sum3d(Xnull)
    sums <- box_sums_from_cs3d(cs, i1_vec, i2_vec, j1_vec, j2_vec, k1_vec, k2_vec)
    means <- sums / vol
    maxima[b] <- max(abs(means))
    if (progress) setTxtProgressBar(pb, b)
  }
  if (progress) close(pb)

  Q <- as.numeric(stats::quantile(maxima, probs = 1 - alpha_q, names = FALSE))
  list(
    Q = Q,
    maxima = maxima,
    alpha_block = alpha,
    alpha_quantile = alpha_q,
    L = c(L1, L2, L3),
    m = c(m1, m2, m3),
    blocks = data.frame(i1 = i1_vec, i2 = i2_vec,
                        j1 = j1_vec, j2 = j2_vec,
                        k1 = k1_vec, k2 = k2_vec,
                        vol = vol)
  )
}

###############################################
# Simulation helpers (3D)
###############################################

# SAR simulation for 3D field
simulate.sar.3d <- function(n, rho = 0.5, spatial.var = 1, tol = 1e-8, max_iter = 1000) {
  n1 <- n[1]; n2 <- n[2]; n3 <- n[3]
  epsilon <- array(rnorm(n1 * n2 * n3, mean = 0, sd = sqrt(spatial.var)), dim = c(n1, n2, n3))
  X <- epsilon

  for (iter in 1:max_iter) {
    X_old <- X
    neighbor_sum <- array(0, dim = c(n1, n2, n3))
    n_neighbors  <- array(0, dim = c(n1, n2, n3))

    # 6-face neighbors (rook in 3D)
    # dim 1: up/down
    if (n1 > 1) {
      neighbor_sum[2:n1, , ] <- neighbor_sum[2:n1, , ] + X[1:(n1-1), , ]
      n_neighbors[2:n1, , ] <- n_neighbors[2:n1, , ] + 1
      neighbor_sum[1:(n1-1), , ] <- neighbor_sum[1:(n1-1), , ] + X[2:n1, , ]
      n_neighbors[1:(n1-1), , ] <- n_neighbors[1:(n1-1), , ] + 1
    }
    # dim 2: left/right
    if (n2 > 1) {
      neighbor_sum[, 2:n2, ] <- neighbor_sum[, 2:n2, ] + X[, 1:(n2-1), ]
      n_neighbors[, 2:n2, ] <- n_neighbors[, 2:n2, ] + 1
      neighbor_sum[, 1:(n2-1), ] <- neighbor_sum[, 1:(n2-1), ] + X[, 2:n2, ]
      n_neighbors[, 1:(n2-1), ] <- n_neighbors[, 1:(n2-1), ] + 1
    }
    # dim 3: front/back
    if (n3 > 1) {
      neighbor_sum[, , 2:n3] <- neighbor_sum[, , 2:n3] + X[, , 1:(n3-1)]
      n_neighbors[, , 2:n3] <- n_neighbors[, , 2:n3] + 1
      neighbor_sum[, , 1:(n3-1)] <- neighbor_sum[, , 1:(n3-1)] + X[, , 2:n3]
      n_neighbors[, , 1:(n3-1)] <- n_neighbors[, , 1:(n3-1)] + 1
    }

    n_neighbors[n_neighbors == 0] <- 1
    X <- epsilon + rho * (neighbor_sum / n_neighbors)

    if (max(abs(X - X_old)) < tol) break
  }
  X
}

# Multiple-patch mean field (3D)
# tau1, tau2: K x 3 matrices with start/end fractions in [0,1]^3
# delta: vector of length K
epi_mean_3d <- function(tau1, tau2, delta, n) {
  if (is.numeric(tau1) && length(tau1) == 3) {
    tau1 <- matrix(tau1, ncol = 3, byrow = TRUE)
  } else { tau1 <- as.matrix(tau1) }
  if (is.numeric(tau2) && length(tau2) == 3) {
    tau2 <- matrix(tau2, ncol = 3, byrow = TRUE)
  } else { tau2 <- as.matrix(tau2) }
  stopifnot(nrow(tau1) == nrow(tau2), ncol(tau1) == 3, ncol(tau2) == 3)

  n1 <- n[1]; n2 <- n[2]; n3 <- n[3]
  P <- nrow(tau1)
  if (length(delta) == 1L) delta <- rep(delta, P)
  stopifnot(length(delta) == P)

  mean_arr <- array(0, dim = c(n1, n2, n3))
  for (p in seq_len(P)) {
    i1 <- max(1L, ceiling(n1 * tau1[p, 1])); i2 <- min(n1, floor(n1 * tau2[p, 1]))
    j1 <- max(1L, ceiling(n2 * tau1[p, 2])); j2 <- min(n2, floor(n2 * tau2[p, 2]))
    k1 <- max(1L, ceiling(n3 * tau1[p, 3])); k2 <- min(n3, floor(n3 * tau2[p, 3]))
    if (i1 <= i2 && j1 <= j2 && k1 <= k2) {
      mean_arr[i1:i2, j1:j2, k1:k2] <- mean_arr[i1:i2, j1:j2, k1:k2] + delta[p]
    }
  }
  mean_arr
}

# Label array from rectangular patches (for Rand Index)
boxes_to_label_3d <- function(tau1, tau2, n) {
  n1 <- n[1]; n2 <- n[2]; n3 <- n[3]
  lab <- array(0, dim = c(n1, n2, n3))
  K <- nrow(tau1)
  for (k in 1:K) {
    r1 <- pmax(1, floor(tau1[k, ] * c(n1, n2, n3)))
    r2 <- pmin(c(n1, n2, n3), ceiling(tau2[k, ] * c(n1, n2, n3)))
    lab[r1[1]:r2[1], r1[2]:r2[2], r1[3]:r2[3]] <- k
  }
  lab
}
