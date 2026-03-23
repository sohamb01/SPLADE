###############################################
# Utilities (3D, base R only)
###############################################

# 3D prefix sums (summed volume table)
# For an array A of dim (n1, n2, n3), computes cs such that
# cs[i,j,k] = sum of A[1:i, 1:j, 1:k].
prefix_sum3d <- function(A) {
  d <- dim(A)
  n1 <- d[1]; n2 <- d[2]; n3 <- d[3]
  cs <- A
  # cumsum along dim 1 (vectorized over dim2 x dim3)
  for (jk in 1:(n2 * n3)) {
    j <- ((jk - 1L) %% n2) + 1L
    k <- ((jk - 1L) %/% n2) + 1L
    cs[, j, k] <- cumsum(cs[, j, k])
  }
  # cumsum along dim 2 (vectorized over dim1 x dim3)
  for (ik in 1:(n1 * n3)) {
    i <- ((ik - 1L) %% n1) + 1L
    k <- ((ik - 1L) %/% n1) + 1L
    cs[i, , k] <- cumsum(cs[i, , k])
  }
  # cumsum along dim 3 (vectorized over dim1 x dim2)
  for (ij in 1:(n1 * n2)) {
    i <- ((ij - 1L) %% n1) + 1L
    j <- ((ij - 1L) %/% n1) + 1L
    cs[i, j, ] <- cumsum(cs[i, j, ])
  }
  cs
}

# Fast box sum from 3D prefix sums; s=(i1,j1,k1), t=(i2,j2,k2), inclusive corners
# Uses inclusion-exclusion with 8 terms.
box_sum3d <- function(cs, s, t) {
  i1 <- s[1]; j1 <- s[2]; k1 <- s[3]
  i2 <- t[1]; j2 <- t[2]; k2 <- t[3]
  # Ensure i1<=i2, j1<=j2, k1<=k2
  if (i1 > i2) { tmp <- i1; i1 <- i2; i2 <- tmp }
  if (j1 > j2) { tmp <- j1; j1 <- j2; j2 <- tmp }
  if (k1 > k2) { tmp <- k1; k1 <- k2; k2 <- tmp }

  # Helper: cs value, returning 0 if any index is < 1
  .cs <- function(a, b, c) {
    if (a < 1 || b < 1 || c < 1) return(0)
    cs[a, b, c]
  }

  # Inclusion-exclusion
  .cs(i2, j2, k2) -
    .cs(i1 - 1L, j2, k2) - .cs(i2, j1 - 1L, k2) - .cs(i2, j2, k1 - 1L) +
    .cs(i1 - 1L, j1 - 1L, k2) + .cs(i1 - 1L, j2, k1 - 1L) + .cs(i2, j1 - 1L, k1 - 1L) -
    .cs(i1 - 1L, j1 - 1L, k1 - 1L)
}

# Vectorized box sums from 3D prefix sums (for quantile simulation)
box_sums_from_cs3d <- function(cs, i1, i2, j1, j2, k1, k2) {
  n <- length(i1)
  ii1 <- pmin(i1, i2); ii2 <- pmax(i1, i2)
  jj1 <- pmin(j1, j2); jj2 <- pmax(j1, j2)
  kk1 <- pmin(k1, k2); kk2 <- pmax(k1, k2)

  dims <- dim(cs)
  n1 <- dims[1]; n2 <- dims[2]; n3 <- dims[3]

  # Helper: vectorized cs lookup, returning 0 where any index < 1
  .csv <- function(a, b, c) {
    valid <- (a >= 1) & (b >= 1) & (c >= 1)
    out <- numeric(n)
    if (any(valid)) {
      # Convert to linear index: (a-1) + (b-1)*n1 + (c-1)*n1*n2 + 1
      idx <- a[valid] + (b[valid] - 1L) * n1 + (c[valid] - 1L) * n1 * n2
      out[valid] <- cs[idx]
    }
    out
  }

  .csv(ii2, jj2, kk2) -
    .csv(ii1 - 1L, jj2, kk2) - .csv(ii2, jj1 - 1L, kk2) - .csv(ii2, jj2, kk1 - 1L) +
    .csv(ii1 - 1L, jj1 - 1L, kk2) + .csv(ii1 - 1L, jj2, kk1 - 1L) + .csv(ii2, jj1 - 1L, kk1 - 1L) -
    .csv(ii1 - 1L, jj1 - 1L, kk1 - 1L)
}

# Connected-components on a logical 3D mask. connectivity = 6, 18, or 26.
label_components_3d <- function(mask, connectivity = 26L) {
  stopifnot(length(dim(mask)) == 3L)
  d <- dim(mask)
  n1 <- d[1]; n2 <- d[2]; n3 <- d[3]
  lab <- array(0L, dim = d)
  visited <- array(FALSE, dim = d)

  # Build neighbor offsets
  if (connectivity == 6L) {
    nbrs <- rbind(c(1,0,0), c(-1,0,0), c(0,1,0), c(0,-1,0), c(0,0,1), c(0,0,-1))
  } else if (connectivity == 18L) {
    nbrs <- as.matrix(expand.grid(dx = -1:1, dy = -1:1, dz = -1:1))
    nbrs <- nbrs[rowSums(abs(nbrs)) <= 2 & rowSums(abs(nbrs)) > 0, , drop = FALSE]
  } else {
    # 26-connectivity
    nbrs <- as.matrix(expand.grid(dx = -1:1, dy = -1:1, dz = -1:1))
    nbrs <- nbrs[!(nbrs[,1] == 0 & nbrs[,2] == 0 & nbrs[,3] == 0), , drop = FALSE]
  }

  comp_sizes <- integer()
  id <- 0L
  max_queue <- min(n1 * n2 * n3, 5e6)  # cap queue size for memory

  for (i in 1:n1) for (j in 1:n2) for (k in 1:n3) {
    if (mask[i, j, k] && !visited[i, j, k]) {
      id <- id + 1L
      # BFS
      qi <- integer(max_queue); qj <- integer(max_queue); qk <- integer(max_queue)
      head <- 1L; tail <- 1L
      qi[tail] <- i; qj[tail] <- j; qk[tail] <- k; tail <- tail + 1L
      visited[i, j, k] <- TRUE; lab[i, j, k] <- id
      size <- 0L
      while (head < tail) {
        ci <- qi[head]; cj <- qj[head]; ck <- qk[head]; head <- head + 1L
        size <- size + 1L
        for (nb in 1:nrow(nbrs)) {
          ni <- ci + nbrs[nb, 1]; nj <- cj + nbrs[nb, 2]; nk <- ck + nbrs[nb, 3]
          if (ni >= 1 && ni <= n1 && nj >= 1 && nj <= n2 && nk >= 1 && nk <= n3) {
            if (mask[ni, nj, nk] && !visited[ni, nj, nk]) {
              visited[ni, nj, nk] <- TRUE
              lab[ni, nj, nk] <- id
              qi[tail] <- ni; qj[tail] <- nj; qk[tail] <- nk; tail <- tail + 1L
            }
          }
        }
      }
      comp_sizes[id] <- size
    }
  }
  list(labels = lab, sizes = comp_sizes)
}

############################################
# Algorithm 1 (3D): Single spatial change-point
############################################
single_sp_changepoint_3d <- function(X, alpha, kappa,  C = 1, connectivity = 26L) {
  stopifnot(length(dim(X)) == 3L)
  d <- dim(X)
  n1 <- d[1]; n2 <- d[2]; n3 <- d[3]
  n_total <- n1 * n2 * n3

  # Step 1: L_k, m_k, Y_k (sub-sample starts)
  L1 <- max(1L, floor(n1^alpha)); L2 <- max(1L, floor(n2^alpha)); L3 <- max(1L, floor(n3^alpha))
  m1 <- ceiling(n1 / L1); m2 <- ceiling(n2 / L2); m3 <- ceiling(n3 / L3)
  Y1 <- ((0:(m1 - 1)) * L1) + 1L
  Y2 <- ((0:(m2 - 1)) * L2) + 1L
  Y3 <- ((0:(m3 - 1)) * L3) + 1L

  # Step 2: Sub-sampled dataset
  Y <- X[Y1, Y2, Y3, drop = FALSE]
  m <- length(Y)

  # Step 3: preliminary naive estimate on coarse grid
  csY <- prefix_sum3d(Y)
  best_val <- -Inf
  a_I <- c(1L, 1L, 1L); b_I <- c(m1, m2, m3)
  m_total <- m1 * m2 * m3
  sumAll <- csY[m1, m2, m3]

  for (i1 in 1:m1) for (j1 in 1:m2) for (k1 in 1:m3) {
    for (i2 in i1:m1) for (j2 in j1:m2) for (k2 in k1:m3) {
      vol <- (i2 - i1 + 1L) * (j2 - j1 + 1L) * (k2 - k1 + 1L)
      if (vol <= (m_total * 0.2) || vol == m_total) next
      sumR <- box_sum3d(csY, c(i1, j1, k1), c(i2, j2, k2))
      meanR <- sumR / vol
      meanC <- (sumAll - sumR) / (m_total - vol)
      val <- sqrt((vol / m_total) * (1 - vol / m_total)) * abs(meanR - meanC)
      if (val > best_val) {
        best_val <- val
        a_I <- c(i1, j1, k1); b_I <- c(i2, j2, k2)
      }
    }
  }

  # Step 4: define bands around L * a_I and L * b_I
  Lvec <- c(L1, L2, L3)
  band_len <- ceiling(C * Lvec * 0.5 * (min(n1, n2, n3))^(kappa) * (log(n_total))^(1/3))
  centerL <- Lvec * a_I
  centerR <- Lvec * b_I
  nvec <- c(n1, n2, n3)

  L_B <- list(); R_B <- list()
  for (dd in 1:3) {
    L_B[[dd]] <- max(1L, centerL[dd] - band_len[dd]) : min(nvec[dd], centerL[dd] + band_len[dd])
    R_B[[dd]] <- max(1L, centerR[dd] - band_len[dd]) : min(nvec[dd], centerR[dd] + band_len[dd])
  }

  # Step 5: refined argmax over restricted search space
  csX <- prefix_sum3d(X)
  sumAllX <- csX[n1, n2, n3]
  best_val <- -Inf
  s_star <- c(1L, 1L, 1L); t_star <- c(n1, n2, n3)
  for (i1 in L_B[[1]]) for (j1 in L_B[[2]]) for (k1 in L_B[[3]]) {
    for (i2 in R_B[[1]]) for (j2 in R_B[[2]]) for (k2 in R_B[[3]]) {
      if (i2 < i1 || j2 < j1 || k2 < k1) next
      vol <- (i2 - i1 + 1L) * (j2 - j1 + 1L) * (k2 - k1 + 1L)
      if (vol == 0L || vol == n_total) next
      sumR <- box_sum3d(csX, c(i1, j1, k1), c(i2, j2, k2))
      meanR <- sumR / vol
      meanC <- (sumAllX - sumR) / (n_total - vol)
      val <- sqrt((vol / n_total) * (1 - vol / n_total)) * abs(meanR - meanC)
      if (val > best_val) {
        best_val <- val
        s_star <- c(i1, j1, k1); t_star <- c(i2, j2, k2)
      }
    }
  }

  list(
    tilde_I = list(s = s_star, t = t_star),
    L = Lvec, m = c(m1, m2, m3),
    Y_idx = list(Y1 = Y1, Y2 = Y2, Y3 = Y3),
    Y = Y, a_I = a_I, b_I = b_I,
    bands = list(L_B = L_B, R_B = R_B)
  )
}

#########################################################
# Algorithm 2 (3D): Multiple spatial change-point estimation
#########################################################
multiple_sp_changepoints_3d <- function(
    X, alpha, kappa, Q = NULL, c = 1, connectivity = 26L,
    C_for_algo1 = 1
) {
  stopifnot(length(dim(X)) == 3L)
  d <- dim(X)
  n1 <- d[1]; n2 <- d[2]; n3 <- d[3]
  n_total <- n1 * n2 * n3

  L1 <- max(1L, floor(n1^alpha)); L2 <- max(1L, floor(n2^alpha)); L3 <- max(1L, floor(n3^alpha))
  m1 <- ceiling(n1 / L1); m2 <- ceiling(n2 / L2); m3 <- ceiling(n3 / L3)

  I_fn <- function(s, Lk, nk) {
    start <- (s - 1L) * Lk + 1L; end <- min(s * Lk, nk); start:end
  }

  # Compute block means using prefix sums
  csX <- prefix_sum3d(X)
  block_mean <- array(NA_real_, dim = c(m1, m2, m3))
  for (s1 in 1:m1) for (s2 in 1:m2) for (s3 in 1:m3) {
    i_rng <- I_fn(s1, L1, n1); j_rng <- I_fn(s2, L2, n2); k_rng <- I_fn(s3, L3, n3)
    sumB <- box_sum3d(csX, c(min(i_rng), min(j_rng), min(k_rng)),
                            c(max(i_rng), max(j_rng), max(k_rng)))
    block_mean[s1, s2, s3] <- sumB / (length(i_rng) * length(j_rng) * length(k_rng))
  }

  if (is.null(Q)) {
    Q <- as.numeric(stats::quantile(as.numeric(block_mean), probs = 1 - alpha, na.rm = TRUE))
  }

  # Threshold: flag blocks
  Mmask <- array(FALSE, dim = c(n1, n2, n3))
  for (s1 in 1:m1) for (s2 in 1:m2) for (s3 in 1:m3) {
    if (abs(block_mean[s1, s2, s3]) > Q) {
      i_rng <- I_fn(s1, L1, n1); j_rng <- I_fn(s2, L2, n2); k_rng <- I_fn(s3, L3, n3)
      Mmask[i_rng, j_rng, k_rng] <- TRUE
    }
  }

  # Connected components
  cc <- label_components_3d(Mmask, connectivity = connectivity)
  labs <- cc$labels; sizes <- cc$sizes
  keep_ids <- which(sizes > c * (n_total^alpha))

  I_hat <- list()
  I_bounds <- matrix(0L, nrow = 0, ncol = 6,
                     dimnames = list(NULL, c("i1","j1","k1","i2","j2","k2")))
  D_bounds <- matrix(0L, nrow = 0, ncol = 6,
                     dimnames = list(NULL, c("i1","j1","k1","i2","j2","k2")))

  message(paste0("First stage done! ", length(keep_ids), " change-points detected"))

  if (length(keep_ids) == 0) {
    return(list(K_hat = 0, I_hat = NA, I_bounds = NA, D_bounds = NA, first_stage_data = NA))
  }

  n_ids <- length(keep_ids)
  pb <- utils::txtProgressBar(min = 0, max = n_ids, style = 3)
  on.exit(close(pb), add = TRUE)
  message("Second stage:")

  Lvec <- c(L1, L2, L3)
  nvec <- c(n1, n2, n3)

  for (j in seq_along(keep_ids)) {
    id <- keep_ids[j]
    pts <- which(labs == id, arr.ind = TRUE)  # Nx3 matrix

    # Map voxels back to block indices
    s1s <- pmin(ceiling(pts[, 1] / L1), m1)
    s2s <- pmin(ceiling(pts[, 2] / L2), m2)
    s3s <- pmin(ceiling(pts[, 3] / L3), m3) 
    l_block <- c(min(s1s), min(s2s), min(s3s))
    r_block <- c(max(s1s), max(s2s), max(s3s))

    # Expand by log factor
    expand <- ceiling(c * Lvec * 0.5 * (log(n_total))^(1/3))
    l_vec <- Lvec * l_block - expand
    r_vec <- Lvec * r_block + expand

    lo <- pmax(1L, floor(l_vec))
    hi <- pmin(nvec, ceiling(r_vec))

    D_bounds <- rbind(D_bounds, c(i1 = lo[1], j1 = lo[2], k1 = lo[3],
                                  i2 = hi[1], j2 = hi[2], k2 = hi[3]))

    D_sub <- X[lo[1]:hi[1], lo[2]:hi[2], lo[3]:hi[3], drop = FALSE]
    a1 <- single_sp_changepoint_3d(D_sub, alpha = alpha, kappa=kappa, C = C_for_algo1,
                                   connectivity = connectivity)
    s_loc <- a1$tilde_I$s; t_loc <- a1$tilde_I$t
    s_glob <- (lo - 1L) + s_loc
    t_glob <- (lo - 1L) + t_loc
    I_hat[[j]] <- list(s = s_glob, t = t_glob)
    I_bounds <- rbind(I_bounds,
                      c(i1 = s_glob[1], j1 = s_glob[2], k1 = s_glob[3],
                        i2 = t_glob[1], j2 = t_glob[2], k2 = t_glob[3]))

    utils::setTxtProgressBar(pb, j)
  }

  list(K_hat = length(I_hat),
       I_hat = I_hat,
       I_bounds = I_bounds,
       D_bounds = D_bounds,
       first_stage_data = list(labs, keep_ids))
}
