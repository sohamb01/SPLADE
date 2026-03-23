#############################
# Utilities (2D, base R only)
#############################

# 2D prefix sums (integral image)
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

# Plot a matrix with matrix row/col coordinates preserved
plot_matrix <- function(M, main = "", col = gray.colors(256), legend = FALSE, zlim = NULL) {
  nr <- nrow(M); nc <- ncol(M)
  if (is.null(zlim)) zlim <- range(M, finite = TRUE)
  image(x = 1:nc, y = 1:nr, z = t(M[nr:1, , drop = FALSE]),
        xlab = "", ylab = "", axes = FALSE, useRaster = TRUE, col = col, zlim = zlim, main = main)
  box()
  if (legend) {
    # simple legend strip
    z <- matrix(seq(zlim[1], zlim[2], length.out = 100), nrow = 1)
    par(new = TRUE, fig = c(0.9, 0.93, 0.15, 0.85), mar = c(1,1,1,1))
    image(z = z, axes = FALSE, col = col, xlab = "", ylab = "", useRaster = TRUE)
    par(new = FALSE)
  }
}

# Draw rectangle using matrix indices (rows i1:i2, cols j1:j2)
draw_rect_indices <- function(i1, i2, j1, j2, nrow_total, col = "red", lwd = 2, lty = 1) {
  # mapping to image() coords used in plot_matrix
  xleft   <- j1 - 0.5
  xright  <- j2 + 0.5
  ybottom <- (nrow_total - i2) + 0.5
  ytop    <- (nrow_total - i1) + 0.5
  rect(xleft, ybottom, xright, ytop, border = col, lwd = lwd, lty = lty)
}

# Connected-components on a logical mask (2D). connectivity = 4 ("rook") or 8 ("queen")
label_components <- function(mask, connectivity = 4L) {
  stopifnot(is.matrix(mask))
  nr <- nrow(mask); nc <- ncol(mask)
  lab <- matrix(0L, nr, nc)
  visited <- matrix(FALSE, nr, nc)
  if (connectivity == 4L) {
    nbrs <- rbind(c( 1, 0), c(-1, 0), c(0, 1), c(0, -1))
  } else {
    nbrs <- rbind(expand.grid(dx = -1:1, dy = -1:1))
    nbrs <- as.matrix(nbrs[!(nbrs[,1] == 0 & nbrs[,2] == 0), ])
  }
  comp_sizes <- integer()
  id <- 0L
  for (i in 1:nr) for (j in 1:nc) {
    if (mask[i, j] && !visited[i, j]) {
      id <- id + 1L
      # BFS queue
      qi <- integer(nr * nc); qj <- integer(nr * nc); head <- 1L; tail <- 1L
      qi[tail] <- i; qj[tail] <- j; tail <- tail + 1L
      visited[i, j] <- TRUE; lab[i, j] <- id
      size <- 0L
      while (head < tail) {
        ci <- qi[head]; cj <- qj[head]; head <- head + 1L
        size <- size + 1L
        for (k in 1:nrow(nbrs)) {
          ni <- ci + nbrs[k, 1]; nj <- cj + nbrs[k, 2]
          if (ni >= 1 && ni <= nr && nj >= 1 && nj <= nc) {
            if (mask[ni, nj] && !visited[ni, nj]) {
              visited[ni, nj] <- TRUE
              lab[ni, nj] <- id
              qi[tail] <- ni; qj[tail] <- nj; tail <- tail + 1L
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
# Algorithm 1: Single spatial change-point
############################################
# Implements your steps exactly (on a 2D array X)
# Returns: list with tilde_I = list(s = c(i1,j1), t = c(i2,j2)), plus internals
single_sp_changepoint <- function(X, alpha, kappa, C = 1, connectivity = 4L) {
  stopifnot(length(dim(X)) == 2L)
  n1 <- nrow(X); n2 <- ncol(X)
  n_total <- n1 * n2
  
  # Step 1: L_k, m_k, Y_k (starts)
  L1 <- max(1L, floor(n1^alpha)); L2 <- max(1L, floor(n2^alpha))
  m1 <- ceiling(n1 / L1);        m2 <- ceiling(n2 / L2)
  Y1 <- ((0:(m1 - 1)) * L1) + 1L
  Y2 <- ((0:(m2 - 1)) * L2) + 1L
  
  # Step 2: Sub-sampled dataset 𝒀
  Y <- X[Y1, Y2, drop = FALSE]
  m <- length(Y)
  
  # Step 3: preliminary naive estimate based on 𝒀
  # We "naively" scan rectangles on the coarse grid (m1 x m2) to get a_I, b_I.
  # (This stays faithful: it's an estimate based on 𝒀; exact form wasn’t specified.)
  csY <- prefix_sum2d(Y)
  best_val <- -Inf
  a_I <- c(1L,1L); b_I <- c(m1, m2)
  for (i1 in 1:m1) for (j1 in 1:m2) {
    for (i2 in i1:m1) for (j2 in j1:m2) {
      area <- (i2 - i1 + 1L) * (j2 - j1 + 1L)
      if (area <= (m1*m2*0.2) || area == (m1 * m2)) next
      sumR <- rect_sum2d(csY, c(i1, j1), c(i2, j2))
      meanR <- sumR / area
      sumAll <- csY[m1, m2]
      meanC <- (sumAll - sumR) / (m1 * m2 - area)
      val <- sqrt( (area/ (m1 * m2)) * (1 - (area/ (m1 * m2))) ) * abs(meanR - meanC)
      if (val > best_val) {
        best_val <- val
        a_I <- c(i1, j1); b_I <- c(i2, j2)
      }
    }
  }
  
  # Step 4: define bands \hat L_B and \hat R_B around L * a_I and L * b_I
  band_len1 <- ceiling(C * L1 * 0.5 * (min(n1, n2))^(kappa) * (log(n_total))^(1/2) )
  band_len2 <- ceiling(C * L2 * 0.5 * (log(n_total))^(1/2) )
  centerL  <- c(L1 * a_I[1], L2 * a_I[2])
  centerR  <- c(L1 * b_I[1], L2 * b_I[2])
  
  L_B1 <- max(1L, centerL[1] - band_len1) : min(n1, centerL[1] + band_len1)
  L_B2 <- max(1L, centerL[2] - band_len2) : min(n2, centerL[2] + band_len2)
  R_B1 <- max(1L, centerR[1] - band_len1) : min(n1, centerR[1] + band_len1)
  R_B2 <- max(1L, centerR[2] - band_len2) : min(n2, centerR[2] + band_len2)
  
  # Step 5: refined argmax over (\bb{i} in \hat L_B, \bb{j} in \hat R_B)
  csX <- prefix_sum2d(X)
  sumAllX <- csX[n1, n2]
  best_val <- -Inf
  s_star <- c(1L,1L); t_star <- c(n1, n2)
  for (i1 in L_B1) for (j1 in L_B2) {
    for (i2 in R_B1) for (j2 in R_B2) {
      if (i2 < i1 || j2 < j1) next  # need a valid rectangle
      area <- (i2 - i1 + 1L) * (j2 - j1 + 1L)
      if (area == 0L || area == n_total) next
      sumR <- rect_sum2d(csX, c(i1, j1), c(i2, j2))
      meanR <- sumR / area
      meanC <- (sumAllX - sumR) / (n_total - area)
      val <- sqrt((area/n_total) * (1- area / n_total)) * abs(meanR - meanC)
      if (val > best_val) {
        paste0("completed", i1, " ", i2, " ", j1, " ", j2)
        best_val <- val
        s_star <- c(i1, j1); t_star <- c(i2, j2)
      }
    }
  }
  
  list(
    tilde_I = list(s = s_star, t = t_star),
    L = c(L1, L2), m = c(m1, m2),
    Y_idx = list(Y1 = Y1, Y2 = Y2),
    Y = Y, a_I = a_I, b_I = b_I,
    bands = list(L_B1 = L_B1, L_B2 = L_B2, R_B1 = R_B1, R_B2 = R_B2)
  )
}

#########################################################
# Algorithm 2: Multiple spatial change-point estimation
#########################################################
# Returns rectangles \hat I_j in global coords, and draws:
#   (1) selected connected components
#   (2) final estimated rectangles overlaid on X
multiple_sp_changepoints <- function(
    X, alpha, kappa, Q = NULL, c = 1, connectivity = 8L,
    C_for_algo1 = 1, plot_components = TRUE, plot_estimates = TRUE
) {
  stopifnot(length(dim(X)) == 2L)
  n1 <- nrow(X); n2 <- ncol(X)
  n_total <- n1 * n2
  
  L1 <- max(1L, floor(n1^alpha)); L2 <- max(1L, floor(n2^alpha))
  m1 <- ceiling(n1 / L1);        m2 <- ceiling(n2 / L2)
  I1 <- function(s) { start <- (s - 1L) * L1 + 1L; end <- min(s * L1, n1); start:end }
  I2 <- function(s) { start <- (s - 1L) * L2 + 1L; end <- min(s * L2, n2); start:end }
  
  csX <- prefix_sum2d(X)
  block_mean <- matrix(NA_real_, m1, m2)
  for (s1 in 1:m1) for (s2 in 1:m2) {
    i_rng <- I1(s1); j_rng <- I2(s2)
    sumB <- rect_sum2d(csX, c(min(i_rng), min(j_rng)), c(max(i_rng), max(j_rng)))
    block_mean[s1, s2] <- sumB / (length(i_rng) * length(j_rng))
  }
  if (is.null(Q)) {
    Q <- as.numeric(stats::quantile(as.numeric(block_mean), probs = 1 - alpha, na.rm = TRUE))
  }
  
  Mmask <- matrix(FALSE, n1, n2)
  for (s1 in 1:m1) for (s2 in 1:m2) {
    if (abs(block_mean[s1, s2]) > Q) Mmask[I1(s1), I2(s2)] <- TRUE
  }
  
  cc <- label_components(Mmask, connectivity = connectivity)
  labs <- cc$labels; sizes <- cc$sizes
  keep_ids <- which(sizes > c * (n_total^alpha))
  
  if (plot_components) {
    show <- matrix(0L, n1, n2)
    for (id in keep_ids) show[labs == id] <- match(id, keep_ids)
    pal <- c("#FFFFFF", grDevices::hcl.colors(max(1, length(keep_ids)), "Set3"))
    plot_matrix(show, main = "Selected connected components", col = pal, zlim = c(0, max(show)))
  }
  
  I_hat <- list()
  I_bounds <- matrix(0L, nrow = 0, ncol = 4, dimnames = list(NULL, c("i1","j1","i2","j2")))
  D_bounds <- matrix(0L, nrow = 0, ncol = 4, dimnames = list(NULL, c("i1","j1","i2","j2")))  # <-- NEW
  
  print(paste0("First stage done! ", length(keep_ids), " change-points detected"))
  
  if(length(keep_ids)==0){
    return(  list(
           K_hat = 0,
           I_hat = NA,
           I_bounds = NA,     # final estimate corners
           D_bounds = NA,     # <-- NEW: search-window corners
           first_stage_data = NA)
    )
  }
  else{
  n_ids <- length(keep_ids)
  pb <- utils::txtProgressBar(min = 0, max = n_ids, style = 3)
  on.exit(close(pb), add = TRUE)
  utils::setTxtProgressBar(pb, 0)
  message("second stage:")
  for (j in seq_along(keep_ids)) {
    id <- keep_ids[j]
    pts <- which(labs == id, arr.ind = TRUE)
    s1s <- pmin(ceiling(pts[, 1] / L1), m1)
    s2s <- pmin(ceiling(pts[, 2] / L2), m2)
    l1_j <- min(s1s); r1_j <- max(s1s)
    l2_j <- min(s2s); r2_j <- max(s2s)
    
    expand1 <- ceiling(c * L1 *  0.5 * (log(n_total))^(1/2) )
    expand2 <- ceiling(c * L2 *  0.5 * (log(n_total))^(1/2) )
    l_vec <- c(L1 * l1_j - expand1, L2 * l2_j - expand2)
    r_vec <- c(L1 * r1_j + expand1, L2 * r2_j + expand2)
    i1 <- max(1L, floor(l_vec[1])); j1 <- max(1L, floor(l_vec[2]))
    i2 <- min(n1, ceiling(r_vec[1])); j2 <- min(n2, ceiling(r_vec[2]))
    
    # Save the search window bounds (D_j)
    D_bounds <- rbind(D_bounds, c(i1 = i1, j1 = j1, i2 = i2, j2 = j2))
    
    D_sub <- X[i1:i2, j1:j2, drop = FALSE]
    a1 <- single_sp_changepoint(D_sub, alpha = alpha, kappa=kappa, C = C_for_algo1, connectivity = connectivity)
    s_loc <- a1$tilde_I$s; t_loc <- a1$tilde_I$t
    s_glob <- c(i1 - 1L, j1 - 1L) + s_loc
    t_glob <- c(i1 - 1L, j1 - 1L) + t_loc
    I_hat[[j]] <- list(s = s_glob, t = t_glob)
    I_bounds <- rbind(I_bounds,
                      c(i1 = s_glob[1], j1 = s_glob[2],
                        i2 = t_glob[1], j2 = t_glob[2]))
    
    
    utils::setTxtProgressBar(pb, j)
  }
  
  if (plot_estimates) {
    plot_matrix(X, main = "Final estimates (rectangles)", col = gray.colors(256))
    for (r in I_hat) {
      s <- r$s; t <- r$t
      draw_rect_indices(i1 = s[1], i2 = t[1], j1 = s[2], j2 = t[2], nrow_total = nrow(X),
                        col = "red", lwd = 2, lty = 1)
    }
  }
  
  list(K_hat = length(I_hat),
       I_hat = I_hat,
       I_bounds = I_bounds,     # final estimate corners
       D_bounds = D_bounds,     # <-- NEW: search-window corners
       first_stage_data = list(labs, keep_ids))
  }
}

