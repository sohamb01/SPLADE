border_band<- function(n1, n2,
                       w_row = ceiling(sqrt(n1)),   # thickness in rows (top/bottom)
                       w_col = ceiling(sqrt(n2))) { # thickness in cols (left/right)
  stopifnot(n1 >= 1L, n2 >= 1L)
  
  # Clip widths so the inner rectangle isn't inverted
  w_row <- max(0L, min(w_row, floor(n1 / 2)))
  w_col <- max(0L, min(w_col, floor(n2 / 2)))
  
  # Build full grid
  i <- rep.int(seq_len(n1), times = n2)
  j <- rep(seq_len(n2), each = n1)
  
  # Points inside the *inner* unshaded rectangle
  inner <- (i > w_row) & (i <= n1 - w_row) &
    (j > w_col) & (j <= n2 - w_col)
  
  # Keep the complement = shaded border band
  keep <- !inner
  cbind(i = i[keep], j = j[keep])
}

rects_to_label <- function(tau1, tau2, N) {
  lab <- matrix(0, nrow = N, ncol = N)
  K <- nrow(tau1)
  for (k in 1:K) {
    r1 <- pmax(1, floor(tau1[k, ] * N))
    r2 <- pmin(N, ceiling(tau2[k, ] * N))
    lab[r1[1]:r2[1], r1[2]:r2[2]] <- k
  }
  lab
}



lrv_estimator_outer_shell <- function(X, alpha, Bn_pow = 1/3, kernel = c("epanechnikov"),
                                      show_progress = TRUE) {
  kernel <- match.arg(kernel)
  stopifnot(is.matrix(X), is.numeric(X))
  n1 <- nrow(X); n2 <- ncol(X)
  n_tot <- n1 * n2
  #Bn <- c((n1)^(Bn_pow),(n2)^(Bn_pow))   # keep configurable per your note
  
  outer_shell <- border_band(n1, n2)
  nloc <- nrow(outer_shell)
  
  Xbar <- mean(X[cbind(outer_shell[, "i"], outer_shell[, "j"])])
  
  
  Bn <- c((nloc)^(Bn_pow/2),(nloc)^(Bn_pow/2)) 
  
  Y <- matrix(0, n1, n2) 
  
  Y[cbind(outer_shell[, "i"], outer_shell[, "j"])]<-
    X[cbind(outer_shell[, "i"], outer_shell[, "j"])] - Xbar
  
  
  #Y <- X - Xbar  # centered residual field
  
  
  
  # ------------------------------
  # 2) Kernel-weighted lag sum
  #    Bounded support: Epanechnikov has omega = 1
  #    => |h_k| <= floor(Bn) in index units
  # ------------------------------
  H1 <- min(floor(Bn), n1 - 1L)
  H2 <- min(floor(Bn), n2 - 1L)
  h1_vals <- (-H1):H1
  h2_vals <- (-H2):H2
  
  # precompute 1D kernel weights
  if (kernel == "epanechnikov") {
    k1 <- K_epanechnikov(h1_vals / Bn[1])
    k2 <- K_epanechnikov(h2_vals / Bn[2])
  } else {
    stop("Only 'epanechnikov' is implemented.")
  }
  
  acc <- 0.0
  total_steps <- length(h1_vals) * length(h2_vals)
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = total_steps, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  step <- 0L
  
  for (a in seq_along(h1_vals)) {
    h1 <- h1_vals[a]
    # valid row indices i with i and i+h1 in [1, n1]
    if (h1 >= 0L) {
      rows <- 1:(n1 - h1)
    } else {
      rows <- (1 - h1):n1
    }
    if (length(rows) == 0L) {
      if (show_progress) { step <- step + length(h2_vals); utils::setTxtProgressBar(pb, step) }
      next
    }
    
    for (b in seq_along(h2_vals)) {
      h2 <- h2_vals[b]
      w <- k1[a] * k2[b]
      step <- step + 1L
      if (w == 0) {
        if (show_progress) utils::setTxtProgressBar(pb, step)
        next
      }
      
      # valid col indices j with j and j+h2 in [1, n2]
      if (h2 >= 0L) {
        cols <- 1:(n2 - h2)
      } else {
        cols <- (1 - h2):n2
      }
      if (length(cols) == 0L) {
        if (show_progress) utils::setTxtProgressBar(pb, step)
        next
      }
      
      # elementwise product of overlapping submatrices
      A <- Y[rows, cols, drop = FALSE]
      B <- Y[rows + h1, cols + h2, drop = FALSE]
      acc <- acc + w * sum(A * B)
      
      if (show_progress) utils::setTxtProgressBar(pb, step)
    }
  }
  
  Sigma_hat <- as.numeric(acc / nloc)
  Sigma_hat
}

source("multiple_cp_helper.R")
source("LRV_helper.R")
source("quantile_calc_simu.R")


norm_rect <- function(bl, tr) {
  # Ensure proper ordering and clip to [0,1]^2
  xmin <- max(0, min(bl[1], tr[1])); xmax <- min(1, max(bl[1], tr[1]))
  ymin <- max(0, min(bl[2], tr[2])); ymax <- min(1, max(bl[2], tr[2]))
  c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}
area_rect <- function(r) {
  w <- max(0, r["xmax"] - r["xmin"]); h <- max(0, r["ymax"] - r["ymin"]); w * h
}
inter_area <- function(a, b) {
  w <- max(0, min(a["xmax"], b["xmax"]) - max(a["xmin"], b["xmin"]))
  h <- max(0, min(a["ymax"], b["ymax"]) - max(a["ymin"], b["ymin"]))
  w * h
}

rect_symdiff_area_bltr <- function(bl1, tr1, bl2, tr2) {
  stopifnot(length(bl1) == 2, length(tr1) == 2, length(bl2) == 2, length(tr2) == 2)
  
  
  
  r1 <- norm_rect(bl1, tr1)
  r2 <- norm_rect(bl2, tr2)
  
  a1 <- area_rect(r1)
  a2 <- area_rect(r2)
  ai <- inter_area(r1, r2)
  
  # Symmetric difference area: |A| + |B| - 2|A ∩ B|
  a1 + a2 - 2 * ai
}


spatial.var <- 1

library(parallel)
library(foreach)
num_core=5

totaliter=100



# tau1 <- rbind(c(0.20, 0.20), c(0.60, 0.20) , c(0.65, 0.65))
# tau2 <- rbind(c(0.45, 0.60) , c(0.80, 0.45), c(0.85, 0.85))

tau1 <- rbind(c(0.20, 0.20), c(0.60, 0.60), c(0.65, 0.15))
tau2 <- rbind(c(0.45, 0.70), c(0.85, 0.85), c(0.85, 0.45))

#tau1 <- tau1[order(tau1[,1], tau1[,2]), ]
#tau2 <- tau2[order(tau1[,1], tau1[,2]), ]



datfin <- data.frame()

for (N in c(250, 500)) {
  print("Gridzise"); print(N)
  true_lab <- rects_to_label(tau1, tau2, N)
  Q_base <- simulate_Q_blockmeans(n = c(N, N), alpha = 0.5, alpha_q=0.5, 
                        sigma2 = 1, nsim = 1000, seed = 123)
  for (jsize in  c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) {
    delta <- c(jsize, jsize, -jsize)
    print("Jumpsize"); print(jsize)
    mu <- epi_mean(tau1, tau2, delta, N)
    for (rho in c(0,0.25, 0.5,0.75)) {
      print("Spatial Correlation"); print(rho)
      
      start <- Sys.time()
      
      cl <- parallel::makeCluster(num_core, type = "PSOCK")
      doSNOW::registerDoSNOW(cl)
      
      # ---- progress bar setup ----
      pb <- utils::txtProgressBar(min = 0, max = totaliter, style = 3)
      progress <- function(n) utils::setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
      # ensure cleanup even if an error occurs
      on.exit({
        try(close(pb), silent = TRUE)
        try(stopCluster(cl), silent = TRUE)
      }, add = TRUE)
      
      T_all_2 <- foreach(rep_number = 1:totaliter,
                         .combine = 'rbind',
                         .options.snow = opts) %dopar% {
                           set.seed(rep_number + 1609)
                           
                           
                           Z  <- simulate.sar(N = N, rho = rho, spatial.var = spatial.var)
                           X  <- Z + mu
                           
                           alpha <- 0.5; kappa <- 0.01
                           
                           lrv_est <- lrv_estimator_outer_shell(X, alpha = 0.5, Bn_pow = 1/3)
                           resQ_est <- Q_base$Q * sqrt(lrv_est)
                           res_est <- multiple_sp_changepoints(
                             X, alpha, kappa, Q = resQ_est,
                             connectivity = 8L, plot_components = FALSE, plot_estimates = FALSE
                           )
                           
                           
                           
                           mval <- NA_real_
                           if (!is.null(res_est$K_hat) && res_est$K_hat == nrow(tau1)) {
                             tau1_hat <- do.call(rbind, lapply(res_est$I_hat, function(x) x$s / N))
                             tau2_hat <- do.call(rbind, lapply(res_est$I_hat, function(x) x$t / N))
                             
                             
                             D <- as.matrix(dist(rbind(tau1, tau1_hat)))[1:nrow(tau1), (nrow(tau1)+1):(2*nrow(tau1))]
                             
                             # Find best matches (minimum distance for each tau1)
                             ord <- apply(D, 1, which.min)
                             
                             # Reorder tau1_hat to match tau1
                             tau1_hat <- tau1_hat[ord, , drop = FALSE]
                             
                             
                             D <- as.matrix(dist(rbind(tau2, tau2_hat)))[1:nrow(tau2), (nrow(tau2)+1):(2*nrow(tau2))]
                             
                             # Find best matches (minimum distance for each tau1)
                             ord <- apply(D, 1, which.min)
                             
                             # Reorder tau1_hat to match tau1
                             tau2_hat <- tau2_hat[ord, , drop = FALSE]
                             diffs <- numeric(nrow(tau1))
                             for (k in 1:nrow(tau1)) {
                               diffs[k] <- rect_symdiff_area_bltr(
                                 tau1[k, ], tau2[k, ],
                                 tau1_hat[k,], tau2_hat[k,]
                               ) #/ area_rect(norm_rect(tau1[k, ], tau2[k, ]))
                             }
                             mval <- max(diffs)
                           }
                           
                           est_lab <- matrix(0, nrow = N, ncol = N)
                           K_hat <- res_est$K_hat
                           
                           if(K_hat ==0){RI <- NA_real_}
                           else{
                           for (k in 1:K_hat) {
                             s <- pmax(1, floor(res_est$I_hat[[k]]$s))
                             t <- pmin(N, ceiling(res_est$I_hat[[k]]$t))
                             est_lab[s[1]:t[1], s[2]:t[2]] <- k
                           }
                           
                           #--- Compute Rand Index ---#
                           vec_true <- as.vector(true_lab)
                           vec_est  <- as.vector(est_lab)
                           RI <- mclust::adjustedRandIndex(vec_true, vec_est)
                           }
                           
                           data.frame(
                             K_hat = as.integer(res_est$K_hat),
                             symdiff_max = mval,
                             RI = RI
                           )
                         }
      
      close(pb)          # close the bar
      stopCluster(cl)    # stop workers
      
      cat("\nCounts of K_hat over", nrow(T_all_2), "replicates:\n")
      print(table(T_all_2$K_hat, useNA = "ifany"))
      
      avg_when_3 <- mean(T_all_2$symdiff_max[T_all_2$K_hat == nrow(tau1)], na.rm = TRUE)
      avg_RI <- mean(T_all_2$RI, na.rm = TRUE)
      
      cat(sprintf("\nAverage of max normalized symmetric-difference when K_hat == 3: %.6f\n", avg_when_3))
      cat(sprintf("\nAverage of RI: %.6f\n", avg_RI))
      
      khat_tab <- table(T_all_2$K_hat, useNA = "ifany")
      khat_summary <- paste(names(khat_tab), khat_tab, sep = ":", collapse = ", ")
      
      datfin <- rbind(datfin, data.frame(
        N = N,
        jump = jsize,
        cor = rho,
        norm_max_sym_diff = avg_when_3,
        RI = avg_RI,
        counts_K_hat = khat_summary
      ))
      
      
      end <- Sys.time()
      print(end - start)
    }
  }
}

print(datfin) # in last column, i:j implies i number of patches is estimated j times.
