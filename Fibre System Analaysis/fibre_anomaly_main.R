source("preprocess_fibre_data.R")
# Source the 3D algorithm files (adjust paths as needed)
source("multiple_cp_helper_3d.R")
source("LRV_helper_3d.R")
source("quantile_calc_simu_3d.R")
source("plot_helper_3d.R")

# 1) Load and preprocess
fdat <- load_fibre_data("PP_LGF60_schnell_100_80_4um_Crop_iass_Radius3_fiberdirections.csv")
N <- fdat$grid_dims

alpha = 0.5; kappa = 0.01
attribute = "y_abs"
nsim_Q = 500; connectivity = 26L

# 2) Select which attribute to analyze
X <- switch(attribute,
            "x"     = fdat$X_x_filled,
            "y"     = fdat$X_y_filled,
            "z"     = fdat$X_z_filled,
            "x_abs" = fdat$X_x_abs_filled,
            "y_abs" = fdat$X_y_abs_filled,
            "z_abs" = fdat$X_z_abs_filled,
            "aniso" = fdat$X_aniso,
            stop("Unknown attribute: ", attribute)
)

cat("\nAnalyzing attribute:", attribute, "\n")
cat("Grid:", N[1], "x", N[2], "x", N[3], "\n")

# 3) LRV estimation
cat("Estimating LRV...\n")
lrv_est <- lrv_estimator_fft_3d(X, alpha = alpha, Bn_pow = 1/3)

# 4) Threshold Q
cat("Simulating Q threshold...\n")
resQ <- simulate_Q_blockmeans_3d(n = N, alpha = alpha, sigma2 = lrv_est,
                                 nsim = nsim_Q, seed = 123)
cat("Q =", resQ$Q, "\n")

# 5) Run Algorithm 2
cat("Running multiple changepoint detection...\n")
start <- Sys.time()
res <- multiple_sp_changepoints_3d(X, alpha, gamma, Q = resQ$Q,
                                   connectivity = connectivity)
elapsed <- Sys.time() - start
cat("Done in", format(elapsed), "\n")
cat("K_hat =", res$K_hat, "\n")

if (res$K_hat > 0) {
  cat("\nEstimated patches (grid indices):\n")
  print(res$I_bounds)
  cat("\nEstimated patches (fractions):\n")
  print(round(res$I_bounds / rep(c(N, N), each = nrow(res$I_bounds)), 3))
}

# 6) Plot
if (res$K_hat > 0) {
  pct <- 0.01
  est_boxes <- res$I_hat
  cat("\nGenerating 3D plot...\n")
  png(paste0("fibre_detection_real_", attribute, " ", pct, ".png"), width = 800, height = 700, res = 100)
  plot_detection_3d_static(X, true_boxes = NULL, est_boxes = est_boxes,
                           subsample = pct, est_col = "black",
                           main = paste0("")
                           )
  dev.off()
  cat("Saved plot.\n")
}

