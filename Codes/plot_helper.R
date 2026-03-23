library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)

# Convert matrix to data frame for ggplot
matrix_to_df <- function(M) {
  nr <- nrow(M)
  nc <- ncol(M)
  
  expand.grid(
    row_frac = seq(0, 1, length.out = nr),
    col_frac = seq(0, 1, length.out = nc)
  ) %>%
    mutate(value = as.vector(M))
}

# Convert index-based rectangle to fraction-based
indices_to_rect <- function(i1, i2, j1, j2, nr, nc) {
  data.frame(
    xmin = (i1 - 1) / nc,
    xmax = i2 / nc,
    ymin = (j1 - 1) / nr,
    ymax = j2 / nr
  )
}

# Main plotting function
plot_detection_summary <- function(X, tau1, tau2, res,
                                   line_width = 2.5) {
  n1 <- nrow(X)
  n2 <- ncol(X)
  
  # Prepare tau matrices
  tau1m <- if (is.numeric(tau1) && length(tau1) == 2) 
    matrix(tau1, ncol = 2, byrow = TRUE) else as.matrix(tau1)
  tau2m <- if (is.numeric(tau2) && length(tau2) == 2) 
    matrix(tau2, ncol = 2, byrow = TRUE) else as.matrix(tau2)
  
  # Convert matrix to data frame
  df_matrix <- matrix_to_df(X)
  
  # Prepare truth rectangles from tau
  df_truth <- data.frame(
    xmin = tau1m[, 1],
    xmax = tau2m[, 1],
    ymin = tau1m[, 2],
    ymax = tau2m[, 2],
    type = "True change-patches"
  )
  
  # Prepare estimated rectangles (orange)
  df_est <- data.frame()
  if (!is.null(res$I_bounds) && nrow(res$I_bounds) > 0) {
    df_est <- do.call(rbind, lapply(seq_len(nrow(res$I_bounds)), function(k) {
      indices_to_rect(res$I_bounds[k, "i1"], res$I_bounds[k, "i2"],
                      res$I_bounds[k, "j1"], res$I_bounds[k, "j2"],
                      nr = n1, nc = n2)
    }))
  } else if (length(res$I_hat) > 0) {
    df_est <- do.call(rbind, lapply(res$I_hat, function(r) {
      indices_to_rect(r$s[1], r$t[1], r$s[2], r$t[2],
                      nr = n1, nc = n2)
    }))
  }
  if (nrow(df_est) > 0) {
    df_est$type <- "Estimated patches"
  }
  
  # # Prepare oracle rectangles (cyan)
  # df_oracle <- data.frame()
  # if (!is.null(res_oracle$I_bounds) && nrow(res_oracle$I_bounds) > 0) {
  #   df_oracle <- do.call(rbind, lapply(seq_len(nrow(res_oracle$I_bounds)), function(k) {
  #     indices_to_rect(res_oracle$I_bounds[k, "i1"], res_oracle$I_bounds[k, "i2"],
  #                     res_oracle$I_bounds[k, "j1"], res_oracle$I_bounds[k, "j2"],
  #                     nr = n1, nc = n2)
  #   }))
  # } else if (length(res_oracle$I_hat) > 0) {
  #   df_oracle <- do.call(rbind, lapply(res_oracle$I_hat, function(r) {
  #     indices_to_rect(r$s[1], r$t[1], r$s[2], r$t[2],
  #                     nr = n1, nc = n2)
  #   }))
  # }
  # if (nrow(df_oracle) > 0) {
  #   df_oracle$type <- "Estimated patches with oracle LRV"
  # }
  # 
  # Combine all rectangles
  df_rects <- rbind(df_truth, df_est)
  df_rects$type <- factor(df_rects$type, 
                          levels = c("True change-patches",
                                     "Estimated patches")
  )
  
  # Define colors
  rect_colors <- c(
    "True change-patches" = "blue",  # magenta
    "Estimated patches" = "red"  # orange
    #"Estimated patches with oracle LRV" = "#00FFFF"  # cyan
  )
  
  # Create plot
  p <- ggplot() +
    geom_raster(data = df_matrix, aes(x = col_frac, y = row_frac, fill = value)) +
    scale_fill_viridis_c(name = "") +
    geom_rect(data = df_rects, 
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                  color = type),
              fill = NA, linewidth = line_width) +
    scale_color_manual(name = "", values = rect_colors) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), 
                       labels = sprintf("%.1f", seq(0, 1, by = 0.1)),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                       labels = sprintf("%.1f", seq(0, 1, by = 0.1)),
                       expand = c(0, 0)) +
    labs(x = "Row fraction",
         y = "Column fraction") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = c(0.99, 0.01),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key.height = unit(0.4, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      aspect.ratio = n1/n2
    ) +
    guides(
      fill = guide_colorbar(order = 1),
      color = guide_legend(order = 2, override.aes = list(linewidth = 2))
    )
  
  return(p)
}

## example usage

#plot_detection_summary(t(X), tau1 = tau1, tau2 = tau2, res = res_est)

