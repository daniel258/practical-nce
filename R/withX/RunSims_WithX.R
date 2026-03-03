# RunSims_WithX.R
# Run simulation replications over a grid of structural parameters (With measured covariates X1, X2).
# Requires DGM_WithX.R and FitModels_WithX.R to be sourced.

RowToPars_WithX <- function(row) {
  list(
    a0 = row$a0, a1 = row$a1, a2 = row$a2,
    b0 = row$b0, b1 = row$b1, b2 = row$b2,
    c0 = row$c0, c1 = row$c1, c2 = row$c2,
    sigma_eY = row$sigma_eY,
    ax_1 = row$ax_1, ax_2 = row$ax_2,
    bx_1 = row$bx_1, bx_2 = row$bx_2,
    p_x2 = row$p_x2
  )
}

# If `save_prefix` is provided, saves raw+agg to disk as `*_sims.rds` and `*_agg.csv`.
RunSims_WithX <- function(grid,
                          n_iters = 1000,
                          n_sample = 1000,
                          seed = 1,
                          intercept = TRUE,
                          noise_dist = c("norm", "t"),
                          df_t = 5,
                          alpha = 0.05,
                          robust_se = FALSE,
                          save_prefix = NULL) {

  noise_dist <- match.arg(noise_dist)
  if (!is.null(seed)) set.seed(seed)

  needed <- c("a0","a1","a2","b0","b1","b2","c0","c1","c2","sigma_eY",
              "ax_1","ax_2","bx_1","bx_2","p_x2")
  miss <- setdiff(needed, names(grid))
  if (length(miss) > 0) stop(paste("Grid is missing columns:", paste(miss, collapse = ", ")))

  n_grid <- nrow(grid)
  n_out  <- n_grid * n_iters

  fit_cols <- c(
    "beta_A_fit1","se_A_fit1","p_A_fit1","power_A_fit1",

    "beta_A_fit2","se_A_fit2","p_A_fit2","power_A_fit2",
    "beta_Atilde_fit2","se_Atilde_fit2","p_Atilde_fit2","power_Atilde_fit2",

    "beta_A_fit3","se_A_fit3","p_A_fit3","power_A_fit3",
    "beta_V_fit3","se_V_fit3","p_V_fit3","power_V_fit3",

    "beta_A_fit4","se_A_fit4","p_A_fit4","power_A_fit4",
    "beta_Atilde_fit4","se_Atilde_fit4","p_Atilde_fit4","power_Atilde_fit4",
    "beta_V_fit4","se_V_fit4","p_V_fit4","power_V_fit4"
  )

  raw <- grid[rep(seq_len(n_grid), each = n_iters), , drop = FALSE]
  raw$iter <- rep(seq_len(n_iters), times = n_grid)
  raw[, fit_cols] <- NA_real_
  raw <- raw[, c(names(grid), "iter", fit_cols), drop = FALSE]

  k <- 1
  for (g in seq_len(n_grid)) {
    pars <- RowToPars_WithX(grid[g, ])
    for (it in seq_len(n_iters)) {
      dat <- DGM_WithX(n_sample = n_sample, pars = pars, noise_dist = noise_dist, df = df_t)
      fit <- FitModels_WithX(dat, alpha = alpha, intercept = intercept, robust_se = robust_se)
      raw[k, fit_cols] <- fit[1, fit_cols]
      k <- k + 1
    }
  }

  mean_cols <- c(
    "beta_A_fit1","se_A_fit1","power_A_fit1",
    "beta_A_fit2","se_A_fit2","power_A_fit2",
    "beta_Atilde_fit2","se_Atilde_fit2","power_Atilde_fit2",

    "beta_A_fit3","se_A_fit3","power_A_fit3",
    "beta_V_fit3","se_V_fit3","power_V_fit3",

    "beta_A_fit4","se_A_fit4","power_A_fit4",
    "beta_Atilde_fit4","se_Atilde_fit4","power_Atilde_fit4",
    "beta_V_fit4","se_V_fit4","power_V_fit4"
  )

  sd_cols <- c("beta_A_fit1","beta_A_fit2","beta_Atilde_fit2",
               "beta_A_fit3","beta_V_fit3",
               "beta_A_fit4","beta_Atilde_fit4","beta_V_fit4")

  id_cols <- names(grid)

  agg_mean <- aggregate(raw[, mean_cols, drop = FALSE],
                        by = raw[, id_cols, drop = FALSE],
                        FUN = mean)

  agg_sd <- aggregate(raw[, sd_cols, drop = FALSE],
                      by = raw[, id_cols, drop = FALSE],
                      FUN = sd)

  names(agg_sd)[names(agg_sd) == "beta_A_fit1"]      <- "sd_beta_A_fit1"
  names(agg_sd)[names(agg_sd) == "beta_A_fit2"]      <- "sd_beta_A_fit2"
  names(agg_sd)[names(agg_sd) == "beta_Atilde_fit2"] <- "sd_beta_Atilde_fit2"
  names(agg_sd)[names(agg_sd) == "beta_A_fit3"]      <- "sd_beta_A_fit3"
  names(agg_sd)[names(agg_sd) == "beta_V_fit3"]      <- "sd_beta_V_fit3"
  names(agg_sd)[names(agg_sd) == "beta_A_fit4"]      <- "sd_beta_A_fit4"
  names(agg_sd)[names(agg_sd) == "beta_Atilde_fit4"] <- "sd_beta_Atilde_fit4"
  names(agg_sd)[names(agg_sd) == "beta_V_fit4"]      <- "sd_beta_V_fit4"

  agg <- merge(agg_mean, agg_sd, by = id_cols, all.x = TRUE)

  # Empirical bias summaries for beta_A estimates (relative to true b2)
  agg$bias_A_fit1 <- agg$beta_A_fit1 - agg$b2
  agg$bias_A_fit2 <- agg$beta_A_fit2 - agg$b2
  agg$bias_A_fit3 <- agg$beta_A_fit3 - agg$b2
  agg$bias_A_fit4 <- agg$beta_A_fit4 - agg$b2

  raw$source   <- "sim"
  agg$source   <- "sim"
  raw$se_type  <- if (robust_se) "HC1" else "classical"
  agg$se_type  <- if (robust_se) "HC1" else "classical"
  raw$noise_dist <- noise_dist
  agg$noise_dist <- noise_dist
  raw$df_t <- if (noise_dist == "t") df_t else NA_real_
  agg$df_t <- if (noise_dist == "t") df_t else NA_real_
  raw$n_sample <- n_sample
  agg$n_sample <- n_sample
  raw$n_iters <- n_iters
  agg$n_iters <- n_iters
  raw$seed <- if (is.null(seed)) NA_real_ else seed
  agg$seed <- if (is.null(seed)) NA_real_ else seed

  if (!is.null(save_prefix)) {
    dir.create(dirname(save_prefix), showWarnings = FALSE, recursive = TRUE)
    saveRDS(list(raw = raw, agg = agg), file = paste0(save_prefix, "_sims.rds"))
    write.csv(agg, file = paste0(save_prefix, "_agg.csv"), row.names = FALSE)
  }

  list(raw = raw, agg = agg)
}
