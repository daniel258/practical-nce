# RunSims.R
# Run simulation replications over a grid of structural parameters.
# Requires DGM.R and FitModels.R to be sourced

RowToPars <- function(row) {
  list(
    a0 = row$a0, a1 = row$a1, a2 = row$a2,
    b0 = row$b0, b1 = row$b1, b2 = row$b2,
    c0 = row$c0, c1 = row$c1, c2 = row$c2,
    sigma_eY = row$sigma_eY
  )
}

# Main simulation function.
# If `save_prefix` is provided (e.g., `"results/noV_example"`), the function saves the raw and aggregated simulation results to disk as `*_sims.rds` and `*_agg.csv`; if it is `NULL`, nothing is saved.
RunSims <- function(grid,
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
  
  needed <- c("a0","a1","a2","b0","b1","b2","c0","c1","c2","sigma_eY")
  miss <- setdiff(needed, names(grid))
  if (length(miss) > 0) stop(paste("Grid is missing columns:", paste(miss, collapse = ", ")))
  
  
  n_grid <- nrow(grid)
  n_out  <- n_grid * n_iters
  # columns to extract from FitModels output
  fit_cols <- c(
    # fit1: Y ~ A
    "beta_A_fit1","se_A_fit1","p_A_fit1","power_A_fit1",
    
    # fit2: Y ~ A + Atilde
    "beta_A_fit2","se_A_fit2","p_A_fit2","power_A_fit2",
    "beta_Atilde_fit2","se_Atilde_fit2","p_Atilde_fit2","power_Atilde_fit2",
    
    # fit3: Y ~ A + V
    "beta_A_fit3","se_A_fit3","p_A_fit3","power_A_fit3",
    "beta_V_fit3","se_V_fit3","p_V_fit3","power_V_fit3",
    
    # fit4: Y ~ A + Atilde + V
    "beta_A_fit4","se_A_fit4","p_A_fit4","power_A_fit4",
    "beta_Atilde_fit4","se_Atilde_fit4","p_Atilde_fit4","power_Atilde_fit4",
    "beta_V_fit4","se_V_fit4","p_V_fit4","power_V_fit4"
  )
  
  # make space for results
  # --- make space for results (preserve grid column types) ---
  raw <- grid[rep(seq_len(n_grid), each = n_iters), , drop = FALSE]
  raw$iter <- rep(seq_len(n_iters), times = n_grid)
  
  # pre-allocate fit outputs
  raw[, fit_cols] <- NA_real_
  
  # optional: enforce column order
  raw <- raw[, c(names(grid), "iter", fit_cols), drop = FALSE]
  
  k <- 1
  for (g in seq_len(n_grid)) {
    pars <- RowToPars(grid[g, ])
    for (it in seq_len(n_iters)) {
      dat <- DGM(n_sample = n_sample, pars = pars, noise_dist = noise_dist, df = df_t)
      fit <- FitModels(dat, alpha = alpha, intercept = intercept, robust_se = robust_se)
      
      # only fill the fit columns; grid + iter already set
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
  
  
  # Bias summaries for beta_A estimates
  agg$bias_A_fit1 <- agg$beta_A_fit1 - agg$b2
  agg$bias_A_fit2 <- agg$beta_A_fit2 - agg$b2
  agg$bias_A_fit3 <- agg$beta_A_fit3 - agg$b2
  agg$bias_A_fit4 <- agg$beta_A_fit4 - agg$b2
  
  # add summaries of DGM and analyses
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
    saveRDS(list(raw = raw, agg = agg),
            file = paste0(save_prefix, "_sims.rds"))
    write.csv(agg, file = paste0(save_prefix, "_agg.csv"), row.names = FALSE)
  }
  
  list(raw = raw, agg = agg)
}

#---- minimal example ----
# 
# source("DGM.R")
# source("FitModels.R")
# all_a1 <- c(0.2, 0.4, 0.6)
# all_b1 <- c(0.2, 0.4)
# all_c1 <- c(0.2, 0.4, 0.6)
# all_b2 <- 0.3
# all_a0 <- 0 ;all_b0 <- 0;all_c0 <- 0
# all_a2 <- 0;all_c2 <- 0
# all_sigma_eY <- 1
# grid0 <- expand.grid(a0 = all_a0, b0 = all_b0, c0 = all_c0,
#                     a1 = all_a1, c1 = all_c1, b1 = all_b1,
#                     a2 = all_a2, b2 = all_b2, c2 = all_c2,
#                     sigma_eY = all_sigma_eY, stringsAsFactors = FALSE)
# 
# 
# #
# res0 <- RunSims(
#   grid = grid0,
#   n_sample = 1000,
#   n_iters = 200,
#   robust_se = FALSE,
#   noise_dist = "norm",
#   seed = 314,
#   save_prefix = "../results/noV_example"
# )
# 
# head(res0$agg) 
