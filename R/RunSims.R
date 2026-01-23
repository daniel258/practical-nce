# RunSims.R
# Run simulation replications over a grid of structural parameters.
# Requires DGM.R and FitModels.R to be sourced
source("DGM.R")
source("FitModels.R")

# Simple grid builder for the "no V" case (set a2=c2=0).
MakeGrid_NoV <- function(a1, c1, b1, b2, sigma_eY = 1) {
  grid <- expand.grid(a1 = a1, c1 = c1, b1 = b1, b2 = b2, stringsAsFactors = FALSE)
  grid$a0 <- 0; grid$b0 <- 0; grid$c0 <- 0
  grid$a2 <- 0; grid$c2 <- 0
  grid$sigma_eY <- sigma_eY
  grid
}

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
    "beta_A_fit1","se_A_fit1","p_A_fit1","power_A_fit1",
    "beta_A_fit2","se_A_fit2","p_A_fit2","power_A_fit2",
    "beta_Atilde_fit2","se_Atilde_fit2","p_Atilde_fit2","power_Atilde_fit2"
  )
  # make space for results
  raw <- data.frame(matrix(NA_real_, nrow = n_out,
                           ncol = ncol(grid) + 1 + length(fit_cols)))
  names(raw) <- c(names(grid), "iter", fit_cols)
  
  k <- 1 # row index for raw results
  for (g in seq_len(n_grid)) {
    pars <- RowToPars(grid[g, ])
    
    for (it in seq_len(n_iters)) {
      dat <- DGM(n_sample = n_sample, pars = pars, noise_dist = noise_dist, df = df_t)
      fit <- FitModels(dat, alpha = alpha, intercept = intercept, robust_se = robust_se)
      
      raw[k, names(grid)] <- grid[g, ]
      raw[k, "iter"] <- it
      raw[k, fit_cols] <- fit[1, fit_cols]
      k <- k + 1
    }
  }
  
  mean_cols <- c("beta_A_fit1","se_A_fit1","power_A_fit1",
                 "beta_A_fit2","se_A_fit2","power_A_fit2",
                 "beta_Atilde_fit2","se_Atilde_fit2","power_Atilde_fit2")
  sd_cols <- c("beta_A_fit1","beta_A_fit2","beta_Atilde_fit2")
  agg_mean <- aggregate(raw[, mean_cols, drop = FALSE],
                        by = raw[, needed, drop = FALSE],
                        FUN = mean)
  
  agg_sd <- aggregate(raw[, sd_cols, drop = FALSE],
                      by = raw[, needed, drop = FALSE],
                      FUN = sd)
  
  names(agg_sd)[names(agg_sd) == "beta_A_fit1"]      <- "sd_beta_A_fit1"
  names(agg_sd)[names(agg_sd) == "beta_A_fit2"]      <- "sd_beta_A_fit2"
  names(agg_sd)[names(agg_sd) == "beta_Atilde_fit2"] <- "sd_beta_Atilde_fit2"
  
  agg <- merge(agg_mean, agg_sd, by = needed, all.x = TRUE)
  
  
  # Bias summaries for beta_A (relative to structural b2 in your DGM)
  agg$bias_A_fit1 <- agg$beta_A_fit1 - agg$b2
  agg$bias_A_fit2 <- agg$beta_A_fit2 - agg$b2
  
  if (!is.null(save_prefix)) {
    saveRDS(list(raw = raw, agg = agg),
            file = paste0(save_prefix, "_sims.rds"))
    write.csv(agg, file = paste0(save_prefix, "_agg.csv"), row.names = FALSE)
  }
  
  list(raw = raw, agg = agg)
}

# ---- minimal example ----

# grid0 <- MakeGrid_NoV(
#   a1 = c(0.2, 0.4, 0.6),
#   c1 = c(0.2, 0.4, 0.6),
#   b1 = c(0.2, 0.4),
#   b2 = 0.3,
#   sigma_eY = 1
# )
# 
# res0 <- RunSims(
#   grid = grid0,
#   n_sample = 1000,
#   n_iters = 200,
#   robust_se = FALSE,
#   noise_dist = "norm",
#   seed = 1,
#   save_prefix = "results/noV_example"
# )
# 
# head(res0$agg)
