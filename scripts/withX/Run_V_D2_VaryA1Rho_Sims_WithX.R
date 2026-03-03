# Run_V_D2_VaryA1Rho_Sims_WithX.R
# Design D2 sweep over (a1, rho_total), with two measured covariates.
# Within each (a1, rho_total): a1 fixed (=> bias fixed up to a constant given X adjustment), rho_total fixed,
# and vary f in [0,1] retaining only feasible f values.

# -------------------- user inputs --------------------
n_sample <- 1000
n_iters  <- 1000
seed     <- 314

noise_dist <- "norm"
df_t       <- 5
alpha      <- 0.05
robust_se  <- FALSE

# Measured covariates affect A and Y, not Atilde
ax_1 <- 0.15
ax_2 <- 0.15
bx_1 <- 0.2
bx_2 <- 0.2
p_x2 <- 0.5

# Full f grid; feasibility handled by gridmaker
f_vec <- seq(0, 1, by = 0.05)

# Sweep a1 values (a2=a1). Must satisfy 2*a1^2 + ax_1^2 + ax_2^2 < 1
a1_vec <- c(0.2, 0.4, 0.6)

# Sweep rho_total values
rho_total_vec <- c(0.2, 0.4, 0.6)

# Outcome model parameters
b0 <- 0
b1 <- 0.30
b2 <- 0.30
sigma_eY <- 1

out_dir <- "results/withX/withV/D2"
design_label <- "D2_varyA1Rho_WithX"

# -------------------- load project code --------------------
source("R/withX/DGM_WithX.R")
source("R/withX/FitModels_WithX.R")
source("R/withX/MakeGrids_WithX.R")
source("R/withX/RunSims_WithX.R")

# -------------------- output prefix  ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
tag <- sprintf("withV_D2_varyA1Rho_withX_%s_n%d_it%d_seed%d", noise_dist, n_sample, n_iters, seed)
save_prefix <- file.path(out_dir, paste0(tag, "_", stamp))

# -------------------- build + run over (a1, rho_total) pairs --------------------
skip_log <- data.frame(a1 = numeric(0), rho_total = numeric(0), reason = character(0), stringsAsFactors = FALSE)
agg_list <- list()
raw_list <- list()
k <- 1

t0 <- proc.time()[3]

for (a1 in a1_vec) {
  a2 <- a1

  # quick feasibility check for A given X
  if ((a1^2 + a2^2 + ax_1^2 + ax_2^2) >= 1) {
    skip_log <- rbind(skip_log, data.frame(a1 = a1, rho_total = NA_real_, reason = "Infeasible A variance with X", stringsAsFactors = FALSE))
    next
  }

  for (rho_total in rho_total_vec) {
    grid_try <- tryCatch({
      MakeGrid_D2_FixBiasFixRho_VaryF_WithX(
        rho_total = rho_total,
        f_vec = f_vec,
        a1 = a1,
        a2 = a2,
        b0 = b0, b1 = b1, b2 = b2,
        ax_1 = ax_1, ax_2 = ax_2, bx_1 = bx_1, bx_2 = bx_2, p_x2 = p_x2,
        a0 = 0, c0 = 0,
        sigma_eY = sigma_eY,
        add_derived = TRUE
      )
    }, error = function(e) e)

    if (inherits(grid_try, "error")) {
      skip_log <- rbind(skip_log,
                        data.frame(a1 = a1, rho_total = rho_total, reason = grid_try$message, stringsAsFactors = FALSE))
      next
    }

    grid <- grid_try
    grid$design <- design_label
    grid$a2_rule <- "a2_equal_a1"

    res <- RunSims_WithX(
      grid = grid,
      n_sample = n_sample,
      n_iters = n_iters,
      seed = seed,
      intercept = TRUE,
      noise_dist = noise_dist,
      df_t = df_t,
      alpha = alpha,
      robust_se = robust_se,
      save_prefix = NULL
    )

    agg_list[[k]] <- res$agg
    raw_list[[k]] <- res$raw
    k <- k + 1
  }
}

elapsed_sec <- proc.time()[3] - t0

if (length(agg_list) == 0) stop("All (a1, rho_total) pairs were skipped; nothing to save.")

agg_all <- do.call(rbind, agg_list)
raw_all <- do.call(rbind, raw_list)

saveRDS(list(raw = raw_all, agg = agg_all), file = paste0(save_prefix, "_sims.rds"))
write.csv(agg_all, file = paste0(save_prefix, "_agg.csv"), row.names = FALSE)

manifest <- list(
  created_at = as.character(Sys.time()),
  script = "Run_V_D2_VaryA1Rho_Sims_WithX.R",
  save_prefix = save_prefix,
  n_sample = n_sample,
  n_iters = n_iters,
  seed = seed,
  noise_dist = noise_dist,
  df_t = df_t,
  alpha = alpha,
  robust_se = robust_se,
  f_vec = f_vec,
  a1_vec = a1_vec,
  rho_total_vec = rho_total_vec,
  b0 = b0, b1 = b1, b2 = b2,
  sigma_eY = sigma_eY,
  ax_1 = ax_1, ax_2 = ax_2, bx_1 = bx_1, bx_2 = bx_2, p_x2 = p_x2,
  grid_dim = dim(agg_all),
  grid_preview = utils::head(agg_all, 10),
  session_info = utils::sessionInfo(),
  runtime_minutes = elapsed_sec / 60
)

saveRDS(manifest, file = paste0(save_prefix, "_manifest.rds"))
writeLines(capture.output(str(manifest, max.level = 2)),
           con = paste0(save_prefix, "_manifest.txt"))
write.csv(skip_log, file = paste0(save_prefix, "_skipped.csv"), row.names = FALSE)

message("Saved: ", paste0(save_prefix, "_sims.rds"))
message("Saved: ", paste0(save_prefix, "_agg.csv"))
message("Saved: ", paste0(save_prefix, "_manifest.rds / _manifest.txt"))
message("Saved: ", paste0(save_prefix, "_skipped.csv"))

invisible(list(raw = raw_all, agg = agg_all, skipped = skip_log))
