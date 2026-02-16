# Run_V_D2_VaryA1Rho_Sims.R
# Design D2:
#   - Within each (a1, rho_total): a1 fixed (=> bias b1*a1 fixed), rho_total fixed
#   - Vary f in [0,1], keep only feasible f values (gridmaker drops infeasible rows)
# Sweep:
#   - multiple a1 values (always a2=a1)
#   - multiple rho_total values

# -------------------- user inputs --------------------
n_sample <- 1000
n_iters  <- 1000
seed     <- 314

noise_dist <- "norm"   # "norm" or "t"
df_t       <- 5
alpha      <- 0.05
robust_se  <- FALSE

# Full f grid; feasibility handled by gridmaker
f_vec <- seq(0, 1, by = 0.05)

# Sweep a1 values (a2=a1). Must satisfy 2*a1^2 < 1 => a1 < 0.7071
a1_vec <- c(0.2, 0.4, 0.6)

# Sweep rho_total values (some f endpoints may be infeasible depending on a1)
rho_total_vec <- c(0.2, 0.4, 0.6)

# Outcome model parameters
b0 <- 0
b1 <- 0.3
b2 <- 0.3
sigma_eY <- 1

# Intercepts (usually 0)
a0 <- 0
c0 <- 0

out_dir <- "results/withV/D2"
design_label <- "D2_varyA1Rho"
# -----------------------------------------------------

# -------------------- load project code --------------------
source("R/DGM.R")
source("R/FitModels.R")
source("R/MakeGrids.R")
source("R/RunSims.R")

# -------------------- output prefix  ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
tag <- sprintf("withV_D2_varyA1Rho_%s_n%d_it%d_seed%d", noise_dist, n_sample, n_iters, seed)
save_prefix <- file.path(out_dir, paste0(tag, "_", stamp))

# Safe wrapper: returns NULL if the D2 gridmaker errors (e.g., empty grid)
SafeMakeGrid_D2 <- function(...) {
  out <- tryCatch(
    MakeGrid_D2_FixBiasFixRho_VaryF(...),
    error = function(e) {
      msg <- conditionMessage(e)
      message("  [skip] ", msg)
      return(NULL)
    }
  )
  out
}
# -------------------- build combined grid --------------------
grid_list <- list()
k <- 1

skip_log <- data.frame(a1 = numeric(0), rho = numeric(0), reason = character(0))

for (a1 in a1_vec) {
  a2 <- a1
  
  if ((a1^2 + a2^2) >= 1) {
    message(sprintf("Skipping a1=%.2f (infeasible: a1^2+a2^2>=1).", a1))
    skip_log <- rbind(skip_log, data.frame(a1=a1, rho=NA, reason="a1^2+a2^2>=1"))
    next
  }
  
  for (rho in rho_total_vec) {
    
    g <- SafeMakeGrid_D2(
      rho_total = rho,
      f_vec     = f_vec,
      a1        = a1,
      a2        = a2,
      b1        = b1,
      b2        = b2,
      a0        = a0,
      b0        = b0,
      c0        = c0,
      sigma_eY  = sigma_eY,
      add_derived = TRUE
    )
    
    if (is.null(g)) {
      skip_log <- rbind(skip_log, data.frame(a1=a1, rho=rho, reason="gridmaker error/empty"))
      next
    }
    
    kept_f <- sort(unique(g$f_target))
    message(sprintf(
      "a1=%.2f, rho=%.2f: kept %d/%d f values; f range [%.2f, %.2f]",
      a1, rho, length(kept_f), length(f_vec), min(kept_f), max(kept_f)
    ))
    
    g$design <- design_label
    grid_list[[k]] <- g
    k <- k + 1
  }
}

# If nothing feasible, still write skip log and stop.
if (length(grid_list) == 0) {
  write.csv(skip_log, file = paste0(save_prefix, "_skipped.csv"), row.names = FALSE)
  stop("All (a1, rho_total) combinations were infeasible; nothing to run.")
}

grid <- do.call(rbind, grid_list)
rownames(grid) <- NULL

# -------------------- run --------------------

t0 <- proc.time()

res <- RunSims(
  grid        = grid,
  n_iters     = n_iters,
  n_sample    = n_sample,
  seed        = seed,
  noise_dist  = noise_dist,
  df_t        = df_t,
  alpha       = alpha,
  robust_se   = robust_se,
  save_prefix = save_prefix
)

elapsed_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("[Run_V_D2_VaryA1Rho_Sims] Total runtime: %.1f sec (%.2f min)", elapsed_sec, elapsed_sec / 60))

# -------------------- manifest --------------------
manifest <- list(
  created_at      = as.character(Sys.time()),
  script          = "Run_V_D2_VaryA1Rho_Sims.R",
  design_label    = design_label,
  save_prefix     = save_prefix,
  n_sample        = n_sample,
  n_iters         = n_iters,
  seed            = seed,
  noise_dist      = noise_dist,
  df_t            = df_t,
  alpha           = alpha,
  robust_se       = robust_se,
  a1_vec          = a1_vec,
  rho_total_vec   = rho_total_vec,
  f_vec           = f_vec,
  a0              = a0,
  c0              = c0,
  b0              = b0,
  b1              = b1,
  b2              = b2,
  sigma_eY        = sigma_eY,
  grid_dim        = dim(grid),
  grid_preview    = utils::head(grid, 10),
  session_info    = utils::sessionInfo(),
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

invisible(res)
