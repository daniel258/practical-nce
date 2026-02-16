# Run_V_D2_HighRho_Sims.R
# Design D2, but for high rho_total (0.70--0.95):
#   - keep a1 moderate (controls confounding/bias magnitude via b1*a1)
#   - choose a2 as large as possible while keeping a1^2 + a2^2 < 1
#   - vary f in [0,1], keep only feasible f's (gridmaker drops infeasible rows)
#   - skip (a1, rho) pairs with empty grids using tryCatch
#   - ALWAYS save a skip log (even if nothing is feasible)

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

# Moderate a1 values (bias ~ b1*a1). Keep these "not too extreme".
a1_vec <- c(0.25, 0.35, 0.45)

# High rho_total values to explore
rho_total_vec <- c(0.70, 0.80, 0.90, 0.95)

# Outcome model parameters
b0 <- 0
b1 <- 0.30
b2 <- 0.30
sigma_eY <- 1

# Intercepts (usually 0)
a0 <- 0
c0 <- 0

# Choose how close to the boundary you want a2 to be.
# a2 is set to sqrt(1 - a1^2) - a2_margin, capped at a2_cap.
a2_margin <- 0.03
a2_cap    <- 0.97

out_dir <- "results/withV/D2"
design_label <- "D2_highRho_unequalA2"
# -----------------------------------------------------

# -------------------- load project code --------------------
source("R/DGM.R")
source("R/FitModels.R")
source("R/MakeGrids.R")
source("R/RunSims.R")

# -------------------- output prefix  ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
tag <- sprintf("withV_D2_highRho_unequalA2_%s_n%d_it%d_seed%d",
               noise_dist, n_sample, n_iters, seed)
save_prefix <- file.path(out_dir, paste0(tag, "_", stamp))
# -----------------------------------------------------------------------------------------------

# Safe wrapper: return NULL if gridmaker errors (e.g., empty grid)
SafeMakeGrid_D2 <- function(...) {
  tryCatch(
    MakeGrid_D2_FixBiasFixRho_VaryF(...),
    error = function(e) {
      message("  [skip] ", conditionMessage(e))
      return(NULL)
    }
  )
}

# Choose a2 given a1 (as large as possible without violating a1^2+a2^2<1)
ChooseA2 <- function(a1, margin = a2_margin, cap = a2_cap) {
  a2_max <- sqrt(max(0, 1 - a1^2)) - margin
  a2 <- min(a2_max, cap)
  if (a2 <= 0) return(NA_real_)
  if ((a1^2 + a2^2) >= 1) {
    # push slightly inward if numerical issues
    a2 <- sqrt(max(0, 1 - a1^2)) - (margin + 1e-6)
  }
  a2
}

grid_list <- list()
k <- 1
skip_log <- data.frame(a1 = numeric(0), a2 = numeric(0), rho = numeric(0), reason = character(0))

for (a1 in a1_vec) {
  
  a2 <- ChooseA2(a1)
  if (is.na(a2) || (a1^2 + a2^2) >= 1) {
    msg <- sprintf("a1=%.2f gives invalid a2 (%.3f).", a1, a2)
    message("Skipping: ", msg)
    skip_log <- rbind(skip_log, data.frame(a1 = a1, a2 = a2, rho = NA, reason = msg))
    next
  }
  
  message(sprintf("a1=%.2f -> using a2=%.3f (Var(eA)=%.3f)",
                  a1, a2, 1 - a1^2 - a2^2))
  
  for (rho in rho_total_vec) {
    
    g <- SafeMakeGrid_D2(
      rho_total   = rho,
      f_vec       = f_vec,
      a1          = a1,
      a2          = a2,
      b1          = b1,
      b2          = b2,
      a0          = a0,
      b0          = b0,
      c0          = c0,
      sigma_eY    = sigma_eY,
      add_derived = TRUE
    )
    
    if (is.null(g) || nrow(g) == 0) {
      skip_log <- rbind(skip_log, data.frame(a1 = a1, a2 = a2, rho = rho, reason = "empty grid"))
      next
    }
    
    # Keep track of a2 explicitly (so it's in agg.csv for plotting/diagnostics)
    g$a2 <- a2
    
    kept_f <- sort(unique(g$f_target))
    message(sprintf(
      "  rho=%.2f: kept %d/%d f values; f range [%.2f, %.2f]",
      rho, length(kept_f), length(f_vec), min(kept_f), max(kept_f)
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
message(sprintf("[Run_V_D2_HighRho_Sims] Total runtime: %.1f sec (%.2f min)",
                elapsed_sec, elapsed_sec / 60))

# -------------------- manifest + skipped combos --------------------
manifest <- list(
  created_at      = as.character(Sys.time()),
  script          = "Run_V_D2_HighRho_Sims.R",
  design_label    = design_label,
  save_prefix     = save_prefix,
  n_sample        = n_sample,
  n_iters         = n_iters,
  seed            = seed,
  noise_dist      = noise_dist,
  df_t            = df_t,
  alpha           = alpha,
  robust_se       = robust_se,
  f_vec           = f_vec,
  a1_vec          = a1_vec,
  rho_total_vec   = rho_total_vec,
  a2_margin       = a2_margin,
  a2_cap          = a2_cap,
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