# Run_V_D1.R
# Design 1 (with V): Fix total corr(A, Atilde) and fix bias; vary decomposition via V vs U.


# ---- user inputs ----
n_sample <- 1000
n_iters  <- 1000
seed     <- 314

noise_dist <- "norm"  # "norm" or "t"
df_t       <- 5
alpha      <- 0.05
robust_se  <- FALSE

# Output
out_dir <- "results/withV/D1"

# ---- Design 1 settings ----
# Fix a1,a2 then choose rho_total targets; vary c1 and solve for c2 in MakeGrid_D1_FixRhoFixBias()
a1 <- 0.5
a2 <- 0.5

rho_total_vec <- c(0.2, 0.4, 0.6)         # targets for rho(A, Atilde)
bias_vec <- c(0.05, 0.10, 0.15, 0.20, 0.25)  # target "bias" parameter used by MakeGrid 
c1_vec <- seq(0.05, 0.95, by = 0.05)       # varying contribution of V to corr

# Outcome model parameters (as in your grids)
b2 <- 0.3
sigma_eY <- 1

# ---- load project functions ----
source("R/DGM.R")
source("R/FitModels.R")
source("R/MakeGrids.R")
source("R/RunSims.R")

# ---- build grid ----
tol <- 1e-12
grid_list <- list()
k <- 1

for (rho_total in rho_total_vec) {
  for (bias in bias_vec) {
    
    g <- MakeGrid_D1_FixRhoFixBias(
      c1_vec      = c1_vec,
      a1          = a1,
      a2          = a2,
      rho_total   = rho_total,
      bias        = bias,
      b2          = b2,
      sigma_eY    = sigma_eY,
      add_derived = TRUE
    )
    
    # Skip infeasible cells (these become blank regions in the slide-style figure)
    if (nrow(g) == 0) next
    
    g$rho_target  <- rho_total
    g$bias_target <- bias
    g$design      <- "D1_FixRhoFixBias"
    
    grid_list[[k]] <- g
    k <- k + 1
  }
}

grid <- do.call(rbind, grid_list)

if (is.null(grid) || nrow(grid) == 0) stop("No feasible parameter combinations in grid.")


# ---- run sims ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
tag <- sprintf("withV_D1_%s_n%d_it%d_seed%d", noise_dist, n_sample, n_iters, seed)
save_prefix <- file.path(out_dir, paste0(tag, "_", stamp))

t0 <- proc.time()
res <- RunSims(
  grid = grid,
  n_iters = n_iters,
  n_sample = n_sample,
  seed = seed,
  noise_dist = noise_dist,
  df_t = df_t,
  alpha = alpha,
  robust_se = robust_se,
  save_prefix = save_prefix
)
elapsed_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("[D1] Total runtime: %.1f seconds (%.2f minutes)", elapsed_sec, elapsed_sec/60))

# ---- manifest ----
manifest <- list(
  created_at = as.character(Sys.time()),
  script = "Run_V_D1.R",
  design = "D1",
  save_prefix = save_prefix,
  n_sample = n_sample,
  n_iters = n_iters,
  seed = seed,
  noise_dist = noise_dist,
  df_t = df_t,
  alpha = alpha,
  robust_se = robust_se,
  a1 = a1,
  a2 = a2,
  rho_total_vec = rho_total_vec,
  bias = bias,
  c1_vec = c1_vec,
  b2 = b2,
  sigma_eY = sigma_eY,
  grid_dim = dim(grid),
  grid_preview = utils::head(grid, 10),
  session_info = utils::sessionInfo(),
  runtime_minutes = elapsed_sec / 60
)

saveRDS(manifest, file = paste0(save_prefix, "_manifest.rds"))
writeLines(capture.output(str(manifest, max.level = 2)),
           con = paste0(save_prefix, "_manifest.txt"))

message("Saved: ", paste0(save_prefix, "_sims.rds"))
message("Saved: ", paste0(save_prefix, "_agg.csv"))
message("Saved: ", paste0(save_prefix, "_manifest.rds / _manifest.txt"))

invisible(res)
