# Run_V_D1.R
# Design D1: Fix corr(A, Atilde) and vary f = rho_V / (rho_U + rho_V)

# ---- simulation inputs ----
n_sample <- 1000
n_iters  <- 1000
seed     <- 314

noise_dist <- "norm"
df_t       <- 5
alpha      <- 0.05
robust_se  <- FALSE

out_dir <- "results/withV/D1"

# ---- design inputs (same as current) ----
rho_total_vec <- 0.5
b1_vec        <- 0.3
f_vec         <- seq(0.05, 0.95, by = 0.05)

r_c      <- 0.7
b2       <- 0.3
sigma_eY <- 1

design_label <- "D1_FixRho"

# ---- load project functions ----
source("R/DGM.R")
source("R/FitModels.R")
source("R/MakeGrids.R")
source("R/RunSims.R")

# ---- build grid ----
grid <- MakeGrid_D1_FixRho(
  rho_total_vec = rho_total_vec,
  b1_vec        = b1_vec,
  f_vec         = f_vec,
  r_c           = r_c,
  b2            = b2,
  sigma_eY      = sigma_eY,
  add_derived   = TRUE
)

if (is.null(grid) || nrow(grid) == 0) stop("No feasible parameter combinations in grid.")
grid$design <- design_label

# ---- run sims ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
tag <- sprintf("withV_D1_%s_n%d_it%d_seed%d", noise_dist, n_sample, n_iters, seed)
save_prefix <- file.path(out_dir, paste0(tag, "_", stamp))

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
message(sprintf("[Run_V_D1] Total runtime: %.1f sec (%.2f min)", elapsed_sec, elapsed_sec / 60))

# ---- manifest ----
manifest <- list(
  created_at     = as.character(Sys.time()),
  script         = "Run_V_D1.R",
  design_label   = design_label,
  save_prefix    = save_prefix,
  n_sample       = n_sample,
  n_iters        = n_iters,
  seed           = seed,
  noise_dist     = noise_dist,
  df_t           = df_t,
  alpha          = alpha,
  robust_se      = robust_se,
  rho_total_vec  = rho_total_vec,
  b1_vec         = b1_vec,
  f_vec          = f_vec,
  r_c            = r_c,
  b2             = b2,
  sigma_eY       = sigma_eY,
  grid_dim       = dim(grid),
  grid_preview   = utils::head(grid, 10),
  session_info   = utils::sessionInfo(),
  runtime_minutes = elapsed_sec / 60
)

saveRDS(manifest, file = paste0(save_prefix, "_manifest.rds"))
writeLines(capture.output(str(manifest, max.level = 2)),
           con = paste0(save_prefix, "_manifest.txt"))

message("Saved: ", paste0(save_prefix, "_sims.rds"))
message("Saved: ", paste0(save_prefix, "_agg.csv"))
message("Saved: ", paste0(save_prefix, "_manifest.rds / _manifest.txt"))

invisible(res)
