# Run_V_D2.R
# Runs simulations for Design D2:
#   - Fix bias via fixed a1 (bias in Y~A is b1*a1)
#   - Fix rho_total = cor(A, Atilde)
#   - Vary f in [0,1], where f = rho_V / rho_total

# -------------------- user inputs --------------------
n_sample <- 1000
n_iters  <- 1000
seed     <- 314

noise_dist <- "norm"  # "norm" or "t"
df_t       <- 5
alpha      <- 0.05
robust_se  <- FALSE

# D2 controls
rho_total <- 0.50
f_vec     <- seq(0, 1, by = 0.05)

# Choose a1,a2 so the whole f range is feasible.
# A simple, safe choice for full f in [0,1]:
a1 <- 0.60
a2 <- 0.60

# Outcome model parameters
b0 <- 0
b1 <- 0.30
b2 <- 0.30
sigma_eY <- 1

a0 <- 0
c0 <- 0

out_dir <- "results/withV/D2"
design_label <- "D2"

# -------------------- load project code --------------------
source("R/DGM.R")
source("R/FitModels.R")
source("R/MakeGrids.R")
source("R/RunSims.R")

# -------------------- grid --------------------
grid <- MakeGrid_D2_FixBiasFixRho_VaryF(
  rho_total = rho_total,
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

grid$design <- design_label

# -------------------- run --------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
tag <- sprintf("withV_D2_%s_n%d_it%d_seed%d", noise_dist, n_sample, n_iters, seed)
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
message(sprintf("[Run_V_D2] Total runtime: %.1f sec (%.2f min)", elapsed_sec, elapsed_sec / 60))

# -------------------- manifest --------------------
manifest <- list(
  created_at      = as.character(Sys.time()),
  script          = "Run_V_D2.R",
  design_label    = design_label,
  save_prefix     = save_prefix,
  n_sample        = n_sample,
  n_iters         = n_iters,
  seed            = seed,
  noise_dist      = noise_dist,
  df_t            = df_t,
  alpha           = alpha,
  robust_se       = robust_se,
  rho_total       = rho_total,
  f_vec           = f_vec,
  a1              = a1,
  a2              = a2,
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

message("Saved: ", paste0(save_prefix, "_sims.rds"))
message("Saved: ", paste0(save_prefix, "_agg.csv"))
message("Saved: ", paste0(save_prefix, "_manifest.rds / _manifest.txt"))

invisible(res)