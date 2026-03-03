# Run_V_D1_Sims_WithX.R
# Runs simulations for Design D1 (with V), with two measured confounders.
# - D1 grid varies only V->A and V->Atilde via a2=c2 to hit rho_target

# -------------------- user inputs --------------------
n_sample <- 1000
n_iters  <- 1000
seed     <- 314

noise_dist <- "norm"  # "norm" or "t"
df_t       <- 5
alpha      <- 0.05
robust_se  <- FALSE

# Measured covariates affect A and Y, not Atilde
ax_1 <- 0.15
ax_2 <- 0.15
bx_1 <- 0.2
bx_2 <- 0.2
p_x2 <- 0.5

rho_target_vec <- seq(0.1, 0.9, by = 0.05)

# Set a1*c1 = 0.1 but with c1/a1 = 0.6 so beta_Atilde crosses ~0 at rho=0.6
a1 <- sqrt(0.1 / 0.6)
c1 <- 0.6 * a1

# Outcome model parameters
b0 <- 0
b1 <- 0.3
b2 <- 0.3
sigma_eY <- 1

out_dir <- "results/withX/withV/D1"
design_label <- "D1_WithX"

# -------------------- load project code --------------------
source("R/withX/DGM_WithX.R")
source("R/withX/FitModels_WithX.R")
source("R/withX/MakeGrids_WithX.R")
source("R/withX/RunSims_WithX.R")

# -------------------- output prefix  ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
tag <- sprintf("withV_D1_withX_%s_n%d_it%d_seed%d", noise_dist, n_sample, n_iters, seed)
save_prefix <- file.path(out_dir, paste0(tag, "_", stamp))

# -------------------- grid --------------------
grid <- MakeGrid_D1_FixU_VaryV_WithX(
  rho_target_vec = rho_target_vec,
  a1 = a1, c1 = c1,
  b0 = b0, b1 = b1, b2 = b2,
  ax_1 = ax_1, ax_2 = ax_2, bx_1 = bx_1, bx_2= bx_2, p_x2 = p_x2,
  a0 = 0, c0 = 0,
  sigma_eY = sigma_eY,
  add_derived = TRUE
)
grid$design <- design_label

# -------------------- run sims --------------------
t0 <- proc.time()[3]
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
  save_prefix = save_prefix
)
elapsed_sec <- proc.time()[3] - t0

# -------------------- manifest --------------------
manifest <- list(
  created_at = as.character(Sys.time()),
  script = "Run_V_D1_Sims_WithX.R",
  save_prefix = save_prefix,
  n_sample = n_sample,
  n_iters = n_iters,
  seed = seed,
  noise_dist = noise_dist,
  df_t = df_t,
  alpha = alpha,
  robust_se = robust_se,
  rho_target_vec = rho_target_vec,
  a1 = a1, c1 = c1,
  b0 = b0, b1 = b1, b2 = b2,
  sigma_eY = sigma_eY,
  ax_1 = ax_1, ax_2 = ax_2, bx_1 = bx_1, bx_2 = bx_2, p_x2 = p_x2,
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
