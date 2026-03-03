# Run_NoV_Sims_WithX.R
# Run simulations for the "no V affects A/Atilde" case (a2=c2=0),
# with two measured covariates (X1 continuous, X2 binary).

# -------------------- user inputs --------------------

n_sample <- 1000
n_iters  <- 1000
seed     <- 314

noise_dist <- "norm"  # "norm" or "t"
df_t       <- 5
alpha      <- 0.05
robust_se  <- FALSE

# Measured covariates: affect A and Y, not Atilde
ax_1 <- 0.15
ax_2 <- 0.15
bx_1 <- 0.2
bx_2 <- 0.2
p_x2 <- 0.5

# Grid choices (noV grid varies a1 and c1)
a1_vec <- seq(0.05, 0.95, by = 0.05)
c1_vec <- c(0.1, 0.3, 0.5, 0.7, 0.9)

# Structural parameters
b0 <- 0
b1 <- 0.5
b2 <- 0.3
sigma_eY <- 1

# Output folder + prefix
out_dir <- "results/withX/noV"
design_label <- "NoV_WithX"

# -------------------- load project code --------------------
source("R/withX/DGM_WithX.R")
source("R/withX/FitModels_WithX.R")
source("R/withX/MakeGrids_WithX.R")
source("R/withX/RunSims_WithX.R")

# -------------------- output prefix  ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
tag <- sprintf("noV_withX_%s_n%d_it%d_seed%d", noise_dist, n_sample, n_iters, seed)
save_prefix <- file.path(out_dir, paste0(tag, "_", stamp))

# ---- build grid ----
grid <- MakeGrid_NoV_WithX(
  a1 = a1_vec,
  c1 = c1_vec,
  b1 = b1,
  b2 = b2,
  ax_1 = ax_1, ax_2 = ax_2,
  bx_1 = bx_1, bx_2 = bx_2,
  p_x2 = p_x2,
  a0 = 0, b0 = b0, c0 = 0,
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
  script = "Run_NoV_Sims_WithX.R",
  save_prefix = save_prefix,
  n_sample = n_sample,
  n_iters = n_iters,
  seed = seed,
  noise_dist = noise_dist,
  df_t = df_t,
  alpha = alpha,
  robust_se = robust_se,
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
