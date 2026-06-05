# Run_V_D1_Sims.R
# Runs simulations for Design D1 (with V), where:
# - D1 grid varies only V->A and V->Atilde via a2=c2 to hit rho_target

# -------------------- user inputs --------------------
n_sample <- 1000
n_iters  <- 1000
seed     <- 314

noise_dist <- "norm"  # "norm" or "t"
df_t       <- 3
alpha      <- 0.05
robust_se  <- FALSE

rho_target_vec <- seq(0.1, 0.9, by = 0.05)

# Set a1*c1 = 0.1 but with c1/a1 = 0.6 so beta_Atilde crosses ~0 at rho=0.6
a1 <- sqrt(0.1 / 0.6)
c1 <- 0.6 * a1       

# Fix everything else
a0 <- 0
c0 <- 0

b0 <- 0
b1 <- 0.3
b2 <- 0.3
sigma_eY <- 1

# Output folder + prefix
out_dir <- "results/withV/D1"
design_label <- "D1"

# -------------------- load project code --------------------
source("R/DGM.R")
source("R/FitModels.R")
source("R/MakeGrids.R")
source("R/RunSims.R")

# -------------------- output prefix  ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
tag <- sprintf("withV_D1_%s_n%d_it%d_seed%d", noise_dist, n_sample, n_iters, seed)
save_prefix <- file.path(out_dir, paste0(tag, "_", stamp))

# -------------------- grid --------------------
grid <- MakeGrid_D1_FixU_VaryV(
  rho_target_vec = rho_target_vec,
  a1 = a1, c1 = c1,
  a0 = a0, c0 = c0,
  b0 = b0, b1 = b1, b2 = b2,
  sigma_eY = sigma_eY
)

grid$design <- design_label

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
message(sprintf("[Run_V_D1_Sims] Total runtime: %.1f sec (%.2f min)", elapsed_sec, elapsed_sec / 60))

# -------------------- manifest --------------------
manifest <- list(
  created_at      = as.character(Sys.time()),
  script          = "Run_V_D1_Sims.R",
  design_label    = design_label,
  save_prefix     = save_prefix,
  n_sample        = n_sample,
  n_iters         = n_iters,
  seed            = seed,
  noise_dist      = noise_dist,
  df_t            = df_t,
  alpha           = alpha,
  robust_se       = robust_se,
  rho_target_vec  = rho_target_vec,
  a1              = a1,
  c1              = c1,
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
