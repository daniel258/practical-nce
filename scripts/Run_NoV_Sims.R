# Run_NoV_Sims.R
# Run simulations for the "no V affects A/Atilde" case by setting a2=c2=0 in the grid.
# Uses your existing DGM(), FitModels(), RunSims() without modification.

# ---- user inputs ----

n_sample <- 1000
n_iters  <- 1000
seed     <- 314

noise_dist <- "norm"  # "norm" or "t"
df_t       <- 5
alpha      <- 0.05
robust_se  <- FALSE

# Grid choices (noV grid varies a1 and c1)
a1_vec <- seq(0.05, 0.95, by = 0.05)
c1_vec <- c(0.1, 0.3, 0.5, 0.7, 0.9)

# Structural parameters
b1 <- 0.5
b2 <- 0.3
sigma_eY <- 1

# Output folder + prefix
out_dir <- "results/noV"
tag <- sprintf("noV_%s_n%d_it%d_seed%d", noise_dist, n_sample, n_iters, seed)

# ---- load project functions ----

source("R/DGM.R")
source("R/FitModels.R")
source("R/MakeGrids.R")
source("R/RunSims.R")

# ---- build NO-V grid  ----
#  MakeGrid_NoV() hard-sets a2=0 and c2=0.
grid <- MakeGrid_NoV(
  a1 = a1_vec,
  c1 = c1_vec,
  b1 = b1,
  b2 = b2,
  sigma_eY = sigma_eY,
  add_derived = TRUE
)

# Make output prefix
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save_prefix <- file.path(out_dir, paste0(tag, "_", stamp))

# ---- run sims ----
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
message(sprintf("Total runtime: %.1f seconds (%.2f minutes)", elapsed_sec, elapsed_sec / 60))

# ---- manifest (for reproducibility) ----
manifest <- list(
  created_at = as.character(Sys.time()),
  script = "Run_NoV_Sims.R",
  save_prefix = save_prefix,
  n_sample = n_sample,
  n_iters = n_iters,
  seed = seed,
  noise_dist = noise_dist,
  df_t = df_t,
  alpha = alpha,
  robust_se = robust_se,
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
