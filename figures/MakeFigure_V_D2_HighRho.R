# MakeFigure_V_D2_HighRho.R
# Design D2 (high rho_total), faceted by a1 (rows) and rho_total (cols).
# Outputs two files:
#   1) Mean NCE coefficient +/-  SD
#   2) Power (reject proportion) for H0: beta_Atilde = 0
# Both shown for Model 2 (Y~A+Atilde) and Model 4 (Y~A+Atilde+V)
#
# Notes:
# - The high-rho run script drops infeasible f values from the grid.
#   This script reconstructs the *full* intended f grid (from the manifest
#   when available) and inserts NA rows so lines/ribbons break where infeasible.

# ---- YOU EDIT THIS ----
run_prefix <- "results/withV/D2/withV_D2_highRho_unequalA2_norm_n1000_it1000_seed314_20260216_110119"

out_dir  <- "figures"
out_stem <- "Fig_D2_highRho"
# -----------------------

library(ggplot2)

# ---- read agg ----
infile <- paste0(run_prefix, "_agg.csv")
if (!file.exists(infile)) stop("Cannot find: ", infile)
agg <- read.csv(infile)

# ---- read manifest (preferred: gives full intended grids incl. skipped combos) ----
manifest_file <- paste0(run_prefix, "_manifest.rds")
man <- NULL
if (file.exists(manifest_file)) man <- readRDS(manifest_file)

# ---- x-axis: f (prefer f_target, fallback to derived) ----
if ("f_target" %in% names(agg)) {
  agg$f_plot <- agg$f_target
} else if ("f_achieved" %in% names(agg)) {
  agg$f_plot <- agg$f_achieved
} else if ("rho_rel_U" %in% names(agg)) {
  agg$f_plot <- 1 - agg$rho_rel_U
} else {
  stop("Can't find f_target/f_achieved/rho_rel_U in agg.csv")
}

# ---- rho value (prefer rho_target, fallback to rho_total) ----
if ("rho_target" %in% names(agg)) {
  agg$rho_plot <- agg$rho_target
} else if ("rho_total" %in% names(agg)) {
  agg$rho_plot <- agg$rho_total
} else {
  stop("Can't find rho_target/rho_total in agg.csv")
}

# ---- required columns for faceting ----
if (!("a1" %in% names(agg))) stop("agg.csv is missing column: a1")
if (!("a2" %in% names(agg))) stop("agg.csv is missing column: a2")

# ---- reconstruct intended grids (use manifest when available) ----
F_full   <- if (!is.null(man) && !is.null(man$f_vec))         as.numeric(man$f_vec) else seq(0, 1, by = 0.05)
A1_full  <- if (!is.null(man) && !is.null(man$a1_vec))        as.numeric(man$a1_vec) else sort(unique(agg$a1))
Rho_full <- if (!is.null(man) && !is.null(man$rho_total_vec)) as.numeric(man$rho_total_vec) else sort(unique(agg$rho_plot))

# ---- a1 labels: include a2 when available; try to recover from skipped log too ----
GetA1Map <- function(run_prefix, agg, A1_full) {
  # Start from agg
  map0 <- unique(agg[, c("a1", "a2")])
  map0 <- map0[order(map0$a1), , drop = FALSE]
  
  # Augment with skipped log if present (helps retain a2 for skipped rho panels)
  skipped_file <- paste0(run_prefix, "_skipped.csv")
  if (file.exists(skipped_file)) {
    sk <- read.csv(skipped_file)
    if (all(c("a1", "a2") %in% names(sk))) {
      map1 <- unique(sk[!is.na(sk$a1) & !is.na(sk$a2), c("a1", "a2")])
      map0 <- unique(rbind(map0, map1))
      map0 <- map0[order(map0$a1), , drop = FALSE]
    }
  }
  
  # Ensure one row per a1 (take the first; they should all agree)
  map0 <- map0[!duplicated(map0$a1), , drop = FALSE]
  
  # Force all intended a1's to appear
  map_full <- data.frame(a1 = A1_full, stringsAsFactors = FALSE)
  map_full <- merge(map_full, map0, by = "a1", all.x = TRUE, sort = FALSE)
  
  map_full$a1_label <- ifelse(
    is.na(map_full$a2),
    paste0("a[1]==", sprintf("%.2f", map_full$a1)),
    paste0("a[1]==", sprintf("%.2f", map_full$a1),
           "*','~a[2]==", sprintf("%.2f", map_full$a2))
  )
  map_full
}

a1_map <- GetA1Map(run_prefix, agg, A1_full)

# ---- facet labels ----
rho_levels  <- sprintf("%.2f", Rho_full)
#rho_labels <- paste0("Corr(A, ", tilde(A), ")=", rho_levels)
rho_labels  <- paste0("Corr(A,~tilde(A))==", rho_levels)
# Attach labels to agg
agg <- merge(agg, a1_map[, c("a1", "a1_label")], by = "a1", all.x = TRUE, sort = FALSE)
#agg$rho_label <- paste0("rho=", sprintf("%.2f", agg$rho_plot))
agg$rho_label <- paste0("Corr(A,~tilde(A))==", sprintf("%.2f", agg$rho_plot))

agg$a1_facet  <- factor(agg$a1_label, levels = a1_map$a1_label)
agg$rho_facet <- factor(agg$rho_label, levels = rho_labels)

# ---- long format: models 2 and 4 only ----
d2 <- data.frame(
  a1    = agg$a1_facet,
  rho   = agg$rho_facet,
  f     = agg$f_plot,
  model = "Y ~ A + Atilde",
  beta  = agg$beta_Atilde_fit2,
  se    = agg$se_Atilde_fit2,
  power = agg$power_Atilde_fit2
)

d4 <- data.frame(
  a1    = agg$a1_facet,
  rho   = agg$rho_facet,
  f     = agg$f_plot,
  model = "Y ~ A + Atilde + V",
  beta  = agg$beta_Atilde_fit4,
  se    = agg$se_Atilde_fit4,
  power = agg$power_Atilde_fit4
)

dat <- rbind(d2, d4)
dat$model <- factor(dat$model, levels = c("Y ~ A + Atilde", "Y ~ A + Atilde + V"))

# ---- insert NA rows so lines/ribbons break where f is infeasible ----
support <- expand.grid(
  a1    = levels(dat$a1),
  rho   = levels(dat$rho),
  f     = F_full,
  model = levels(dat$model),
  stringsAsFactors = FALSE
)

dat2 <- merge(support, dat, by = c("a1", "rho", "f", "model"), all.x = TRUE, sort = FALSE)
dat2$model <- factor(dat2$model, levels = levels(dat$model))
dat2$lo <- dat2$beta - dat2$se
dat2$hi <- dat2$beta + dat2$se

# Order for clean plotting
ord <- order(dat2$a1, dat2$rho, dat2$model, dat2$f)
dat2 <- dat2[ord, ]

BaseTheme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

XlabExpr <-  expression(pi[V])  


# -------- NCE coef +/- SD --------
p_nce <- ggplot(dat2, aes(x = f, y = beta, color = model, shape = model, group = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model, group = model),
              alpha = 0.18, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.9, linetype = 2, na.rm = FALSE) +
  geom_point(size = 1.5, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_grid(a1 ~ rho, labeller = labeller(a1 = label_parsed, rho = label_parsed)) +  labs(
    title = "Mean NCE coefficient \u00B1 SD",
    x = XlabExpr,
    y = expression(paste("NCE coefficient (", hat(beta)[tilde(A)], ")")),
    color = "Model:",
    shape = "Model:"
  ) +
  BaseTheme() +
  scale_color_discrete(labels = c(
    expression(Y %~% A + tilde(A)),
    expression(Y %~% A + tilde(A) + V)
  )) +
  scale_shape_discrete(labels = c(
    expression(Y %~% A + tilde(A)),
    expression(Y %~% A + tilde(A) + V)
  )) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

# -------- Power --------
p_pow <- ggplot(dat2, aes(x = f, y = power, color = model, shape = model, group = model)) +
  geom_line(linewidth = 0.9, linetype = 2, na.rm = FALSE) +
  geom_point(size = 1.5, na.rm = TRUE) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_grid(a1 ~ rho, labeller = labeller(a1 = label_parsed, rho = label_parsed)) +
  labs(
    title = "Power",
    x = XlabExpr,
    y = "Power",
    color = "Model:",
    shape = "Model:"
  ) +
  BaseTheme() +
  scale_color_discrete(labels = c(
    expression(Y %~% A + tilde(A)),
    expression(Y %~% A + tilde(A) + V)
  )) +
  scale_shape_discrete(labels = c(
    expression(Y %~% A + tilde(A)),
    expression(Y %~% A + tilde(A) + V)
  )) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

# ---- draw to device ----
print(p_nce)
print(p_pow)

# ---- save (size scales with facets) ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_col <- length(levels(dat2$rho))
n_row <- length(levels(dat2$a1))
width_in  <- max(9, 2.4 * n_col + 3.0)
height_in <- max(6, 2.0 * n_row + 2.5)

# png(file.path(out_dir, paste0(out_stem, "_NCE.png")), width = width_in, height = height_in, units = "in", res = 300)
# print(p_nce)
# dev.off()

pdf(file.path(out_dir, paste0(out_stem, "_NCE.pdf")), width = width_in, height = height_in)
print(p_nce)
dev.off()

# png(file.path(out_dir, paste0(out_stem, "_Power.png")), width = width_in, height = height_in, units = "in", res = 300)
# print(p_pow)
# dev.off()

pdf(file.path(out_dir, paste0(out_stem, "_Power.pdf")), width = width_in, height = height_in)
print(p_pow)
dev.off()