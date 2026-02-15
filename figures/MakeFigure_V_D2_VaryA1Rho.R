# MakeFigure_V_D2_VaryA1Rho.R
# Creates two figures (separate files):
#   1) Mean NCE coefficient +/- mean SE
#   2) Power
# Facets: rows = a1, cols = rho_total
# Handles infeasible f by inserting NA rows so lines/ribbons break.

# ---- YOU EDIT THIS ----
run_prefix <- "results/withV/D2/withV_D2_varyA1Rho_norm_n1000_it1000_seed314_20260213_172543"
out_dir    <- "figures"
out_stem   <- "Fig_D2_varyA1Rho"
# -----------------------

library(ggplot2)

infile <- paste0(run_prefix, "_agg.csv")
agg <- read.csv(infile)

# x-axis f: prefer f_target, else f_achieved, else 1-rho_rel_U
if ("f_target" %in% names(agg)) {
  agg$f_plot <- agg$f_target
} else if ("f_achieved" %in% names(agg)) {
  agg$f_plot <- agg$f_achieved
} else if ("rho_rel_U" %in% names(agg)) {
  agg$f_plot <- 1 - agg$rho_rel_U
} else {
  stop("Can't find f_target/f_achieved/rho_rel_U in agg.csv")
}

stopifnot("rho_target" %in% names(agg))
stopifnot("a1" %in% names(agg))

# Recover full intended f grid from manifest if present
f_full <- NULL
manifest_file <- paste0(run_prefix, "_manifest.rds")
if (file.exists(manifest_file)) {
  man <- readRDS(manifest_file)
  if (!is.null(man$f_vec)) f_full <- as.numeric(man$f_vec)
}
if (is.null(f_full)) f_full <- seq(0, 1, by = 0.05)

# Facet labels (formatted factors)
rho_levels <- sort(unique(agg$rho_target))
a1_levels  <- sort(unique(agg$a1))

agg$rho_facet <- factor(sprintf("%.2f", agg$rho_target),
                        levels = sprintf("%.2f", rho_levels))
agg$a1_facet  <- factor(sprintf("%.2f", agg$a1),
                        levels = sprintf("%.2f", a1_levels))

# long format: models 2 and 4 only
d2 <- data.frame(
  rho   = agg$rho_facet,
  a1    = agg$a1_facet,
  f     = agg$f_plot,
  model = "Y ~ A + Atilde",
  beta  = agg$beta_Atilde_fit2,
  se    = agg$se_Atilde_fit2,
  power = agg$power_Atilde_fit2
)

d4 <- data.frame(
  rho   = agg$rho_facet,
  a1    = agg$a1_facet,
  f     = agg$f_plot,
  model = "Y ~ A + Atilde + V",
  beta  = agg$beta_Atilde_fit4,
  se    = agg$se_Atilde_fit4,
  power = agg$power_Atilde_fit4
)

dat <- rbind(d2, d4)
dat$model <- factor(dat$model, levels = c("Y ~ A + Atilde", "Y ~ A + Atilde + V"))

# Insert NA rows on the full f grid so lines/ribbons break where f is infeasible
support <- expand.grid(
  rho   = levels(dat$rho),
  a1    = levels(dat$a1),
  f     = f_full,
  model = levels(dat$model),
  stringsAsFactors = FALSE
)

dat2 <- merge(support, dat, by = c("rho", "a1", "f", "model"), all.x = TRUE, sort = FALSE)
dat2$model <- factor(dat2$model, levels = levels(dat$model))
dat2$lo <- dat2$beta - dat2$se
dat2$hi <- dat2$beta + dat2$se

# Order for clean plotting
dat2 <- dat2[order(dat2$a1, dat2$rho, dat2$model, dat2$f), ]

BaseTheme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

FacetLabeller <- labeller(
  rho = function(x) paste0("rho = ", x),
  a1  = function(x) paste0("a1 = ", x)
)

# -------- Figure 1: NCE coef +/- SE --------
p_nce <- ggplot(dat2, aes(x = f, y = beta, color = model, shape = model, group = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model, group = model),
              alpha = 0.18, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.9, linetype = 2, na.rm = FALSE) +
  geom_point(size = 1.5, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_grid(a1 ~ rho, labeller = FacetLabeller) +
  labs(
    title = "Mean NCE coefficient \u00B1 SE",
    x = expression(rho[V] / rho),
    y = expression(paste("NCE coefficient (", hat(beta)[tilde(A)], ")")),
    color = "Model:",
    shape = "Model:"
  ) +
  BaseTheme() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

# -------- Figure 2: Power --------
p_pow <- ggplot(dat2, aes(x = f, y = power, color = model, shape = model, group = model)) +
  geom_line(linewidth = 0.9, linetype = 2, na.rm = FALSE) +
  geom_point(size = 1.5, na.rm = TRUE) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_grid(a1 ~ rho, labeller = FacetLabeller) +
  labs(
    title = "Power",
    x = expression(rho[V] / rho),
    y = "Power",
    color = "Model:",
    shape = "Model:"
  ) +
  BaseTheme() +
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

png(file.path(out_dir, paste0(out_stem, "_NCE.png")), width = width_in, height = height_in, units = "in", res = 300)
print(p_nce)
dev.off()

pdf(file.path(out_dir, paste0(out_stem, "_NCE.pdf")), width = width_in, height = height_in)
print(p_nce)
dev.off()

png(file.path(out_dir, paste0(out_stem, "_Power.png")), width = width_in, height = height_in, units = "in", res = 300)
print(p_pow)
dev.off()

pdf(file.path(out_dir, paste0(out_stem, "_Power.pdf")), width = width_in, height = height_in)
print(p_pow)
dev.off()
