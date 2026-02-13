# MakeFigure_V_D1.R
# 3-row figure: mean coef / mean SE / power vs f, comparing models with vs without V adjustment

# ---- edit these ----
run_prefix <- "results/withV/D1/withV_D1_norm_n1000_it1000_seed314_20260213_105406"
out_dir    <- "figures"
out_stem   <- "Fig_withV_D1"
# --------------------

library(ggplot2)
library(grid)

# ---- read agg ----
infile <- paste0(run_prefix, "_agg.csv")
agg <- read.csv(infile)

# ---- enforce: this figure is for ONE rho_target only ----
rho_vals <- unique(agg$rho_target)
if (length(rho_vals) != 1) {
  stop(paste0(
    "This figure expects exactly one rho_target in agg, but found: ",
    paste(rho_vals, collapse = ", ")
  ))
}

# ---- compute f for plotting ----
if ("f_target" %in% names(agg)) {
  agg$f_plot <- agg$f_target
} else if (all(c("rho_U", "rho_V") %in% names(agg))) {
  agg$f_plot <- with(agg, rho_V / (rho_U + rho_V))
} else if ("rho_rel_U" %in% names(agg)) {
  agg$f_plot <- 1 - agg$rho_rel_U
} else {
  stop("Cannot compute f: need f_target, or (rho_U & rho_V), or rho_rel_U.")
}

# IMPORTANT for presentation: ensure rows are ordered by f
agg <- agg[order(agg$f_plot), , drop = FALSE]

# ---- required columns ----
need <- c(
  "beta_Atilde_fit2", "beta_Atilde_fit4",
  "se_Atilde_fit2",   "se_Atilde_fit4",
  "power_Atilde_fit2","power_Atilde_fit4"
)
miss <- setdiff(need, names(agg))
if (length(miss) > 0) stop(paste("Agg is missing columns:", paste(miss, collapse = ", ")))

# ---- labels ----
model_levels <- c("No V adjustment", "Adjusted for V")

# ---- reshape helper ----
MakeLong <- function(df, col_noV, col_withV) {
  rbind(
    data.frame(f = df$f_plot, model = model_levels[1], value = df[[col_noV]]),
    data.frame(f = df$f_plot, model = model_levels[2], value = df[[col_withV]])
  )
}

df_coef <- MakeLong(agg, "beta_Atilde_fit2",  "beta_Atilde_fit4")
df_se   <- MakeLong(agg, "se_Atilde_fit2",    "se_Atilde_fit4")
df_pow  <- MakeLong(agg, "power_Atilde_fit2", "power_Atilde_fit4")

# ---- shared styling ----
ScaleModel <- function() {
  list(
    scale_linetype_manual(values = c("solid", "dashed"), breaks = model_levels),
    scale_shape_manual(values = c(16, 17), breaks = model_levels),
    scale_color_manual(values = c("black", "blue"), breaks = model_levels)
  )
}

BaseTheme <- function(legend_pos = "none") {
  theme_minimal(base_size = 12) +
    theme(
      legend.position = legend_pos,
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

MakePlot <- function(df, ylab, title, xlab = NULL, ylim01 = FALSE, legend_pos = "none") {
  p <- ggplot(df, aes(x = f, y = value,
                      linetype = model, color = model, shape = model,
                      group = model)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3) +
    ScaleModel() +
    scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    labs(title = title, x = xlab, y = ylab,
         linetype = NULL, color = NULL, shape = NULL) +
    BaseTheme(legend_pos)
  
  if (ylim01) p <- p + scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))
  p
}

p_coef <- MakePlot(
  df_coef,
  ylab = "Mean coefficient",
  title = "(A) Mean NCE coefficient",
  legend_pos = "none"
)

p_se <- MakePlot(
  df_se,
  ylab = "Mean SE",
  title = "(B) SE of NCE coefficient",
  legend_pos = "none"
) + scale_y_continuous(breaks = seq(0.03, 0.05, 0.005), limits = c(0.03, 0.05))

p_pow <- MakePlot(
  df_pow,
  ylab = "Power",
  title = "(C) Power of NCE test",
  xlab = "Share of corr(A, Ã) attributable to V",
  ylim01 = TRUE,
  legend_pos = "bottom"
)

DrawCombined <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 1)))
  vp <- function(r) viewport(layout.pos.row = r, layout.pos.col = 1)
  
  print(p_coef, vp = vp(1))
  print(p_se,   vp = vp(2))
  print(p_pow,  vp = vp(3))
}

# ---- display + save ----
DrawCombined()

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

png(file.path(out_dir, paste0(out_stem, ".png")), width = 12, height = 8.5, units = "in", res = 300)
DrawCombined(); dev.off()

pdf(file.path(out_dir, paste0(out_stem, ".pdf")), width = 12, height = 8.5)
DrawCombined(); dev.off()
