# MakeFigure_NoV_NaiveBias_OneLegend.R
# Requires: ggplot2

# ---- YOU EDIT THIS ----
run_prefix <- "results/noV/noV_norm_n1000_it1000_seed314_20260206_114532"
# Reads: paste0(run_prefix, "_agg.csv")

# ---- output ----
out_dir  <- "figures"
out_stem <- "Fig_NoV"

library(ggplot2)

# ---- read agg ----
infile <- paste0(run_prefix, "_agg.csv")
agg <- read.csv(infile)

agg$c1_f <- factor(agg$c1)

# ---- theme ----
BaseTheme <- function(legend_pos = "none") {
  theme_minimal(base_size = 12) +
    theme(
      legend.position = legend_pos,
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# ---- extract legend grob ----
GetLegendGrob <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx[1]]]
}

# ---- panels ----
p_bias <- ggplot(agg, aes(x = a1, y = bias_A_fit1, color = c1_f)) +
  geom_line(linewidth = 0.8, linetype = 2) +
  geom_point(size = 1.6) +
  labs(title = "(A) Bias in Exposure Effect", x = NULL, y = "Bias", color = expression(c[1])) +
  BaseTheme(legend_pos = "bottom")

p_coef <- ggplot(agg, aes(x = a1, y = beta_Atilde_fit2, color = c1_f)) +
  geom_line(linewidth = 0.8, linetype = 2) +
  geom_point(size = 1.6) +
  labs(title = "(B) NCE Coefficient", x = NULL, y = "Coefficient", color = expression(c[1])) +
  BaseTheme()

p_pow <- ggplot(agg, aes(x = a1, y = power_Atilde_fit2, color = c1_f)) +
  geom_line(linewidth = 0.8, linetype = 2) +
  geom_point(size = 1.6) +
  labs(title = "(C) Power of NCE Test", x = expression(a[1]), y = "Power", color = expression(c[1])) +
  BaseTheme()

p_se <- ggplot(agg, aes(x = a1, y = se_Atilde_fit2, color = c1_f)) +
  geom_line(linewidth = 0.8, linetype = 2) +
  geom_point(size = 1.6) +
  labs(title = "(D) SE of NCE Coefficient", x = expression(a[1]), y = "SE", color = expression(c[1])) +
  BaseTheme()

legend_grob <- GetLegendGrob(p_bias)
p_bias_noleg <- p_bias + theme(legend.position = "none")

# ---- draw combined 2x2 + legend to current device ----
DrawCombined <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = 3, ncol = 2,
    heights = grid::unit(c(1, 1, 0.20), "null")
  )))
  vp <- function(r, c) grid::viewport(layout.pos.row = r, layout.pos.col = c)
  
  print(p_bias_noleg, vp = vp(1, 1))
  print(p_coef,       vp = vp(1, 2))
  print(p_pow,        vp = vp(2, 1))
  print(p_se,         vp = vp(2, 2))
  
  if (!is.null(legend_grob)) {
    grid::grid.draw(grid::editGrob(legend_grob, vp = vp(3, 1:2)))
  }
}

DrawCombined()

# ---- save PNG + PDF ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

png(file.path(out_dir, paste0(out_stem, ".png")), width = 8.5, height = 6.5, units = "in", res = 300)
DrawCombined()
dev.off()

pdf(file.path(out_dir, paste0(out_stem, ".pdf")), width = 8.5, height = 6.5)
DrawCombined()
dev.off()
