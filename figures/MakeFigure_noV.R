# MakeFigure_NoV.R
# Requires: ggplot2

# ---- file prefix name ----
run_prefix <- "results/noV/noV_norm_n1000_it1000_seed314_20260216_080300"
# Reads below: paste0(run_prefix, "_agg.csv")

# ---- output ----
out_dir  <- "figures"
out_stem <- "Fig_NoV"

library(ggplot2)

# ---- read agg ----
infile <- paste0(run_prefix, "_agg.csv")
agg <- read.csv(infile)

agg$c1_f <- factor(agg$c1, levels = sort(unique(agg$c1)))

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
p_bias <- ggplot(agg, aes(x = a1, y = bias_A_fit1, group = c1_f)) +
  geom_line(linewidth = 0.8, linetype = 2) +
  geom_point(size = 1.6) +
  labs(title = "(A) Bias in Exposure Effect", x = NULL, y = "Bias", color = "black") +
  BaseTheme() + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

p_coef <- ggplot(agg, aes(x = a1, y = beta_Atilde_fit2, color = c1_f)) +
  geom_line(linewidth = 0.8, linetype = 2) +
  geom_point(size = 1.6) +
  labs(title = "(B) NCE Coefficient", x = NULL, y = "NCE Coefficient", color = expression(c[1])) +
  BaseTheme(legend_pos = "bottom") + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

p_pow <- ggplot(agg, aes(x = a1, y = power_Atilde_fit2, color = c1_f)) +
  geom_line(linewidth = 0.8, linetype = 2) +
  geom_point(size = 1.6) +
  labs(title = "(C) Power of NCE Test", x = expression(a[1]), y = "Power", color = expression(c[1])) +
  BaseTheme() + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

p_sd <- ggplot(agg, aes(x = a1, y = sd_beta_Atilde_fit2, color = c1_f)) +
  geom_line(linewidth = 0.8, linetype = 2) +
  geom_point(size = 1.6) +
  labs(title = "(D) SD of NCE Coefficient", x = expression(a[1]), y = "SD", color = expression(c[1])) +
  BaseTheme() + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

legend_grob <- GetLegendGrob(p_coef)
p_coef_noleg <- p_coef + theme(legend.position = "none")

# ---- draw combined 2x2 + legend to current device ----
DrawCombined <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = 3, ncol = 2,
    heights = grid::unit(c(1, 1, 0.20), "null")
  )))
  vp <- function(r, c) grid::viewport(layout.pos.row = r, layout.pos.col = c)
  
  print(p_bias, vp = vp(1, 1))
  print(p_coef_noleg,       vp = vp(1, 2))
  print(p_pow,        vp = vp(2, 1))
  print(p_sd,         vp = vp(2, 2))
  
  if (!is.null(legend_grob)) {
    grid::grid.draw(grid::editGrob(legend_grob, vp = vp(3, 1:2)))
  }
}

DrawCombined()

# ---- save PNG + PDF ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# png(file.path(out_dir, paste0(out_stem, ".png")), width = 8.5, height = 6.5, units = "in", res = 300)
# DrawCombined()
# dev.off()

pdf(file.path(out_dir, paste0(out_stem, ".pdf")), width = 8.5, height = 6.5)
DrawCombined()
dev.off()
