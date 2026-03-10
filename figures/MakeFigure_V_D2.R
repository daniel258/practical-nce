# MakeFigure_V_D2.R
# Two-panel figure for Design D2 (vary f in [0,1]):
#   (A) mean NCE coefficient +/-  SD
#   (B) power for H0: beta_Atilde = 0
# Shown for Model 2 (Y~A+Atilde) and Model 4 (Y~A+Atilde+V)

# ---- file prefix name ----
run_prefix <- "results/withV/D2/withV_D2_norm_n1000_it1000_seed314_20260216_082059"
# reads: paste0(run_prefix, "_agg.csv")

# ---- output ----
out_dir  <- "figures"
out_stem <- "Fig_D2"

library(ggplot2)

# ---- read agg ----
infile <- paste0(run_prefix, "_agg.csv")
agg <- read.csv(infile)

# x-axis: f (prefer f_target, fallback to derived)
if ("f_target" %in% names(agg)) {
  agg$f_plot <- agg$f_target
} else if ("f_achieved" %in% names(agg)) {
  agg$f_plot <- agg$f_achieved
} else if ("rho_rel_U" %in% names(agg)) {
  agg$f_plot <- 1 - agg$rho_rel_U
} else {
  stop("Can't find f_target/f_achieved/rho_rel_U in agg.csv")
}

# sort for clean lines
agg <- agg[order(agg$f_plot), ]

# ---- long format for models 2 and 4 ----
d2 <- data.frame(
  f     = agg$f_plot,
  model = "Y ~ A + Atilde",
  beta  = agg$beta_Atilde_fit2,
  se    = agg$se_Atilde_fit2,
  power = agg$power_Atilde_fit2
)

d4 <- data.frame(
  f     = agg$f_plot,
  model = "Y ~ A + Atilde + V",
  beta  = agg$beta_Atilde_fit4,
  se    = agg$se_Atilde_fit4,
  power = agg$power_Atilde_fit4
)

dat <- rbind(d2, d4)
dat$model <- factor(dat$model, levels = c("Y ~ A + Atilde", "Y ~ A + Atilde + V"))
dat$lo <- dat$beta - dat$se
dat$hi <- dat$beta + dat$se

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

# ---- panel (A): mean coef +/- mean SD (not a CI band) ----
p_coef <- ggplot(dat, aes(x = f, y = beta, color = model, shape = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = 0.18, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.9, linetype = 2) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(
    title = "(A) Mean NCE coefficient \u00B1 SD",
    x = expression(pi[V]),
    # x = expression(paste("Share of Corr(A, ", tilde(A), ") from V (", pi[V],")")), 
    y = expression(paste("NCE coefficient")),# (", hat(beta)[tilde(A)], ")")),
    color = "Model:",
    shape = "Model:"
  ) +
  BaseTheme(legend_pos = "bottom") +
  scale_color_discrete(labels = c(
    expression(Y %~% A + tilde(A)),
    expression(Y %~% A + tilde(A) + V)
  )) +
  scale_shape_discrete(labels = c(
    expression(Y %~% A + tilde(A)),
    expression(Y %~% A + tilde(A) + V)
  )) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

# ---- panel (B): power ----
p_pow <- ggplot(dat, aes(x = f, y = power, color = model, shape = model)) +
  geom_line(linewidth = 0.9, linetype = 2) +
  geom_point(size = 1.6) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "(B) Power",
   # x = expression(paste("Share of correlation from V (", rho[V]/rho, ")")),    
    x = expression(pi[V]),    
    y = "Power",
    color = "Model:",
    shape = "Model:"
  ) +
  BaseTheme(legend_pos = "bottom") +
  scale_color_discrete(labels = c(
    expression(Y %~% A + tilde(A)),
    expression(Y %~% A + tilde(A) + V)
  )) +
  scale_shape_discrete(labels = c(
    expression(Y %~% A + tilde(A)),
    expression(Y %~% A + tilde(A) + V)
  )) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

# ---- one legend, two panels ----
legend_grob <- GetLegendGrob(p_coef)
p_coef_noleg <- p_coef + theme(legend.position = "none")
p_pow_noleg  <- p_pow  + theme(legend.position = "none")

DrawCombined <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = 2, ncol = 2,
    heights = grid::unit(c(1, 0.18), "null")
  )))
  vp <- function(r, c) grid::viewport(layout.pos.row = r, layout.pos.col = c)
  
  print(p_coef_noleg, vp = vp(1, 1))
  print(p_pow_noleg,  vp = vp(1, 2))
  
  if (!is.null(legend_grob)) {
    grid::grid.draw(grid::editGrob(legend_grob, vp = vp(2, 1:2)))
  }
}

# draw to R device
DrawCombined()

# ---- save PNG + PDF ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# png(file.path(out_dir, paste0(out_stem, ".png")), width = 8.5, height = 4.5, units = "in", res = 300)
# DrawCombined()
# dev.off()

pdf(file.path(out_dir, paste0(out_stem, ".pdf")), width = 8.5, height = 4.5)
DrawCombined()
dev.off()
