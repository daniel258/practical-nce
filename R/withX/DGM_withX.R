# DGM_WithX.R
# DGMs for Practical Considerations NCEs
# Extension: two measured confounders X1 (continuous) and X2 (scaled binary),
# and include them in the outcome and exposure models. 
#   - X affects A and Y (confounding that is controlled by adjustment)
#   - X does not affect N (so Corr(A, N) is still driven by U and V only)
#
#  DGM
#   U ~ N(0,1),  V ~ N(0,1), independent
#   A      = a0 + a1*U + a2*V + ax_1*X1 + ax_2*X2 + eA
#   Y      = b0 + b1*U + b2*A + bx_1*X1 + bx_2*X2 + eY
#   N = c0 + c1*U + c2*V + eN
#
# Var(eA) and Var(eN) are chosen so that Var(A)=Var(N)=1.
# Because X enters A, feasibility requires:
#   a1^2 + a2^2 + ax_1^2 + ax_2^2 < 1  and  c1^2 + c2^2 < 1

# Standardize a Bernoulli(p) to mean 0, var 1
GenBernStd <- function(n_sample, p = 0.5) {
  if (p <= 0 || p >= 1) stop("p must be in (0,1).")
  x_raw <- rbinom(n_sample, size = 1, prob = p)
  (x_raw - p) / sqrt(p * (1 - p))
}

# Compute residual SD to enforce Var=1
CalcResidSD_WithX <- function(a1, a2, c1, c2, ax_1, ax_2) {
  v_eA  <- 1 - (a1^2 + a2^2 + ax_1^2 + ax_2^2)
  v_eN <- 1 - (c1^2 + c2^2)
  if (v_eA <= 0 | v_eN <= 0) {
    stop(sprintf(
      "Residual variances must be positive. Computed Var(eA)=%.6f and Var(eN)=%.6f.",
      v_eA, v_eN
    ))
  }
  c(sqrt(v_eA), sqrt(v_eN))
}

# Sample mean-zero noise with chosen variance, optionally heavy tails
rnoise <- function(n_sample, dist = c("norm", "t"), df = 5, sigma) {
  dist <- match.arg(dist)
  if (dist == "norm") {
    noise <- rnorm(n_sample, mean = 0, sd = sigma)
  } else {
    if (df <= 2) stop("For t noise, df must be > 2.")
    noise_temp <- rt(n_sample, df = df)
    # scale to have Var = sigma^2
    noise <- sigma * noise_temp / sqrt(df / (df - 2))
  }
  noise
}

# Main DGM function with measured covariates X1 and X2
# pars is a named list with elements:
#   a0,a1,a2,  b0,b1,b2,  c0,c1,c2,  sigma_eY,
#   ax_1,ax_2, bx_1,bx_2, p_x2
DGM_WithX <- function(n_sample, pars, noise_dist = c("norm", "t"), df = 5) {
  noise_dist <- match.arg(noise_dist)
  
  # core
  a0 <- pars$a0; a1 <- pars$a1; a2 <- pars$a2
  b0 <- pars$b0; b1 <- pars$b1; b2 <- pars$b2
  c0 <- pars$c0; c1 <- pars$c1; c2 <- pars$c2
  sigma_eY <- pars$sigma_eY
  
  # measured X
  ax_1 <- pars$ax_1; ax_2 <- pars$ax_2
  bx_1 <- pars$bx_1; bx_2 <- pars$bx_2
  p_x2 <- pars$p_x2

  sd_resids <- CalcResidSD_WithX(a1, a2, c1, c2, ax_1, ax_2)
  sigma_eA  <- sd_resids[1]
  sigma_eN <- sd_resids[2]
  
  # noises
  eA  <- rnoise(n_sample, dist = noise_dist, df = df, sigma = sigma_eA)
  eN <- rnoise(n_sample, dist = noise_dist, df = df, sigma = sigma_eN)
  eY  <- rnoise(n_sample, dist = noise_dist, df = df, sigma = sigma_eY)
  
  # latent + measured drivers
  U  <- rnorm(n_sample)
  V  <- rnorm(n_sample)
  X1 <- rnorm(n_sample)
  X2 <- GenBernStd(n_sample, p = p_x2)
  
  A <- a0 + a1 * U + a2 * V + ax_1 * X1 + ax_2 * X2 + eA
  N <- c0 + c1 * U + c2 * V + eN
  Y <- b0 + b1 * U + b2 * A + bx_1 * X1 + bx_2 * X2 + eY
  
  data.frame(Y = Y, A = A, N = N, U = U, V = V, X1 = X1, X2 = X2)
}
