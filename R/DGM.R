# dgm.R
# Data-generating mechanisms (DGMs) for Practical Considerations NCEs
#
# DGM
#   U ~ N(0,1),  V ~ N(0,1), independent
#   A      = a0 + a1*U + a2*V + eA
#   Y      = b0 + b1*U + b2*A + eY
#   Atilde = c0 + c1*U + c2*V + eAt
#
# Var(eA) and Var(eAt) are chosen so that Var(A)=Var(Atilde)=1.
# For first DGM a2 = c2 = 0
# This is achieved by setting Var(eA) = 1 - a1^2 - a2^2 and Var(eAt) = 1 - c1^2 - c2^2.

## The main function is DGM() below. The rest of the functions are internal helpers.
# Compute residual SD to enforce Var(linear combo + residual)=1
CalcResidSD <- function(a1, a2, c1, c2) {
  v_eA  <- 1 - (a1^2 + a2^2)
  v_eAt <- 1 - (c1^2 + c2^2)
  if (v_eA <= 0 | v_eAt <= 0 ) {
    stop(sprintf(
      "Residual variances must be positive. Computed Var(eA)=%.6f and Var(eAt)=%.6f.",
      v_eA, v_eAt
    ))
  }
  c(sqrt(v_eA), sqrt(v_eAt))
}

# Sample mean-zero, and given variance level noise with optional heavy tails
#' dist: "norm" or "t"
rnoise <- function(n_sample, dist = c("norm", "t"), df = 5, sigma) {
  dist <- match.arg(dist)
  if (dist == "norm") 
  noise <- rnorm(n_sample, mean = 0, sd = sigma)
  else
  {noise_temp <- rt(n_sample, df = df)
  noise <- sigma * noise_temp / sqrt(df / (df - 2))}
  return(noise)
}

# I create a general DGM function for the extended model, and then the no V-case is obtained by setting a2=c2=0. 
# pars is a named list with elements a0,a1,a2,b0,b1,b2,c0,c1,c2,sigma_eY
DGM <- function(n_sample, pars, noise_dist  = c("norm", "t"), df = 5) {
  noise_dist  <- match.arg(noise_dist)
  a0 <- pars$a0; a1  <- pars$a1;  a2 <- pars$a2
  b0 <- pars$b0; b1  <- pars$b1;  b2 <- pars$b2
  c0 <- pars$c0; c1  <- pars$c1;  c2 <- pars$c2
  sigma_eY <- pars$sigma_eY
  sd_resids <- CalcResidSD(a1, a2, c1, c2)
  sigma_eA <- sd_resids[1]
  sigma_eAt <- sd_resids[2]
  eA  <- rnoise(n_sample, dist = noise_dist, df = df, sigma = sigma_eA)  
  eAt <- rnoise(n_sample, dist = noise_dist, df = df, sigma = sigma_eAt)
  eY  <- rnoise(n_sample, dist = noise_dist, df = df, sigma = sigma_eY) 
  U <- rnorm(n_sample)
  V <- rnorm(n_sample)
  A <- a0 + a1 * U + a2 * V + eA
  Y <- b0 + b1 * U + b2 * A + eY
  Atilde <- c0 + c1 * U + c2 * V + eAt
  
  out <- data.frame(Y = Y, A = A, Atilde = Atilde, U, V) # return U,V for checking purposes only
  out
}

#### Here I will add a version with covariates

#####

# ---- quick tests ----
# pars0 <- list(a0=0, a1=0.5, a2=0, b0=0, b1=0.3, b2=0.4, c0=0, c1=0.6, c2=0, sigma_eY=1)
# dat0 <-DGM(n_sample=1e6, pars=pars0, noise_dist="norm")
# ff <- (lm(formula = Y ~ A + U + V + Atilde, data = dat0))
# var(dat0$A)      
# var(dat0$Atilde)

