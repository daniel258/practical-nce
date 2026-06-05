# dgm.R
# Data-generating mechanisms (DGMs) for Practical Considerations NCEs
#
# DGM
#   U ~ N(0,1),  V ~ N(0,1), independent
#   A = a0 + a1*U + a2*V + eA
#   Y = b0 + b1*U + b2*A + eY
#   N = c0 + c1*U + c2*V + eN
#
# Var(eA) and Var(eN) are chosen so that Var(A)=Var(N)=1.
# For first DGM a2 = c2 = 0
# This is achieved by setting Var(eA) = 1 - a1^2 - a2^2 and Var(eN) = 1 - c1^2 - c2^2.

## The main function is DGM() below. The rest of the functions are internal helpers.
# Compute residual SD to enforce Var(linear combo + residual)=1
CalcResidSD <- function(a1, a2, c1, c2) {
  v_eA  <- 1 - (a1^2 + a2^2)
  v_eN <- 1 - (c1^2 + c2^2)
  if (v_eA <= 0 | v_eN <= 0 ) {
    stop(sprintf(
      "Residual variances must be positive. Computed Var(eA)=%.6f and Var(eN)=%.6f.",
      v_eA, v_eN
    ))
  }
  c(sqrt(v_eA), sqrt(v_eN))
}

# Sample mean-zero, and given variance level noise with optional heavy tails
rnoise <- function(n_sample, dist = c("norm", "t"), df = 5, sigma) {
  dist <- match.arg(dist)
  if (dist == "norm") {
    noise <- rnorm(n_sample, mean = 0, sd = sigma)
  } else {
    if (dist == "t" && df <= 2) stop("For t noise, df must be > 2.")
    noise_temp <- rt(n_sample, df = df)
    noise <- sigma * noise_temp / sqrt(df / (df - 2))
  }
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
  sigma_eN <- sd_resids[2]
  eA  <- rnoise(n_sample, dist = noise_dist, df = df, sigma = sigma_eA)  
  eN <- rnoise(n_sample, dist = noise_dist, df = df, sigma = sigma_eN)
  eY  <- rnoise(n_sample, dist = noise_dist, df = df, sigma = sigma_eY) 
  U <- rnorm(n_sample)
  V <- rnorm(n_sample)
  A <- a0 + a1 * U + a2 * V + eA
  Y <- b0 + b1 * U + b2 * A + eY
  N <- c0 + c1 * U + c2 * V + eN
  
  out <- data.frame(Y = Y, A = A, N = N, U, V) # return U,V for checking purposes only
  out
}



# ---- quick tests ----
# pars0 <- list(a0=0, a1=0.5, a2=0, b0=0, b1=0.3, b2=0.4, c0=0, c1=0.6, c2=0, sigma_eY=1)
# dat0 <-DGM(n_sample=1e6, pars=pars0, noise_dist="norm")
# ff <- (lm(formula = Y ~ A + U + V + N, data = dat0))
# var(dat0$A)      
# var(dat0$N)

