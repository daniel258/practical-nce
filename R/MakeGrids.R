# MakeGrids.R
# Grid builders for Practical Considerations NCE simulations
#
# Conventions (consistent with DGM.R):
#   A      = a0 + a1*U + a2*V + eA
#   Atilde = c0 + c1*U + c2*V + eAt
#   Y      = b0 + b1*U + b2*A + eY
#
# Var(A)=Var(Atilde)=1 enforced in DGM.R by choosing Var(eA), Var(eAt).
# Therefore must have a1^2 + a2^2 < 1 and c1^2 + c2^2 < 1.

# ---- internal helpers ----

AssertParsFeasible <- function(a1, a2, c1, c2, tol = 1e-12) {
  vA  <- 1 - (a1^2 + a2^2)
  vAt <- 1 - (c1^2 + c2^2)
  if (any(vA <= tol) || any(vAt <= tol)) {
    stop("Infeasible parameters: need a1^2+a2^2 < 1 and c1^2+c2^2 < 1")
  }
  invisible(TRUE)
}

AddDerivedCols <- function(grid) {
  # With independent U,V and Var(A)=Var(Atilde)=1:
  # corr(A, Atilde) = Cov(A, Atilde) = a1*c1 + a2*c2
  grid$rho_U     <- with(grid, a1 * c1)
  grid$rho_V     <- with(grid, a2 * c2)
  grid$rho_total <- with(grid, rho_U + rho_V)
  # Avoid division by zero
  grid$rho_rel_U <- NA_real_
  ok <- abs(grid$rho_total) > 1e-12
  grid$rho_rel_U[ok] <- grid$rho_U[ok] / grid$rho_total[ok]
  grid$bias     <- with(grid, b1 * a1)  # bias in the Y ~ A model
  grid
}

FinalizeGrid <- function(grid,
                         a0 = 0, b0 = 0, c0 = 0,
                         sigma_eY = 1,
                         add_derived = TRUE) {
  grid$a0 <- a0; grid$b0 <- b0; grid$c0 <- c0
  grid$sigma_eY <- sigma_eY
  
  # Feasibility checks row-wise
  AssertParsFeasible(grid$a1, grid$a2, grid$c1, grid$c2)
  
  if (add_derived) grid <- AddDerivedCols(grid)
  grid
}

# ---- basic (no V) ----

MakeGrid_NoV <- function(a1, c1, b1, b2,
                         a0 = 0, b0 = 0, c0 = 0,
                         sigma_eY = 1,
                         add_derived = TRUE) {
  grid <- expand.grid(a1 = a1, c1 = c1, b1 = b1, b2 = b2, stringsAsFactors = FALSE)
  grid$a2 <- 0
  grid$c2 <- 0
  FinalizeGrid(grid, a0 = a0, b0 = b0, c0 = c0, sigma_eY = sigma_eY, add_derived = add_derived)
}

# ---- Design 1 ----
# Fix total correlation rho and fix confounding strength bias = b1*a1.
# Default "symmetric" choice: c1=a1 and c2=a2 (you can change that later if needed).
#
# With c1=a1 and c2=a2:
#   rho = a1^2 + a2^2  =>  a2 = sqrt(rho - a1^2)
# Constraints require: 0 < rho < 1 and a1^2 < rho.
MakeGrid_D1_FixRhoFixBias <- function(a1_vec,
                                       rho,
                                       bias,
                                       b2,
                                       a0 = 0, b0 = 0, c0 = 0,
                                       sigma_eY = 1,
                                       add_derived = TRUE) {
  if (rho <= 0 || rho >= 1) stop("rho must be in (0,1).")
  
  a1 <- a1_vec
  if (any(a1 <= 0)) stop("a1_vec must be > 0 (needed for b1 = bias/a1).")
  if (any(a1^2 >= rho)) stop("Need a1^2 < rho for a2 = sqrt(rho - a1^2) to be real/positive.")
  
  a2 <- sqrt(rho - a1^2)
  
  # Symmetric choice (easy to interpret; matches rho decomposition cleanly)
  c1 <- a1
  c2 <- a2
  
  # Fix confounding strength: bias = b1*a1 -> b1 = bias/a1
  b1 <- bias / a1
  
  grid <- data.frame(
    a1 = a1,
    a2 = a2,
    c1 = c1,
    c2 = c2,
    b1 = b1,
    b2 = b2,
    stringsAsFactors = FALSE
  )
  
  FinalizeGrid(grid, a0 = a0, b0 = b0, c0 = c0, sigma_eY = sigma_eY, add_derived = add_derived)
}

# ---- Design 2 ----
# Fix confounding strength bias (via b1 = bias/a1) and vary rho_total.

# D2: fix (a1,c1) and vary (a2,c2) symmetrically to change rho.
# With c1=a1 and c2=a2, rho = a1^2 + a2^2, so varying a2 varies rho.
MakeGrid_D2_FixBiasVaryRho_ByA2 <- function(a1,
                                             a2_vec,
                                             bias,
                                             b2,
                                             a0 = 0, b0 = 0, c0 = 0,
                                             sigma_eY = 1,
                                             add_derived = TRUE) {
  if (a1 <= 0) stop("a1 must be > 0.")
  if (any(a2_vec < 0)) stop("a2_vec must be >= 0.")
  if (any(a1^2 + a2_vec^2 >= 1)) stop("Need a1^2 + a2^2 < 1 for all a2.")
  
  a2 <- a2_vec
  c1 <- rep(a1, length(a2))
  c2 <- a2
  
  b1 <- rep(bias / a1, length(a2))
  
  grid <- data.frame(
    a1 = rep(a1, length(a2)),
    a2 = a2,
    c1 = c1,
    c2 = c2,
    b1 = b1,
    b2 = b2,
    stringsAsFactors = FALSE
  )
  
  FinalizeGrid(grid, a0 = a0, b0 = b0, c0 = c0, sigma_eY = sigma_eY, add_derived = add_derived)
}

# D3: fix (a2,c2) and vary (a1,c1) symmetrically to change rho.
# With c1=a1 and c2=a2 fixed: rho = a1^2 + a2^2.
MakeGrid_D3_FixBiasVaryRho_ByA1 <- function(a1_vec,
                                             a2,
                                             bias,
                                             b2,
                                             a0 = 0, b0 = 0, c0 = 0,
                                             sigma_eY = 1,
                                             add_derived = TRUE) {
  if (a2 < 0) stop("a2 must be >= 0.")
  if (any(a1_vec <= 0)) stop("a1_vec must be > 0.")
  if (any(a1_vec^2 + a2^2 >= 1)) stop("Need a1^2 + a2^2 < 1 for all a1.")
  
  a1 <- a1_vec
  c1 <- a1
  c2 <- rep(a2, length(a1))
  
  b1 <- bias / a1
  
  grid <- data.frame(
    a1 = a1,
    a2 = rep(a2, length(a1)),
    c1 = c1,
    c2 = c2,
    b1 = b1,
    b2 = b2,
    stringsAsFactors = FALSE
  )
  
  FinalizeGrid(grid, a0 = a0, b0 = b0, c0 = c0, sigma_eY = sigma_eY, add_derived = add_derived)
}
