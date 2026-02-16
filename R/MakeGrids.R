# MakeGrids.R
# Grid builders for Practical Considerations NCE simulations
#
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
# Fix bias and fix U->A (a1) and U->Atilde (c1). 
# rho_U = a1*c1, then vary rho_total via rho_V = a2*c2 by varying  V->A and V->Atilde with constraint a2=c2, chosen to hit rho(A,Atilde)=rho_target.
# Interpretation: corr(A, Atilde) increases via V.

MakeGrid_D1_FixU_VaryV <- function(rho_target_vec,
                                   a1, c1,
                                   a0 = 0, c0 = 0,
                                   b0 = 0, b1 = 0, b2,
                                   sigma_eY = 1,
                                   add_derived = TRUE,
                                   tol = 1e-12) {
  
  rho_U <- a1 * c1
  
  grid_list <- lapply(rho_target_vec, function(rho_target) {
    
    # Need rho_target >= rho_U so that t^2 = rho_target - rho_U is nonnegative
    t2 <- rho_target - rho_U
    if (t2 < -tol) return(NULL)
    
    t2 <- max(t2, 0)  # guard tiny negative from floating point
    t  <- sqrt(t2)
    
    a2 <- t
    c2 <- t
    
    # Basic feasibility
    if ((a1^2 + a2^2) >= (1 - tol)) return(NULL)
    if ((c1^2 + c2^2) >= (1 - tol)) return(NULL)
    
    data.frame(
      design     = "D1",
      rho_target = rho_target,
      a1 = a1, a2 = a2,
      c1 = c1, c2 = c2,
      b1 = b1, b2 = b2,
      stringsAsFactors = FALSE
    )
  })
  
  grid <- do.call(rbind, grid_list)
  
  if (is.null(grid) || nrow(grid) == 0) {
    stop("D1 grid is empty: rho_target_vec may be below rho_U=a1*c1 or violates feasibility constraints.")
  }
  
  # Add intercepts, sigma, feasibility checks, and derived columns (rho_U/rho_V/rho_total/rho_rel_U/bias)
  grid <- FinalizeGrid(
    grid,
    a0 = a0, b0 = b0, c0 = c0,
    sigma_eY = sigma_eY,
    add_derived = add_derived
  )
  
  # Make D1 look like D2 object by adding the V-share column
  grid$f_target <- NA_real_
  ok <- abs(grid$rho_total) > tol
  grid$f_target[ok] <- grid$rho_V[ok] / grid$rho_total[ok]
  
  rownames(grid) <- NULL
  grid
}


# ---- Design 2  ----
#   - Fix bias in the naive exposure model Y ~ A by holding a1 fixed (bias = b1 * a1).
#   - Fix total corr(A, Atilde) = rho_total.
#   - Vary f = rho_V / rho_total, i.e., the share of corr(A, Atilde) attributable to V.
#
# Construction:
#   rho_total is fixed.
#   rho_U = (1 - f) * rho_total   and   rho_V = f * rho_total
#   Choose c1 and c2 so that:
#       a1*c1 = rho_U   and   a2*c2 = rho_V
#   i.e. c1 = rho_U / a1 and c2 = rho_V / a2
#
# Requirements:
#   - a1^2 + a2^2 < 1
#   - c1^2 + c2^2 < 1 for all f in [0,1]

MakeGrid_D2_FixBiasFixRho_VaryF <- function(
    rho_total,
    f_vec,
    a1,
    a2 = a1,
    b1,
    b2,
    a0 = 0, b0 = 0, c0 = 0,
    sigma_eY = 1,
    add_derived = TRUE,
    tol = 1e-12
) {
  if (length(rho_total) != 1) stop("rho_total must be a single value.")
  if (any(f_vec < 0 | f_vec > 1)) stop("Need f in [0,1].")
  if (abs(a1) <= tol) stop("Need a1 != 0.")
  if (abs(a2) <= tol) stop("Need a2 != 0.")
  if ((a1^2 + a2^2) >= (1 - tol)) {
    stop("Infeasible: need a1^2 + a2^2 < 1.")
  }
  
  grid <- expand.grid(
    f_target = f_vec,
    b1 = b1,
    stringsAsFactors = FALSE
  )
  
  grid$design     <- "D2"
  grid$rho_target <- rho_total
  
  # Fix a's and b2
  grid$a1 <- a1
  grid$a2 <- a2
  grid$b2 <- b2
  
  # Induce fixed rho_total with varying decomposition
  grid$c1 <- ((1 - grid$f_target) * rho_total) / a1
  grid$c2 <- (grid$f_target * rho_total) / a2
  
  # Drop infeasible rows (FinalizeGrid will also check)
  ok <- (grid$c1^2 + grid$c2^2) < (1 - tol)
  grid <- grid[ok, , drop = FALSE]
  
  if (nrow(grid) == 0) {
    stop("D2 grid is empty: try increasing a1/a2 or decreasing |rho_total|.")
  }
  
  grid <- FinalizeGrid(
    grid,
    a0 = a0, b0 = b0, c0 = c0,
    sigma_eY = sigma_eY,
    add_derived = add_derived
  )
  
  # Convenience: achieved share from derived columns
  grid$f_achieved <- NA_real_
  ok2 <- abs(grid$rho_total) > tol
  grid$f_achieved[ok2] <- grid$rho_V[ok2] / grid$rho_total[ok2]
  
  rownames(grid) <- NULL
  grid
}