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
MakeGrid_D1_FixRho <- function(
    rho_total_vec,
    b1_vec,
    f_vec,
    r_c = 0.92,          # must satisfy max(rho_total_vec) < r_c < 1
    b2 = 0.3,
    a0 = 0, b0 = 0, c0 = 0,
    sigma_eY = 1,
    add_derived = TRUE
) {
  tol <- 1e-12
  if (r_c <= 0 || r_c >= 1) stop("Need 0 < r_c < 1.")
  if (any(f_vec <= 0 | f_vec >= 1)) stop("Need f in (0,1).")
  if (any(abs(rho_total_vec) >= (r_c - 1e-8))) {
    stop("Need |rho_total| < r_c (otherwise a1^2+a2^2=(rho/r_c)^2 >= 1).")
  }
  
  grids <- list()
  idx <- 1
  
  for (rho_total in rho_total_vec) {
    for (b1 in b1_vec) {
      for (f in f_vec) {
        
        # Construct c on a fixed-radius circle: c1^2+c2^2 = r_c^2
        c1 <- r_c * sqrt(1 - f)
        c2 <- r_c * sqrt(f)
        
        # Construct a so that corr(A, Atilde)=rho_total and rho_V/rho = f
        a1 <- (rho_total * sqrt(1 - f)) / r_c
        a2 <- (rho_total * sqrt(f)) / r_c
        
        # Safety (should hold by construction, but keep tol margin consistent with AssertParsFeasible)
        if ((c1^2 + c2^2) >= (1 - tol)) next
        if ((a1^2 + a2^2) >= (1 - tol)) next
        
        grid <- data.frame(
          a1 = a1, a2 = a2,
          c1 = c1, c2 = c2,
          b1 = b1, b2 = b2,
          rho_target = rho_total,
          f_target = f,
          r_c = r_c,
          stringsAsFactors = FALSE
        )
        
        # Add required columns + derived columns using your existing helper
        grid <- FinalizeGrid(
          grid,
          a0 = a0, b0 = b0, c0 = c0,
          sigma_eY = sigma_eY,
          add_derived = add_derived
        )
        
        grids[[idx]] <- grid
        idx <- idx + 1
      }
    }
  }
  
  if (length(grids) == 0) return(data.frame())
  do.call(rbind, grids)
}


# ---- Design 2 ----
# Fix bias and fix  Fix U->A (a1) and U->Atilde (c1). 
# rho_U = a1*c1, then vary rho_total via rho_V = a2*c2 by varying  V->A and V->Atilde with constraint a2=c2, chosen to hit rho(A,Atilde)=rho_target.
# Interpretation: corr(A, Atilde) increases "via V" (not via U).

MakeGrid_D2_FixU_VaryV <- function(rho_target_vec,
                                   a1, c1,
                                   a0 = 0, c0 = 0,
                                   b0 = 0, b1 = 0, b2,
                                   sigma_eY = 1) {

  
  rho_U <- a1 * c1
  
  grid_list <- lapply(rho_target_vec, function(rho_target) {
    
    # Need rho_target >= rho_U so that t^2 = rho_target - rho_U is nonnegative
    t2 <- rho_target - rho_U
    if (t2 < 0) return(NULL)
    
    t <- sqrt(t2)
    a2 <- t
    c2 <- t
    
    # Feasibility: a1^2 + a2^2 < 1 and c1^2 + c2^2 < 1
    if (a1^2 + a2^2 >= 1) return(NULL)
    if (c1^2 + c2^2 >= 1) return(NULL)
    
    data.frame(
      design = "D2",
      rho_target = rho_target,
      rho_U = rho_U,
      rho_V = a2 * c2,
      a0 = a0, a1 = a1, a2 = a2,
      c0 = c0, c1 = c1, c2 = c2,
      b0 = b0, b1 = b1, b2 = b2,
      sigma_eY = sigma_eY,
      stringsAsFactors = FALSE
    )
  })
  
  grid <- do.call(rbind, grid_list)
  if (is.null(grid) || nrow(grid) == 0) {
    stop("D2 grid is empty: rho_target_vec may be below rho_U=a1*c1 or violates feasibility constraints.")
  }
  rownames(grid) <- NULL
  grid
}