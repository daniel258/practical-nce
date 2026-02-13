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
# Fix bias and fix rho_U = a1*c1, then vary rho_total via rho_V = a2*c2 by varying a2.
# Interpretation: corr(A, Atilde) increases "via V" (not via U).

MakeGrid_D2_FixBiasVaryRho_ByA2 <- function(a1,
                                            a2_vec,
                                            c1,
                                            c2,
                                            bias,
                                            b2,
                                            a0 = 0, b0 = 0, c0 = 0,
                                            sigma_eY = 1,
                                            add_derived = TRUE) {
  if (length(a1) != 1 || !is.finite(a1) || a1 == 0) stop("a1 must be a single finite nonzero value.")
  if (any(!is.finite(a2_vec))) stop("a2_vec must be finite.")
  if (length(c1) != 1 || !is.finite(c1)) stop("c1 must be a single finite value.")
  if (length(c2) != 1 || !is.finite(c2)) stop("c2 must be a single finite value.")
  
  # Basic feasibility (FinalizeGrid will also check)
  if (any(a1^2 + a2_vec^2 >= 1)) stop("Need a1^2 + a2^2 < 1 for all a2.")
  if (c1^2 + c2^2 >= 1) stop("Need c1^2 + c2^2 < 1.")
  
  # Fix confounding strength: bias = b1*a1 -> b1 = bias/a1
  b1 <- rep(bias / a1, length(a2_vec))
  
  # Allow b2 scalar or length(a2_vec)
  if (length(b2) == 1) {
    b2_use <- rep(b2, length(a2_vec))
  } else if (length(b2) == length(a2_vec)) {
    b2_use <- b2
  } else {
    stop("b2 must be length 1 or length(a2_vec).")
  }
  
  grid <- data.frame(
    a1 = rep(a1, length(a2_vec)),
    a2 = a2_vec,
    c1 = rep(c1, length(a2_vec)),
    c2 = rep(c2, length(a2_vec)),
    b1 = b1,
    b2 = b2_use,
    stringsAsFactors = FALSE
  )
  
  FinalizeGrid(grid, a0 = a0, b0 = b0, c0 = c0, sigma_eY = sigma_eY, add_derived = add_derived)
}



# ---- Design 3 ----
# Fix bias and fix rho_V = a2*c2, then vary rho_total via rho_U = a1*c1 by varying c1.
# Interpretation: corr(A, Atilde) increases "via U" (contrast to Design 2).

MakeGrid_D3_FixBiasVaryRho_ByC1 <- function(c1_vec,
                                            a1,
                                            a2,
                                            c2,
                                            bias,
                                            b2,
                                            a0 = 0, b0 = 0, c0 = 0,
                                            sigma_eY = 1,
                                            add_derived = TRUE) {
  if (length(a1) != 1 || !is.finite(a1) || a1 == 0) stop("a1 must be a single finite nonzero value.")
  if (length(a2) != 1 || !is.finite(a2)) stop("a2 must be a single finite value.")
  if (length(c2) != 1 || !is.finite(c2)) stop("c2 must be a single finite value.")
  if (any(!is.finite(c1_vec))) stop("c1_vec must be finite.")
  
  # Basic feasibility (FinalizeGrid will also check)
  if (a1^2 + a2^2 >= 1) stop("Need a1^2 + a2^2 < 1.")
  if (any(c1_vec^2 + c2^2 >= 1)) stop("Need c1^2 + c2^2 < 1 for all c1.")
  
  # Fix confounding strength: bias = b1*a1 -> b1 = bias/a1
  b1 <- rep(bias / a1, length(c1_vec))
  
  # Allow b2 scalar or length(c1_vec)
  if (length(b2) == 1) {
    b2_use <- rep(b2, length(c1_vec))
  } else if (length(b2) == length(c1_vec)) {
    b2_use <- b2
  } else {
    stop("b2 must be length 1 or length(c1_vec).")
  }
  
  grid <- data.frame(
    a1 = rep(a1, length(c1_vec)),
    a2 = rep(a2, length(c1_vec)),
    c1 = c1_vec,
    c2 = rep(c2, length(c1_vec)),
    b1 = b1,
    b2 = b2_use,
    stringsAsFactors = FALSE
  )
  
  FinalizeGrid(grid, a0 = a0, b0 = b0, c0 = c0, sigma_eY = sigma_eY, add_derived = add_derived)
}
