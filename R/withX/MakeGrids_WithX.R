# MakeGrids_WithX.R
# Grid builders for Practical Considerations NCE simulations (With measured covariates).
#
#  corr(A, Atilde) remains:
#   rho_total = a1*c1 + a2*c2   (since Var(A)=Var(Atilde)=1 and drivers independent)
#
#  feasibility for Var(A)=1 now requires:
#   a1^2 + a2^2 + ax_1^2 + ax_2^2 < 1
# while Atilde feasibility is unchanged:
#   c1^2 + c2^2 < 1

AssertParsFeasible_WithX <- function(a1, a2, c1, c2, ax_1, ax_2, tol = 1e-12) {
  vA  <- 1 - (a1^2 + a2^2 + ax_1^2 + ax_2^2)
  vAt <- 1 - (c1^2 + c2^2)
  if (any(vA <= tol) || any(vAt <= tol)) {
    stop("Infeasible parameters: need a1^2+a2^2+ax_c^2+ax_b^2 < 1 and c1^2+c2^2 < 1")
  }
  invisible(TRUE)
}

AddDerivedCols_WithX <- function(grid, tol = 1e-12) {
  # corr(A, Atilde) is still a1*c1 + a2*c2 (X omitted from Atilde)
  grid$rho_U     <- with(grid, a1 * c1)
  grid$rho_V     <- with(grid, a2 * c2)
  grid$rho_total <- with(grid, rho_U + rho_V)

  grid$rho_rel_U <- NA_real_
  ok <- abs(grid$rho_total) > tol
  grid$rho_rel_U[ok] <- grid$rho_U[ok] / grid$rho_total[ok]

  # "Structural" confounding strength via U in Y~A (ignoring conditioning)
  grid$bias_struct <- with(grid, b1 * a1)

  # Bias (population) in Fit 1 and Fit 3 due to U, given the regressors included:
  # Fit 1 includes X, so residualized A has Var = 1 - (ax_c^2 + ax_b^2)
  grid$var_A_given_X <- with(grid, 1 - (ax_1^2 + ax_2^2))
  grid$bias_fit1_pop <- with(grid, (b1 * a1) / var_A_given_X)

  # Fit 3 includes X and V, so residualized A has Var = 1 - (ax_c^2 + ax_b^2 + a2^2)
  grid$var_A_given_XV <- with(grid, 1 - (ax_1^2 + ax_2^2 + a2^2))
  grid$bias_fit3_pop  <- with(grid, (b1 * a1) / var_A_given_XV)

  grid
}

FinalizeGrid_WithX <- function(grid,
                               a0 = 0, b0 = 0, c0 = 0,
                               sigma_eY = 1,
                               ax_1, ax_2, bx_1, bx_2, p_x2 = 0.5,
                               add_derived = TRUE) {

  grid$a0 <- a0; grid$b0 <- b0; grid$c0 <- c0
  grid$sigma_eY <- sigma_eY

  # Attach X coefficients as grid columns (constants within a given study)
  grid$ax_1 <- ax_1
  grid$ax_2 <- ax_2
  grid$bx_1 <- bx_1
  grid$bx_2 <- bx_2
  grid$p_x2 <- p_x2

  # Feasibility checks row-wise
  AssertParsFeasible_WithX(grid$a1, grid$a2, grid$c1, grid$c2, grid$ax_1, grid$ax_2)

  if (add_derived) grid <- AddDerivedCols_WithX(grid)
  grid
}

# ---- basic (no V) ----

MakeGrid_NoV_WithX <- function(a1, c1, b1, b2,
                               ax_1, ax_2, bx_1, bx_2, p_x2 = 0.5,
                               a0 = 0, b0 = 0, c0 = 0,
                               sigma_eY = 1,
                               add_derived = TRUE) {
  grid <- expand.grid(a1 = a1, c1 = c1, b1 = b1, b2 = b2, stringsAsFactors = FALSE)
  grid$a2 <- 0
  grid$c2 <- 0
  FinalizeGrid_WithX(grid,
                     a0 = a0, b0 = b0, c0 = c0,
                     sigma_eY = sigma_eY,
                     ax_1 = ax_1, ax_2 = ax_2, bx_1 = bx_1, bx_2 = bx_2, p_x2 = p_x2,
                     add_derived = add_derived)
}

# ---- Design 1 ----
# Fix U->A and U->Atilde; vary V->A and V->Atilde via a2=c2 to hit rho_target.
MakeGrid_D1_FixU_VaryV_WithX <- function(rho_target_vec,
                                         a1, c1,
                                         b0 = 0, b1 = 0, b2,
                                         ax_1, ax_2, bx_1, bx_2, p_x2 = 0.5,
                                         a0 = 0, c0 = 0,
                                         sigma_eY = 1,
                                         add_derived = TRUE,
                                         tol = 1e-12) {

  rho_U <- a1 * c1

  grid_list <- lapply(rho_target_vec, function(rho_target) {

    t2 <- rho_target - rho_U
    if (t2 < -tol) return(NULL)

    t2 <- max(t2, 0)
    t  <- sqrt(t2)

    a2 <- t
    c2 <- t

    # feasibility for A includes X
    if ((a1^2 + a2^2 + ax_1^2 + ax_2^2) >= (1 - tol)) return(NULL)
    if ((c1^2 + c2^2) >= (1 - tol)) return(NULL)

    data.frame(
      design     = "D1",
      rho_target = rho_target,
      a1 = a1, a2 = a2,
      c1 = c1, c2 = c2,
      b0 = b0, b1 = b1, b2 = b2,
      stringsAsFactors = FALSE
    )
  })

  grid <- do.call(rbind, grid_list)
  if (is.null(grid) || nrow(grid) == 0) {
    stop("D1 grid is empty: rho_target_vec may be below rho_U=a1*c1 or violates feasibility constraints.")
  }

  grid <- FinalizeGrid_WithX(grid,
                             a0 = a0, b0 = b0, c0 = c0,
                             sigma_eY = sigma_eY,
                             ax_1 = ax_1, ax_2 = ax_2, bx_1 = bx_1, bx_2 = bx_2, p_x2 = p_x2,
                             add_derived = add_derived)

  # share of rho attributable to V (targeted)
  grid$f_target <- NA_real_
  ok <- abs(grid$rho_total) > tol
  grid$f_target[ok] <- grid$rho_V[ok] / grid$rho_total[ok]

  rownames(grid) <- NULL
  grid
}

# ---- Design 2 ----
# Fix rho_total and vary f = rho_V/rho_total; fix a1 (and a2) so bias is fixed up to a constant (given X adjustment).
MakeGrid_D2_FixBiasFixRho_VaryF_WithX <- function(
    rho_total,
    f_vec,
    a1,
    a2 = a1,
    b0 = 0, b1, b2,
    ax_1, ax_2, bx_1, bx_2, p_x2 = 0.5,
    a0 = 0, c0 = 0,
    sigma_eY = 1,
    add_derived = TRUE,
    tol = 1e-12
) {
  if (length(rho_total) != 1) stop("rho_total must be a single value.")
  if (any(f_vec < 0 | f_vec > 1)) stop("Need f in [0,1].")
  if (abs(a1) <= tol) stop("Need a1 != 0.")
  if (abs(a2) <= tol) stop("Need a2 != 0.")

  # feasibility for A includes X
  if ((a1^2 + a2^2 + ax_1^2 + ax_2^2) >= (1 - tol)) {
    stop("Infeasible: need a1^2 + a2^2 + ax_c^2 + ax_b^2 < 1.")
  }

  grid <- expand.grid(
    f_target = f_vec,
    b1 = b1,
    stringsAsFactors = FALSE
  )

  grid$design     <- "D2"
  grid$rho_target <- rho_total

  grid$a1 <- a1
  grid$a2 <- a2
  grid$b0 <- b0
  grid$b2 <- b2

  # Induce fixed rho_total with varying decomposition
  grid$c1 <- ((1 - grid$f_target) * rho_total) / a1
  grid$c2 <- (grid$f_target * rho_total) / a2

  ok <- (grid$c1^2 + grid$c2^2) < (1 - tol)
  grid <- grid[ok, , drop = FALSE]
  if (nrow(grid) == 0) stop("D2 grid is empty: try increasing a1/a2 or decreasing |rho_total|.")

  grid <- FinalizeGrid_WithX(grid,
                             a0 = a0, b0 = b0, c0 = c0,
                             sigma_eY = sigma_eY,
                             ax_1 = ax_1, ax_2 = ax_2, bx_1 = bx_1, bx_2 = bx_2, p_x2 = p_x2,
                             add_derived = add_derived)

  grid$f_achieved <- NA_real_
  ok2 <- abs(grid$rho_total) > tol
  grid$f_achieved[ok2] <- grid$rho_V[ok2] / grid$rho_total[ok2]

  rownames(grid) <- NULL
  grid
}
