# FitModels.R
# Fit working regressions and extract SEs + (empirical) power indicators (reject at alpha).
# - Fit 1: Y ~ A
# - Fit 2: Y ~ A + Atilde
# Robust SEs use sandwich::vcovHC (HC1).

FitModels <- function(dat, alpha = 0.05, intercept = TRUE, robust_se = FALSE) {
  
  if (intercept) {
    f1 <- Y ~ A
    f2 <- Y ~ A + Atilde
  } else {
    f1 <- Y ~ A - 1
    f2 <- Y ~ A + Atilde - 1
  }
  
  fit1 <- lm(f1, data = dat)
  fit2 <- lm(f2, data = dat)
  
  # point estimates
  beta_A_fit1      <- unname(coef(fit1)["A"])
  beta_A_fit2      <- unname(coef(fit2)["A"])
  beta_Atilde_fit2 <- unname(coef(fit2)["Atilde"])
  
  if (!robust_se) {
    s1 <- summary(fit1)$coefficients
    s2 <- summary(fit2)$coefficients
    
    se_A_fit1      <- unname(s1["A", "Std. Error"])
    t_A_fit1       <- unname(s1["A", "t value"])
    p_A_fit1       <- unname(s1["A", "Pr(>|t|)"])
    
    se_A_fit2      <- unname(s2["A", "Std. Error"])
    t_A_fit2       <- unname(s2["A", "t value"])
    p_A_fit2       <- unname(s2["A", "Pr(>|t|)"])
    
    se_Atilde_fit2 <- unname(s2["Atilde", "Std. Error"])
    t_Atilde_fit2  <- unname(s2["Atilde", "t value"])
    p_Atilde_fit2  <- unname(s2["Atilde", "Pr(>|t|)"])
    
  } else {
    if (!requireNamespace("sandwich", quietly = TRUE)) {
      stop("Package 'sandwich' is required for robust_se=TRUE. Please install it.", call. = FALSE)
    }
    
    Vcov1 <- sandwich::vcovHC(fit1, type = "HC1")
    Vcov2 <- sandwich::vcovHC(fit2, type = "HC1")
    
    se_A_fit1 <- sqrt(Vcov1["A", "A"])
    t_A_fit1  <- beta_A_fit1 / se_A_fit1
    df1 <- fit1$df.residual
    p_A_fit1  <- 2 * pt(abs(t_A_fit1), df = df1, lower.tail = FALSE)
    
    se_A_fit2 <- sqrt(Vcov2["A", "A"])
    t_A_fit2  <- beta_A_fit2 / se_A_fit2
    df2 <- fit2$df.residual
    p_A_fit2  <- 2 * pt(abs(t_A_fit2), df = df2, lower.tail = FALSE)
    
    se_Atilde_fit2 <- sqrt(Vcov2["Atilde", "Atilde"])
    t_Atilde_fit2  <- beta_Atilde_fit2 / se_Atilde_fit2
    p_Atilde_fit2  <- 2 * pt(abs(t_Atilde_fit2), df = df2, lower.tail = FALSE)
  }
  
  data.frame(
    beta_A_fit1       = beta_A_fit1,
    se_A_fit1         = se_A_fit1,
    p_A_fit1          = p_A_fit1,
    power_A_fit1      = as.numeric(p_A_fit1 < alpha),
    
    beta_A_fit2       = beta_A_fit2,
    se_A_fit2         = se_A_fit2,
    p_A_fit2          = p_A_fit2,
    power_A_fit2      = as.numeric(p_A_fit2 < alpha),
    
    beta_Atilde_fit2  = beta_Atilde_fit2,
    se_Atilde_fit2    = se_Atilde_fit2,
    p_Atilde_fit2     = p_Atilde_fit2,
    power_Atilde_fit2 = as.numeric(p_Atilde_fit2 < alpha)
  )
}
