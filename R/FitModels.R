# FitModels.R
# Fit working regressions and extract SEs + (empirical) power indicators (reject at alpha).
# - Fit 1: Y ~ A
# - Fit 2: Y ~ A + Atilde
# Robust SEs use sandwich::vcovHC (HC1).

FitModels <- function(dat, alpha = 0.05, intercept = TRUE, robust_se = FALSE) {
  
  if (intercept) {
    f1 <- Y ~ A
    f2 <- Y ~ A + Atilde
    f3 <- Y ~ A + V 
    f4 <- Y ~ A + Atilde + V
  } else {
    f1 <- Y ~ A - 1
    f2 <- Y ~ A + Atilde - 1
    f3 <- Y ~ A + V - 1
    f4 <- Y ~ A + Atilde + V - 1
  }
  
  fit1 <- lm(f1, data = dat)
  fit2 <- lm(f2, data = dat)
  fit3 <- lm(f3, data = dat)
  fit4 <- lm(f4, data = dat)
  
  # point estimates
  beta_A_fit1      <- unname(coef(fit1)["A"])
  beta_A_fit2      <- unname(coef(fit2)["A"])
  beta_Atilde_fit2 <- unname(coef(fit2)["Atilde"])
  beta_A_fit3      <- unname(coef(fit3)["A"])
  beta_V_fit3      <- unname(coef(fit3)["V"])
  beta_A_fit4      <- unname(coef(fit4)["A"])
  beta_Atilde_fit4 <- unname(coef(fit4)["Atilde"])
  beta_V_fit4      <- unname(coef(fit4)["V"])
  
  if (!robust_se) {
    s1 <- summary(fit1)$coefficients
    s2 <- summary(fit2)$coefficients
    s3 <- summary(fit3)$coefficients
    s4 <- summary(fit4)$coefficients
    
    se_A_fit1      <- unname(s1["A", "Std. Error"])
    t_A_fit1       <- unname(s1["A", "t value"])
    p_A_fit1       <- unname(s1["A", "Pr(>|t|)"])
    
    se_A_fit2      <- unname(s2["A", "Std. Error"])
    t_A_fit2       <- unname(s2["A", "t value"])
    p_A_fit2       <- unname(s2["A", "Pr(>|t|)"])
    
    se_Atilde_fit2 <- unname(s2["Atilde", "Std. Error"])
    t_Atilde_fit2  <- unname(s2["Atilde", "t value"])
    p_Atilde_fit2  <- unname(s2["Atilde", "Pr(>|t|)"])
    
    se_A_fit3 <- unname(s3["A", "Std. Error"])
    p_A_fit3  <- unname(s3["A", "Pr(>|t|)"])
    se_V_fit3 <- unname(s3["V", "Std. Error"])
    p_V_fit3  <- unname(s3["V", "Pr(>|t|)"])
    
    se_A_fit4      <- unname(s4["A", "Std. Error"])
    p_A_fit4       <- unname(s4["A", "Pr(>|t|)"])
    se_Atilde_fit4 <- unname(s4["Atilde", "Std. Error"])
    p_Atilde_fit4  <- unname(s4["Atilde", "Pr(>|t|)"])
    se_V_fit4      <- unname(s4["V", "Std. Error"])
    p_V_fit4       <- unname(s4["V", "Pr(>|t|)"])
    
  } else {
    if (!requireNamespace("sandwich", quietly = TRUE)) {
      stop("Package 'sandwich' is required for robust_se=TRUE. Please install it.", call. = FALSE)
    }
    
    Vcov1 <- sandwich::vcovHC(fit1, type = "HC1")
    Vcov2 <- sandwich::vcovHC(fit2, type = "HC1")
    Vcov3 <- sandwich::vcovHC(fit3, type = "HC1")
    Vcov4 <- sandwich::vcovHC(fit4, type = "HC1")
    
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
    
    df3 <- fit3$df.residual
    se_A_fit3 <- sqrt(Vcov3["A", "A"])
    p_A_fit3  <- 2 * pt(abs(beta_A_fit3 / se_A_fit3), df = df3, lower.tail = FALSE)
    se_V_fit3 <- sqrt(Vcov3["V", "V"])
    p_V_fit3  <- 2 * pt(abs(beta_V_fit3 / se_V_fit3), df = df3, lower.tail = FALSE)
    
    df4 <- fit4$df.residual
    se_A_fit4 <- sqrt(Vcov4["A", "A"])
    p_A_fit4  <- 2 * pt(abs(beta_A_fit4 / se_A_fit4), df = df4, lower.tail = FALSE)
    se_Atilde_fit4 <- sqrt(Vcov4["Atilde", "Atilde"])
    p_Atilde_fit4  <- 2 * pt(abs(beta_Atilde_fit4 / se_Atilde_fit4), df = df4, lower.tail = FALSE)
    se_V_fit4 <- sqrt(Vcov4["V", "V"])
    p_V_fit4  <- 2 * pt(abs(beta_V_fit4 / se_V_fit4), df = df4, lower.tail = FALSE)
    
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
    power_Atilde_fit2 = as.numeric(p_Atilde_fit2 < alpha),
    
    beta_A_fit3      = beta_A_fit3,
    se_A_fit3        = se_A_fit3,
    p_A_fit3         = p_A_fit3,
    power_A_fit3     = as.numeric(p_A_fit3 < alpha),
    
    beta_V_fit3      = beta_V_fit3,
    se_V_fit3        = se_V_fit3,
    p_V_fit3         = p_V_fit3,
    power_V_fit3     = as.numeric(p_V_fit3 < alpha),
    
    beta_A_fit4      = beta_A_fit4,
    se_A_fit4        = se_A_fit4,
    p_A_fit4         = p_A_fit4,
    power_A_fit4     = as.numeric(p_A_fit4 < alpha),
    
    beta_Atilde_fit4  = beta_Atilde_fit4,
    se_Atilde_fit4    = se_Atilde_fit4,
    p_Atilde_fit4     = p_Atilde_fit4,
    power_Atilde_fit4 = as.numeric(p_Atilde_fit4 < alpha),
    
    beta_V_fit4      = beta_V_fit4,
    se_V_fit4        = se_V_fit4,
    p_V_fit4         = p_V_fit4,
    power_V_fit4     = as.numeric(p_V_fit4 < alpha)
  )
}
