has_intercept <- function(covariates, n){
  int_cov <- as.matrix(rep(1, n))

  if (!is.null(covariates)){
    if(all(as.matrix(covariates)[, 1] == rep(1, n))){
      int_cov <- as.matrix(covariates)
    } else {
      int_cov <- as.matrix(cbind(1, covariates))
    }
  }
  return(int_cov)
}

has_weights <- function(i_weights, n){

  if(is.null(i_weights)) {
    i_weights <- as.vector(rep(1, n))
  } else if(min(i_weights) < 0) stop("i.weights must be non-negative")

  return(i_weights)
}


out_wols <- function(y_, x_, filter, w){
  reg_coeff <- lm(y_ ~ -1 + x_, subset = filter, weights = w) |>
    coef()

  if(anyNA(reg_coeff)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity of covariates is probably the reason for it.")
  }
  out_y <- tcrossprod(reg_coeff, x_) |> as.vector()
}

w_tc_val <- function(w, d, pst = 1){
  wtc <- w * d * pst
  return(wtc)
}

eta_val <- function(reg_att, w_tc){
  eta_r <- mean(reg_att) / mean(w_tc)
  return(eta_r)
}

asy_lin_rep_olsf <- function(w, d, pst = 1, x_, y, out_y){
  n <- length(y)
  weights_ols <- w * d * pst
  wols_x <- weights_ols * x_
  wols_ex <- weights_ols * (y - out_y) * x_
  xpx_inv <- qr.solve(crossprod(wols_x, x_) / n)
  asy_lin <- wols_ex %*% xpx_inv
  return(asy_lin)
}

inf_treatf <- function(reg_att, w_treat, eta_treat){
  inf_f_u <- reg_att - w_treat * eta_treat
  inf_f <- inf_f_u / mean(w_treat)
  return(inf_f)
}

bstrap_se <- function(inf_func, bstr, nboot = NULL, boot_type = "multiplier", setup, cb){

  n <- setup$n
  y <- setup$y
  d <- setup$d
  int_cov <- setup$int_cov
  i_weights <- setup$i_weights
  reg_att <- setup$reg_att

  if (!bstr) {
    # Estimate of standard error
    se.reg.att <- sd(inf_func)/sqrt(n)
    # Estimate of upper boudary of 95% CI
    uci <- reg_att + 1.96 * se.reg.att
    # Estimate of lower doundary of 95% CI
    lci <- reg_att - 1.96 * se.reg.att
    #Create this null vector so we can export the bootstrap draws too.
    reg.boot <- NULL
  }

  if (bstr){
    if (is.null(nboot)) nboot = 999
    if(boot_type == "multiplier"){
      # do multiplier bootstrap
      reg.boot <- mboot.did(inf_func, nboot)
      # get bootstrap std errors based on IQR
      se.reg.att <- stats::IQR(reg.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs(reg.boot/se.reg.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- reg_att + cv * se.reg.att
      # Estimate of lower doundary of 95% CI
      lci <- reg_att - cv * se.reg.att
    } else {
      # do weighted bootstrap
      boot_type <- 'weighted'
      reg.boot <- unlist(lapply(1:nboot, cb,
                                n = n, deltaY = y, D = d, int.cov = int_cov, i.weights = i_weights))
      # get bootstrap std errors based on IQR
      se.reg.att <- stats::IQR((reg.boot - reg_att)) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs((reg.boot - reg_att)/se.reg.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- reg_att + cv * se.reg.att
      # Estimate of lower doundary of 95% CI
      lci <- reg_att - cv * se.reg.att
    }
  }
  ref_se <- list(
    uci = uci,
    lci = lci,
    se = se.reg.att,
    boots = reg.boot,
    nboot = nboot
  )
  return(ref_se)
}

