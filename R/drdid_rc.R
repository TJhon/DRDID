#' @import stats
NULL
###################################################################################
#  Locally Efficient Doubly Robust DiD estimator with Repeated Cross Section Data
#' Locally efficient doubly robust DiD estimator for the ATT, with repeated cross-section data
#' @description \code{drdid_rc} is used to compute the locally efficient doubly robust estimators for the ATT
#'  in difference-in-differences (DiD) setups with stationary repeated cross-sectional data.
#
#' @param y An \eqn{n} x \eqn{1} vector of outcomes from the both pre and post-treatment periods.
#' @param post An \eqn{n} x \eqn{1} vector of Post-Treatment dummies (post = 1 if observation belongs to post-treatment period,
#'             and post = 0 if observation belongs to pre-treatment period.)
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score and regression estimation.
#' If covariates = NULL, this leads to an unconditional DiD estimator.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#' \item{ATT}{The DR DiD point estimate}
#' \item{se}{ The DR DiD standard error}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = TRUE, estMethod = "trad", boot, boot.type, nboot, type="dr")}
#'
#' @references{
#'
#' \cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
#' "Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Vol. 219 (1), pp. 101-122,
#' \doi{10.1016/j.jeconom.2020.06.003}}
#' }
#'
#'
#' @details
#'
#' The \code{drdid_rc} function implements the locally efficient doubly robust difference-in-differences (DiD)
#' estimator for the average treatment effect on the treated (ATT) defined in equation (3.4)
#' in Sant'Anna and Zhao (2020). This estimator makes use of a logistic propensity score model for the probability
#' of being in the treated group, and of (separate) linear regression models for the outcome of both treated and comparison units,
#' in both pre and post-treatment periods.
#'
#'
#' The propensity score parameters are estimated using maximum
#' likelihood, and the outcome regression coefficients are estimated using ordinary least squares;
#' see Sant'Anna and Zhao (2020) for details.
#'
#' @examples
#' # use the simulated data provided in the package
#' covX = as.matrix(sim_rc[,5:8])
#' # Implement the 'traditional' locally efficient DR DiD estimator
#' drdid_rc(y = sim_rc$y, post = sim_rc$post, D = sim_rc$d,
#'          covariates= covX)
#'
#' @export

drdid_rc <-function(y, post, D, covariates, i_weights = NULL,
                    boot = FALSE, boot_type =  "weighted", nboot = NULL,
                    inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int_cov <- has_intercept(covariates, n)

  # Weights
  i_weights <- has_weights(i_weights)
  #-----------------------------------------------------------------------------
  #Compute the Pscore by MLE
  pscore_tr <- fit_ps(D, int_cov, i_weights, T)
  ps_fit <- pscore_tr$fit

  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  pre_filter <- ((D == 0) & (post == 0))
  out_y_cont_pre <- out_wols(y, int_cov, pre_filter, i_weights)
  #Compute the Outcome regression for the control group at the post-treatment period, using ols.
  post_filter <- ((D == 0) & (post == 1))
  out_y_cont_post <- out_wols(y, int_cov, post_filter, i_weights)
  # Combine the ORs for control group
  out_y_cont <- post * out_y_cont_post + (1 - post) * out_y_cont_pre

  #Compute the Outcome regression for the treated group at the pre-treatment period, using ols.
  treat_pre <- ((D == 1) & (post == 0))
  out_y_treat_pre <- out_wols(y, int_cov, treat_pre, i_weights)
  #Compute the Outcome regression for the treated group at the post-treatment period, using ols.
  treat_post <- ((D == 1) & (post == 1))
  out_y_treat_post <- out_wols(y, int_cov, treat_post, i_weights)

  #-----------------------------------------------------------------------------
  # First, the weights
  w_d <- i_weights * D
  w_treat_pre <- w_d * (1 - post)
  w_treat_post <- w_d * post
  i_p_d <- i_weights * ps_fit * (1 - D)
  w_cont_pre  <- i_p_d * (1 - post) / (1 - ps_fit)
  w_cont_post <- i_p_d * post / (1 - ps_fit)

  w_dt1 <- w_treat_post
  w_dt0 <- w_treat_pre

  # Elements of the influence function (summands)
  y1 <- y - out_y_cont

  eta_treat_pre  <- eta_val(w_treat_pre, 1, y1)#w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta_treat_post <- eta_val(w_treat_post, 1, y1)#w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta_cont_pre   <- eta_val(w_cont_pre, 1, y1)#w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta_cont_post  <- eta_val(w_cont_post, 1, y1)#w.cont.post * (y - out.y.cont) / mean(w.cont.post)

  # extra elements for the locally efficient DRDID
  y2 <- out_y_treat_post - out_y_cont_post
  y3 <- out_y_treat_pre - out_y_cont_pre

  eta_d_post   <- eta_val(w_d, 1, y2) #w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta_dt1_post <- eta_val(w_dt1, 1, y2) #w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta_d_pre    <- eta_val(w_d, 1, y3) #w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta_dt0_pre  <- eta_val(w_dt1, 1, y3) #w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)


  # Estimator of each component
  att_treat_pre <- mean(eta_treat_pre)
  att_treat_post <- mean(eta_treat_post)
  att_cont_pre <- mean(eta_cont_pre)
  att_cont_post <- mean(eta_cont_post)


  att_d_post <- mean(eta_d_post)
  att_dt1_post <- mean(eta_dt1_post)
  att_d_pre <- mean(eta_d_pre)
  att_dt0_pre <- mean(eta_dt0_pre)


  # ATT estimator
  dr_att <- (att_treat_post - att_treat_pre) - (att_cont_post - att_cont_pre) +
    (att_d_post - att_dt1_post) - (att_d_pre - att_dt0_pre)
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  # First, the influence function of the nuisance functions

  # Asymptotic linear representation of OLS parameters in pre-period, control group
  d <- (1 - D)
  p1 <- (1 - post)

  asy_lin_rep_ols_pre <- asy_lin_rep_olsf(
    i_weights, d, p1, int_cov, y, out_y_cont_pre
  )

  # Asymptotic linear representation of OLS parameters in post-period, control group
  asy_lin_rep_ols_post <- asy_lin_rep_olsf(
    i_weights, d, post, int_cov, y, out_y_cont_post
  )

  # Asymptotic linear representation of OLS parameters in pre-period, treated
  asy_lin_rep_ols_pre_treat <- asy_lin_rep_olsf(
    i_weights, D, p1, int_cov, y, out_y_treat_pre
  )

  # Asymptotic linear representation of OLS parameters in post-period, treated
  asy_lin_rep_ols_post_treat <- asy_lin_rep_olsf(
    i_weights, D, post, int_cov, y, out_y_treat_post
  )
  # Asymptotic linear representation of logit's beta's
  asy_lin_rep_ps <- asy_lin_rep_psf(i_weights, D, pscore_tr, int_cov)
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf_treat_pre <- inf_treatf1(eta_treat_pre, w_treat_pre, att_treat_pre)
  inf_treat_post <- inf_treatf1(eta_treat_post, w_treat_post, att_treat_post)
  # Estimation effect from beta hat from post and pre-periods
  # Derivative matrix (k x 1 vector)
  M1_post <- - base::colMeans(w_treat_post * post * int_cov) / mean(w_treat_post)
  M1_pre  <- - base::colMeans(w_treat_pre * p1 * int_cov) / mean(w_treat_pre)

  # Now get the influence function related to the estimation effect related to beta's
  inf_treat_or_post <- asy_lin_rep_ols_post %*% M1_post
  inf_treat_or_pre <- asy_lin_rep_ols_pre %*% M1_pre
  inf_treat_or <- inf_treat_or_post + inf_treat_or_pre

  # Influence function for the treated component
  inf_treat <- inf_treat_post - inf_treat_pre + inf_treat_or
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect from nuisance parameters
  inf_cont_pre  <- inf_treatf1(eta_cont_pre, w_cont_pre, att_cont_pre)
  inf_cont_post <- inf_treatf1(eta_cont_post, w_cont_post, att_cont_post)

  # Estimation effect from gamma hat (pscore)
  # Derivative matrix (k x 1 vector)
  # M2.pre <- base::colMeans(w.cont.pre *(y - out.y.cont - att.cont.pre) * int.cov)/mean(w.cont.pre)
  # M2.post <- base::colMeans(w.cont.post *(y - out.y.cont - att.cont.post) * int.cov)/mean(w.cont.post)
  M2_pre <- m2_f(w_cont_pre, y, out_y_cont, int_cov, T, att_cont_pre)
  M2_post <- m2_f(w_cont_post, y, out_y_cont, int_cov, T, att_cont_post)
  # Now the influence function related to estimation effect of pscores
  inf_cont_ps <- asy_lin_rep_ps %*% (M2_post - M2_pre)

  # Estimation effect from beta hat from post and pre-periods
  # Derivative matrix (k x 1 vector)
  M3_post <- - base::colMeans(w_cont_post * post * int_cov) / mean(w_cont_post)
  M3_pre <- - base::colMeans(w_cont_pre * p1 * int_cov) / mean(w_cont_pre)

  # Now get the influence function related to the estimation effect related to beta's
  inf_cont_or_post <- asy_lin_rep_ols_post %*% M3_post
  inf_cont_or_pre  <- asy_lin_rep_ols_pre %*% M3_pre
  inf_cont_or <- inf_cont_or_post + inf_cont_or_pre

  # Influence function for the control component
  inf_cont <- inf_cont_post - inf_cont_pre + inf_cont_ps + inf_cont_or
  #-----------------------------------------------------------------------------
  #get the influence function of the inefficient DR estimator (put all pieces together)
  dr_att_inf_func1 <- inf_treat - inf_cont
  #-----------------------------------------------------------------------------
  # Now, we only need to get the influence function of the adjustment terms
  # First, the terms as if all OR parameters were known
  inf_eff1 <- inf_treatf1(eta_d_post, w_d, att_d_post)#eta_d_post - w_d * att_d_post / mean(w.d)
  inf_eff2 <- inf_treatf1(eta_dt1_post, w_dt1, att_dt1_post)#eta_dt1_post - w_dt1 * att.dt1.post/mean(w.dt1)
  inf_eff3 <- inf_treatf1(eta_d_pre, w_d, att_d_pre)#eta_d_pre - w_d * att.d.pre/mean(w.d)
  inf_eff4 <- inf_treatf1(eta_dt0_pre, w_dt0, att_dt0_pre)#eta_dt0_pre - w_dt0 * att.dt0.pre/mean(w.dt0)
  inf_eff <- (inf_eff1 - inf_eff2) - (inf_eff3 - inf_eff4)

  # Now the estimation effect of the OR coefficients
  mom_post<- base::colMeans((w_d / mean(w_d) -  w_dt1 / mean(w_dt1)) * int_cov)
  mom_pre <- base::colMeans((w_d / mean(w_d) -  w_dt0 / mean(w_dt0)) * int_cov)
  inf_or_post <- (asy_lin_rep_ols_post_treat - asy_lin_rep_ols_post) %*% mom_post
  inf_or_pre <-  (asy_lin_rep_ols_pre_treat - asy_lin_rep_ols_pre) %*% mom_pre
  inf_or <- inf_or_post - inf_or_pre
  #-----------------------------------------------------------------------------
  #get the influence function of the locally efficient DR estimator (put all pieces together)
  dr_att_inf_func <- dr_att_inf_func1 + inf_eff + inf_or

  setup <- list(
    n = n, y = y, d = D, int_cov = int_cov,
    i_weights = i_weights, reg_att = dr_att
  )

  #-----------------------------------------------------------------------------
  ref_se <- bstrap_se(
    dr_att_inf_func, boot, nboot, boot_type, setup, wboot_drdid_rc
  )

  if(inffunc == FALSE) dr.att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot_type <- ifelse(argu$boot_type=="multiplier", "multiplier", "weighted")
  boot <- ifelse(argu$boot == TRUE, TRUE, FALSE)
  argu <- list(
    panel = FALSE,
    estMethod = "trad",
    boot = boot,
    boot.type = boot_type,
    nboot = nboot,
    type = "dr"
  )

  ret <- (list(ATT = dr_att,
               se =  ref_se$se,
               uci = ref_se$uci,
               lci = ref_se$lci,
               boots = ref_se$boots,
               att.inf.func = dr_att_inf_func,
               call.param = call.param,
               argu = argu))
  # Define a new class
  class(ret) <- "drdid"

  # return the list
  return(ret)
}
