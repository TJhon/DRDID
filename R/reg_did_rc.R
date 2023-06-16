#' @import stats
NULL
###################################################################################
#' Outcome regression DiD estimator for the ATT, with repeated cross-section data
#' @description \code{reg_did_rc} computes the outcome regressions estimators for the average treatment effect on the
#' treated in difference-in-differences (DiD) setups with stationary repeated cross-sectional data.

#' @param y An \eqn{n} x \eqn{1} vector of outcomes from the both pre and post-treatment periods.
#' @param post An \eqn{n} x \eqn{1} vector of Post-Treatment dummies (post = 1 if observation belongs to post-treatment period,
#'             and post = 0 if observation belongs to pre-treatment period.)
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the regression estimation.
#' If covariates = NULL, this leads to an unconditional DiD estimator.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#'  \item{ATT}{The OR DiD point estimate}
#'  \item{se}{The OR DiD standard error}
#'  \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#'  \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#'  \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = FALSE, boot, boot.type, nboot, type="or")}
#'
#' @details
#'
#' The \code{reg_did_rc} function implements
#' outcome regression difference-in-differences (DiD) estimator for the average treatment effect
#' on the treated (ATT) defined in equation (2.2) of Sant'Anna and Zhao (2020) when stationary repeated cross-sectional
#' data are available.  The estimator follows the same spirit of the nonparametric estimators proposed by Heckman, Ichimura and Todd (1997),
#' though here the the outcome regression models are assumed to be linear in covariates (parametric),
#'
#' The nuisance parameters (outcome regression coefficients) are estimated via ordinary least squares.

#' @references
#' \cite{Heckman, James J., Ichimura, Hidehiko, and Todd, Petra E. (1997),"Matching as an Econometric Evaluation Estimator: Evidence from Evaluating a Job Training Programme",
#' Review of Economic Studies, vol. 64(4), p. 605–654, \doi{10.2307/2971733}.
#' }
#'
#'
#' \cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
#' "Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Vol. 219 (1), pp. 101-122,
#' \doi{10.1016/j.jeconom.2020.06.003}}
#'
#'
#'
#' @examples
#' # use the simulated data provided in the package
#' covX = as.matrix(sim_rc[,5:8])
#' # Implement OR DiD estimator
#' reg_did_rc(y = sim_rc$y, post = sim_rc$post, D = sim_rc$d,
#'            covariates= covX)
#'
#' @export

reg_did_rc <-function(y, post, D, covariates, i_weights = NULL,
                      boot = FALSE, boot.type = "weighted", nboot = NULL,
                      inffunc = FALSE){
  #-----------------------------------------------------------------------------
  print('regdid_rc')
  # D as vector
  D <- as.vector(D)
  # post as vector
  post <- as.vector(post)
  # Sample size
  n <- length(D)
  # outcome of interested
  y <- as.vector(y)
  # Add constant to covariate vector

  int_cov <- has_intercept(covariates, n)

  # Weights

  i_weights <- has_weights(i_weights, n)

  #-----------------------------------------------------------------------------
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.

  pre_filter <- ((D == 0) & (post == 0))
  out_y_pre <- out_wols(y, int_cov, pre_filter, i_weights)
  #-----------------------------------------------------------------------------
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  post_filter <- ((D == 0) & (post == 1))
  out_y_post <- out_wols(y, int_cov, post_filter, i_weights)

  #-----------------------------------------------------------------------------
  #Compute the OR DiD estimators
  # First, the weights
  w_treat_pre <- w_tc_val(i_weights, D, 1 - post)
  w_treat_post <- w_tc_val(i_weights, D, post)
  w_cont <- i_weights * D

  reg_att_treat_pre <- w_treat_pre * y
  reg_att_treat_post <- w_treat_post * y
  reg_att_cont <- w_cont * (out_y_post - out_y_pre)

  eta_treat_pre  <- eta_val(reg_att_treat_pre, w_treat_pre)
  eta_treat_post <- eta_val(reg_att_treat_post, w_treat_post)
  eta_cont <- eta_val(reg_att_cont, w_cont)

  reg_att <- (eta_treat_post - eta_treat_pre) - eta_cont

  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  # First, the influence function of the nuisance functions
  # Asymptotic linear representation of OLS parameters in pre-period

  d <- (1 - D)
  asy_lin_rep_ols_pre <- asy_lin_rep_olsf(
    i_weights, d, (1 - post), int_cov, y, out_y_pre
  )

  # Asymptotic linear representation of OLS parameters in post-period
  asy_lin_rep_ols_post <- asy_lin_rep_olsf(
    i_weights, d, post, int_cov, y, out_y_post
  )
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function
  inf_treat_pre  <- inf_treatf(reg_att_treat_pre, w_treat_pre, eta_treat_pre)
  inf_treat_post <- inf_treatf(reg_att_treat_post, w_treat_post, eta_treat_post)
  inf_treat <- inf_treat_post - inf_treat_pre
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect
  inf_cont_1 <- (reg_att_cont - w_cont * eta_cont)
  # Estimation effect from beta hat (OLS using only controls)
  # Derivative matrix (k x 1 vector)
  M1 <- base::colMeans(w_cont * int_cov)
  # Now get the influence function related to the estimation effect related to beta's in post-treatment
  inf_cont_2_post <- asy_lin_rep_ols_post %*% M1
  # Now get the influence function related to the estimation effect related to beta's in pre-treatment
  inf_cont_2_pre <- asy_lin_rep_ols_pre %*% M1
  # Influence function for the control component
  inf_control <- (inf_cont_1 + inf_cont_2_post - inf_cont_2_pre) / mean(w_cont)
  #-----------------------------------------------------------------------------
  #get the influence function of the DR estimator (put all pieces together)
  reg_att_inf_func <- (inf_treat - inf_control)
  #-----------------------------------------------------------------------------
  setup <- list(
    n = n, y = y,
    d = D, int_cov = int_cov,
    i_weights = i_weights, reg_att = reg_att
  )

  ref_se <- bstrap_se(reg_att_inf_func, boot, nboot, boot.type, setup, wboot_reg_rc)



  if(inffunc == FALSE) reg_att_inf_func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot.type <- ifelse(argu$boot.type=="multiplier", "multiplier", "weighted")
  argu <- list(
    panel = FALSE,
    boot = boot,
    boot.type = boot.type,
    nboot = nboot,
    type = "or"
  )
  ret <- (list(ATT = reg_att,
               se = ref_se$se,
               uci = ref_se$uci,
               lci = ref_se$lci,
               boots = ref_se$boots,
               att.inf.func = reg_att_inf_func,
               call.param = call.param,
               argu = argu))
  # Define a new class
  class(ret) <- "drdid"
  # return the list
  return(ret)
}
