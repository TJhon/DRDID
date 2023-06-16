#' @import stats
NULL
###################################################################################
# Standardized version of Abadie's IPW DID estimator

#' Standardized inverse probability weighted DiD estimator, with repeated cross-section data
#' @description \code{std_ipw_did_rc} is used to compute inverse probability weighted (IPW) estimators for the ATT
#'  in DID setups with stationary repeated cross-sectional data. IPW weights are normalized to sum up to one, that is,
#'  the estimator is of the Hajek type.
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcomes from the both pre and post-treatment periods.
#' @param post An \eqn{n} x \eqn{1} vector of Post-Treatment dummies (post = 1 if observation belongs to post-treatment period,
#'             and post = 0 if observation belongs to pre-treatment period.)
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score estimation.
#' If covariates = NULL, this leads to an unconditional DID estimator.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#' \item{ATT}{The IPW DID point estimate.}
#' \item{se}{ The IPW DID standard error}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#' \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = FALSE, normalized = TRUE, boot, boot.type, nboot, type="ipw")}

#' @references
#' \cite{Abadie, Alberto (2005), "Semiparametric Difference-in-Differences Estimators",
#' Review of Economic Studies, vol. 72(1), p. 1-19, \doi{10.1111/0034-6527.00321}.
#' }
#'
#'
#' \cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
#' "Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Vol. 219 (1), pp. 101-122,
#' \doi{10.1016/j.jeconom.2020.06.003}}
#'
#'
#' @examples
#' # use the simulated data provided in the package
#' covX = as.matrix(sim_rc[,5:8])
#' # Implement normalized IPW DID estimator
#' std_ipw_did_rc(y = sim_rc$y, post = sim_rc$post, D = sim_rc$d,
#'                covariates= covX)
#'
#' @export

std_ipw_did_rc <-function(y, post, D, covariates, i_weights = NULL,
                          boot = FALSE, boot_type = "weighted", nboot = NULL,
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
  int_cov <- has_intercept(covariates)
  # Weights
  i_weights <- has_weights(i_weights)
  #-----------------------------------------------------------------------------
  #Pscore estimation (logit) and also its fitted values
  ps <- fit_ps(D, int_cov, i_weights)
  ps_fit <- ps$fit
  #-----------------------------------------------------------------------------
  #Compute IPW estimator
  # First, the weights
  w_left0 <- i_weights * D
  w_treat_pre <- w_left0 * (1 - post)
  w_treat_post <- w_left0 * post

  w_left1 <- i_weights * ps_fit * (1 - D)
  w_cont_pre <- w_left1 * (1 - post) / (1 - ps_fit)
  w_cont_post <- w_left1 * post / (1 - ps_fit)

  # Elements of the influence function (summands)
  eta_treat_pre  <- eta_val(w_treat_pre, y = y)#w.treat.pre * y / mean(w.treat.pre)
  eta_treat_post <- eta_val(w_treat_post, y = y)#w.treat.post * y / mean(w.treat.post)
  eta_cont_pre   <- eta_val(w_cont_pre, y = y)#w.cont.pre * y / mean(w.cont.pre)
  eta_cont_post  <- eta_val(w_cont_post, y = y)#w.cont.post * y / mean(w.cont.post)

  # Estimator of each component
  att_treat_pre <- mean(eta_treat_pre)
  att_treat_post <- mean(eta_treat_post)
  att_cont_pre <- mean(eta_cont_pre)
  att_cont_post <- mean(eta_cont_post)

  # ATT estimator
  ipw_att <- (att_treat_post - att_treat_pre) - (att_cont_post - att_cont_pre)
  # print(ipw_att)
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  # Asymptotic linear representation of logit's beta's

  asy_lin_rep_ps <- asy_lin_rep_psf(i_weights, D, ps, int_cov)
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf_treat_pre  <- inf_treatf1(eta_treat_pre, w_treat_pre, att_treat_pre)# eta.treat.pre - w.treat.pre * att.treat.pre/mean(w.treat.pre)
  inf_treat_post <- inf_treatf1(eta_treat_post, w_treat_post, att_treat_post)# eta.treat.post - w.treat.post * att.treat.post/mean(w.treat.post)
  inf_treat <- inf_treat_post - inf_treat_pre
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect
  inf_cont_pre  <- inf_treatf1(eta_cont_pre, w_cont_pre, att_cont_pre) #eta.cont.pre - w.cont.pre * att.cont.pre/mean(w.cont.pre)
  inf_cont_post <- inf_treatf1(eta_cont_post, w_cont_post, att_cont_post) #eta.cont.post - w.cont.post * att.cont.post/mean(w.cont.post)
  inf_cont <- inf_cont_post - inf_cont_pre

  # Estimation effect from gamma hat (pscore)
  # Derivative matrix (k x 1 vector)
  # M2_pre <- base::colMeans(w_cont_pre *(y - att_cont_pre) * int_cov)/mean(w_cont_pre)
  # M2_post <- base::colMeans(w_cont_post *(y - att_cont_post) * int_cov)/mean(w_cont_post)
  M2_pre <-  m2_f(w_cont_pre, y, att_cont_pre, int_cov, T)
  M2_post <- m2_f(w_cont_post, y, att_cont_post, int_cov, T)

  # Now the influence function related to estimation effect of pscores
  inf_cont_ps <- asy_lin_rep_ps %*% (M2_post - M2_pre)

  # Influence function for the control component
  inf_cont <- inf_cont + inf_cont_ps

  #get the influence function of the DR estimator (put all pieces together)
  att_inf_func <- inf_treat - inf_cont
  setup <- list(
    n = n, y = y,
    d = D, int_cov = int_cov,
    i_weights = i_weights, reg_att = ipw_att
  )
  #-----------------------------------------------------------------------------
  ref_se <-
    bstrap_se(att_inf_func, boot, nboot, boot_type, setup, wboot_std_ipw_rc)

  if(inffunc == FALSE) att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot_type <- ifelse(argu$boot_type=="multiplier", "multiplier", "weighted")
  boot <- ifelse(argu$boot == TRUE, TRUE, FALSE)
  argu <- list(
    panel = FALSE,
    normalized = TRUE,
    boot = boot,
    boot.type = boot_type,
    nboot = nboot,
    type = "ipw"
  )
  ret <- (list(ATT = ipw_att,
               se =  ref_se$se,
               uci = ref_se$uci,
               lci = ref_se$lci,
               boots = ref_se$boots,
               att.inf.func = att_inf_func,
               call.param = call.param,
               argu = argu))
  # Define a new class
  class(ret) <- "drdid"

  # return the list
  return(ret)
}
