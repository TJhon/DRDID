#' @import stats
NULL
###################################################################################
#' Outcome regression DiD estimator for the ATT, with panel data
#'
#' @description \code{reg_did_panel} computes the outcome regressions estimators for the average treatment effect on the
#' treated in difference-in-differences (DiD) setups with panel data.
#'
#' @param y1 An \eqn{n} x \eqn{1} vector of outcomes from the post-treatment period.
#' @param y0 An \eqn{n} x \eqn{1} vector of outcomes from the pre-treatment period.
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
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = TRUE, boot, boot.type, nboot, type="or")}
#'
#' @details
#'
#' The \code{reg_did_panel} function implements
#' outcome regression difference-in-differences (DiD) estimator for the average treatment effect
#' on the treated (ATT) defined in equation (2.2) of Sant'Anna and Zhao (2020) when panel data are available.
#' The estimator follows the same spirit of the nonparametric estimators proposed by Heckman, Ichimura and Todd (1997),
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
#' @examples
#' # Form the Lalonde sample with CPS comparison group
#' eval_lalonde_cps <- subset(nsw, nsw$treated == 0 | nsw$sample == 2)
#' # Further reduce sample to speed example
#' set.seed(123)
#' unit_random <- sample(1:nrow(eval_lalonde_cps), 5000)
#' eval_lalonde_cps <- eval_lalonde_cps[unit_random,]
#' # Select some covariates
#' covX = as.matrix(cbind(eval_lalonde_cps$age, eval_lalonde_cps$educ,
#'                        eval_lalonde_cps$black, eval_lalonde_cps$married,
#'                        eval_lalonde_cps$nodegree, eval_lalonde_cps$hisp,
#'                        eval_lalonde_cps$re74))
#' # Implement OR DiD with panel data
#' reg_did_panel(y1 = eval_lalonde_cps$re78, y0 = eval_lalonde_cps$re75,
#'                D = eval_lalonde_cps$experimental,
#'                covariates = covX)
#'


#' @export

reg_did_panel <-function(y1, y0, D, covariates, i_weights = NULL,
                         boot = FALSE, boot.type = "weighted", nboot = NULL,
                         inffunc = FALSE){
  #-----------------------------------------------------------------------------
  print("reg_did_panel")
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # generate deltaY
  deltaY <- as.vector(y1 - y0)
  # Add constant to covariate vector
  int_cov <- has_intercept(covariates)

  # Weights
  i_weights <- has_weights(i_weights)
  #-----------------------------------------------------------------------------
  #Compute the Outcome regression for the control group using ols.

  out_delta <- out_wols(deltaY, int_cov, D == 0, i_weights)
  #-----------------------------------------------------------------------------
  #Compute the OR-DiD estimator
  # First, the weights
  w_treat <- w_cont <- i_weights * D

  reg_att_treat <- w_treat * deltaY
  reg_att_cont  <- w_cont * out_delta

  eta_treat <- eta_val(reg_att_treat, w_treat)
  eta_cont  <- eta_val(reg_att_cont, w_cont)

  reg_att <-   eta_treat - eta_cont
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  # First, the influence function of the nuisance functions
  # Asymptotic linear representation of OLS parameters
  asy_lin_rep_ols <- asy_lin_rep_olsf(
    i_weights, (1 - D), 1, int_cov, deltaY, out_delta
  )
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function
  inf_treat <- inf_treatf(reg_att, w_treat, eta_treat)
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect
  inf_cont_1 <- (reg_att_cont - w_cont * eta_cont)
  # Estimation effect from beta hat (OLS using only controls)
  # Derivative matrix (k x 1 vector)
  M1 <- base::colMeans(w_cont * int_cov)
  # Now get the influence function related to the estimation effect related to beta's
  inf_cont_2 <- asy_lin_rep_ols %*% M1
  # Influence function for the control component
  inf_control <- (inf_cont_1 + inf_cont_2) / mean(w_cont)
  #-----------------------------------------------------------------------------
  #get the influence function of the DR estimator (put all pieces together)
  reg_att_inf_func <- (inf_treat - inf_control)
  #-----------------------------------------------------------------------------

  setup <- list(
    n = n,
    y = deltaY,
    d = D,
    int_cov = int_cov,
    i_weights = i_weights,
    reg_att = reg_att
  )

  ref_se <- bstrap_se(reg_att_inf_func, boot, nboot, boot.type, setup, wboot.reg.panel)



  if(inffunc == FALSE) reg.att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot.type <- ifelse(argu$boot.type=="multiplier", "multiplier", "weighted")
  argu <- list(
    panel = TRUE,
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
