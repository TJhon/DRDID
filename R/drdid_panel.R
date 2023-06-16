#' @import stats
NULL
###################################################################################
#  Locally Efficient Doubly Robust DiD estimator with panel Data
#' Locally efficient doubly robust DiD estimator for the ATT, with panel data
#'
#' @description \code{drdid_panel} is used to compute the locally efficient doubly robust estimators for the ATT
#'  in difference-in-differences (DiD) setups with panel data.
#'
#' @param y1 An \eqn{n} x \eqn{1} vector of outcomes from the post-treatment period.
#' @param y0 An \eqn{n} x \eqn{1} vector of outcomes from the pre-treatment period.
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
#' \item{ATT}{The DR DiD point estimate.}
#' \item{se}{ The DR DiD standard error.}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT.}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT.}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL.}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL.}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = TRUE, estMethod = "trad", boot, boot.type, nboot, type="dr")}
#'
#' @references{
#'
#' \cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
#' "Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Vol. 219 (1), pp. 101-122,
#' \doi{10.1016/j.jeconom.2020.06.003}}
#' }
#' @details
#'
#' The \code{drdid_panel} function implements the locally efficient doubly robust difference-in-differences (DiD)
#' estimator for the average treatment effect on the treated (ATT) defined in equation (3.1)
#' in Sant'Anna and Zhao (2020). This estimator makes use of a logistic propensity score model for the probability
#' of being in the treated group, and of a linear regression model for the outcome evolution among the comparison units.
#'
#'
#' The propensity score parameters are estimated using maximum
#' likelihood, and the outcome regression coefficients are estimated using ordinary least squares.
#'
#'
#' @examples
#' # Form the Lalonde sample with CPS comparison group (data in wide format)
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
#'
#' # Implement traditional DR locally efficient DiD with panel data
#' drdid_panel(y1 = eval_lalonde_cps$re78, y0 = eval_lalonde_cps$re75,
#'              D = eval_lalonde_cps$experimental,
#'              covariates = covX)
#'
#' @export

drdid_panel <-function(y1, y0, D, covariates, i_weights = NULL,
                       boot = FALSE, boot_type =  "weighted", nboot = NULL,
                       inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # generate deltaY
  deltaY <- as.vector(y1 - y0)
  # Add constant to covariate vector
  int_cov <- has_intercept(covariates, n)

  # Weights
  i_weights <- has_weights(i_weights, n)
  #-----------------------------------------------------------------------------
  #Compute the Pscore by MLE
  ps <- fit_ps(D, int_cov, i_weights, T)
  print('ps')
  ps_fit <- ps[['fit']]
  print('ps_fit')

  #Compute the Outcome regression for the control group using wols
  out_delta <- out_wols(deltaY, int_cov, D == 0, i_weights)
  #-----------------------------------------------------------------------------
  #Compute Traditional Doubly Robust DiD estimators
  # First, the weights
  w_treat <- i_weights * D
  w_cont <- i_weights * ps_fit * (1 - D) / (1 - ps_fit)
  dr_att_treat <- w_treat * (deltaY - out_delta)
  dr_att_cont <- w_cont * (deltaY - out_delta)

  eta_treat <- eta_val(dr_att_treat, w_treat) #mean(dr.att.treat) / mean(w.treat)
  eta_cont  <- eta_val(dr_att_cont, w_cont) #mean(dr.att.cont) / mean(w.cont)

  dr_att <-   eta_treat - eta_cont
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  # First, the influence function of the nuisance functions
  # Asymptotic linear representation of OLS parameters
  asy_lin_rep_wols <- asy_lin_rep_olsf(i_weights, (1 - D), 1, int_cov, deltaY, out_delta)

  # Asymptotic linear representation of logit's beta's
  asy_lin_rep_ps <- asy_lin_rep_psf(i_weights, D, ps, int_cov)

  # Now, the influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf_treat_1 <- (dr_att_treat - w_treat * eta_treat)
  # Estimation effect from beta hat
  # Derivative matrix (k x 1 vector)
  M1 <- base::colMeans(w_treat * int_cov)

  # Now get the influence function related to the estimation effect related to beta's
  inf_treat_2 <- asy_lin_rep_wols %*% M1

  # Influence function for the treated component
  inf_treat <- (inf_treat_1 - inf_treat_2) / mean(w_treat)
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect
  inf_cont_1 <- (dr_att_cont - w_cont * eta_cont)
  # Estimation effect from gamma hat (pscore)
  # Derivative matrix (k x 1 vector)
  M2 <- m2_f(w_cont, deltaY, out_delta, int_cov, F, eta_cont)
  # Now the influence function related to estimation effect of pscores
  inf_cont_2 <- asy_lin_rep_ps %*% M2
  # Estimation Effect from beta hat (weighted OLS)
  M3 <-  base::colMeans(w_cont * int_cov)
  # Now the influence function related to estimation effect of regressions
  inf_cont_3 <- asy_lin_rep_wols %*% M3

  # Influence function for the control component
  inf_control <- (inf_cont_1 + inf_cont_2 - inf_cont_3) / mean(w_cont)

  #get the influence function of the DR estimator (put all pieces together)
  dr_att_inf_func <- inf_treat - inf_control

  setup <- list(
    n = n, y = deltaY, d = D,
    int_cov = int_cov, i_weights = i_weights,
    reg_att = dr_att
  )
  ref_se <- bstrap_se(
    dr_att_inf_func, boot, nboot, boot_type, setup, wboot.dr.tr.panel
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
    panel = TRUE,
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
