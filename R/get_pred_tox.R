#function get_MAPtox
#'Compute predictive DLT rates for a new trial from historical data
#'@rdname get_pred_tox
#'@description The function \code{get_pred_tox()} can be used to derive predictive
#'distributions of the DLT rates of a set of dose levels considered in a new
#'dose finding trial based on historical data from one or more dose-finding trials
#'based on a joint BLRM. The function supports both monotherapy and two-drug
#'combination therapy and is essentially a wrapper for \code{\link[decider:get_MAP]{get_MAP}()}
#'but additionally computes samples from the predictive DLT rates from the
#'drawn samples from the MAP prior. Note that this is equivalent to the posterior
#'distribution of the DLT rates of a trial that has not recorded
#'data yet from a joint BLRM including the specified co-data from other trials.
#'
#'The function \code{get_pred_tox()} can be used for both monotherapy and
#'combination therapy.
#'For the former, the arguments \code{dose2} and \code{dose.ref2} can
#'be left out. The function \code{get_pred_tox_mono()} uses a slightly different
#'syntax and output structure that is specifically designed for settings consisting purely of
#'monotherapy trials.
#'@param doses.of.interest Numerical that defines the set of dose levels for
#'which the meta-analytic predictive DLT rate is to be computed. Its specification
#'differs between \code{get_pred_tox()} and \code{get_pred_tox_mono()}.
#'
#'For \code{get_pred_tox()}: Can be either
#'a vector of length 2 that describes a single dose combination (where one of
#'the doses can be set to 0 for monotherapy doses), or a matrix with 2 rows in
#'which each column defines a dose combination.
#'
#'For \code{get_pred_tox_mono()}: Can be either a vector giving a set of monotherapy
#'dose levels, or a matrix with exactly one row or exactly one column. In both cases,
#'each number in the argument is interpreted as a single dose.
#'@param dose1 Numeric vector with non-negative entries. Describes the dose of
#'compound 1 administed to the patients of a cohort in the data, where cohort \code{i}
#'is assumed to have received \code{dose1[i]} of compound 1. Note that \code{dose1[i]}
#'can be set to 0 to indicate that a cohort received monotherapy with compound 2.
#'@param dose2 Numeric vector with non-negative entries. Describes the dose of
#'compound 2 administed to the patients of a cohort in the data, where cohort \code{i}
#'is assumed to have received \code{dose2[i]} of compound 2. Note that \code{dose2[i]}
#'can be set to 0 to indicate that a cohort received monotherapy with compound 1.
#' @param dose.ref1 Positive number, reference dose for compound 1.
#' @param dose.ref2 Positive number, reference dose for compound 2.
#'@param n.pat Numeric vector with non-negative whole numbers as entries. \code{n.pat[i]}
#'indicates the number of patients that were treated in cohort \code{i}.
#'@param n.dlt Numeric vector with non-negative whole numbers as entries. \code{n.dlt[i]}
#'indicates the number of patients that experienced DLT in cohort \code{i}.
#' @param trial Numerical or character vector that indicates to which trial a cohort is assigned to.
#' The entries can either be numericals or strings, where \code{1} and \code{"1"}
#' would be interpreted as the same trial. Note that the given entries are internally
#' converted to numbers from 1 to the number of different studies if not already
#' given in this form.
#' @param prior.mu Named list that specifies the distribution of the
#' hypermeans of the parameters of the combination therapy BLRM. Same format,
#' default, and requirements as documented for the function
#' \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()}.
#' @param prior.tau Named list that specifies the distribution of the
#' between-trial heterogeneities of the parameters of the combination therapy BLRM.
#' Same format, default, and requirements as documented for the function
#' \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()}.
#' @param saturating Optional logical, defaults to \code{FALSE}. If \code{TRUE}, the BLRM will be using a saturating interaction term as described in
#' \code{\link[OncoBayes2:blrm_formula_saturating]{OncoBayes2::blrm_formula_saturating}()}. Also refer to the Details section in the documentation
#' of \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()}.
#' @param dosing.intervals Optional numeric ascending positive entries between 0 and 1.
#' Defines the limits of the dosing intervals to be evaluated, similarly to
#' \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()}. The default value
#' is \code{c(0.16, 0.33)}.
#' @param return.samples Optional logical, defaults to \code{FALSE}. If \code{TRUE},
#' the function will return a matrix with samples from the MAP prior additionally
#' to a prior summary.
#' @param probs Optional numerical, defaults to \code{c(0.025, 0.5,  0.975)}.
#' Indicates which quantiles should be included in the posterior summaries.
#' @param ... Optional additional arguments that are passed to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()}. These
#' include the settings for MCMC and the random seed. Refer
#' to the \code{\link[rstan:rstan]{rstan-package}} for an overview. Note that
#' the function will automatically use \code{iter=10000} and a \code{control} argument
#' containing \code{adapt_delta=0.9} to ensure a relatively good fit if the
#' arguments \code{iter} and \code{control} for \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()}
#' are not provided manually in \code{...}.
#'
#'@details Refer to the documentation of \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()}
#'for a description of the underlying model used to compute predictive distributions
#'of the DLT rates.
#'
#'The function \code{get_pred_tox()} can optionally return the samples from the
#'predictive distribution, which can in turn be used to perform single-dose analyses
#'using a binomial-beta-mixture model. Refer to the documentation of \code{\link[decider:binomial_betamix]{binomial_betamix}()}
#'and the vignette on single dose confirmation for details.
#'@usage
#'get_pred_tox(
#'  doses.of.interest,
#'  dose1,
#'  dose2,
#'  dose.ref1,
#'  dose.ref2,
#'  n.pat,
#'  n.dlt,
#'  trial,
#'  prior.mu = list(mu_a1 =  c(logit(0.33), 2),
#'                  mu_b1 =  c(0,           1),
#'                  mu_a2 =  c(logit(0.33), 2),
#'                  mu_b2 =  c(0,           1),
#'                  mu_eta = c(0,           1.121)
#'                  ),
#'  prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
#'                   tau_b1 =  c(log(0.125), log(2)/1.96),
#'                   tau_a2 =  c(log(0.25),  log(2)/1.96),
#'                   tau_b2 =  c(log(0.125), log(2)/1.96),
#'                   tau_eta = c(log(0.125), log(2)/1.96)
#'                   ),
#'  saturating = FALSE,
#'  return.samples = FALSE,
#'  probs = c(0.025, 0.5, 0.975),
#'  dosing.intervals = c(0.16, 0.33),
#'  ...
#')
#' @seealso  \code{\link[decider:get_MAP]{get_MAP}()},\code{\link[decider:binomial_betamix]{binomial_betamix}()},
#' \code{\link[decider:fit_jointBLRM]{fit_jointBLRM}()},
#' \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()},
#' \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()}, \code{\link[rstan:rstan]{rstan-package}}.
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. \url{https://mc-stan.org}.
#'
#' Neuenschwander, B., Branson, M., & Gsponer, T. (2008). Critical aspects of the Bayesian approach to phase I cancer trials.
#' Statistics in medicine, 27(13), 2420-2439, <doi:10.1002/sim.3230>.
#'
#' Schmidli, H., Gsteiger, S., Roychoudhury, S., O'Hagan, A., Spiegelhalter, D., & Neuenschwander B. (2014).
#' Robust meta-analytic-predictive priors in clinical trials with historical control information.
#'
#' Neuenschwander, B., Matano, A., Tang, Z., Roychoudhury, S., Wandel, S., & Bailey, S. (2014). A Bayesian Industry Approach to Phase I Combination Trials in Oncology.
#' In: Zhao. W & Yang, H. (editors). Statistical methods in drug combination studies. Chapman and Hall/CRC, 95-135, <doi:10.1201/b17965>.
#'
#' Neuenschwander, B., Roychoudhury, S., & Schmidli, H. (2016). On the use of co-data in clinical trials.
#' Statistics in Biopharmaceutical Research, 8(3), 345-354, <doi:10.1080/19466315.2016.1174149>.
#'
#'@examples
#'\dontrun{
#'#derive the predictive distribution of the DLT rate of the dose 10+24
#'#from historical data of three previous trials
#'pred <- get_pred_tox(
#'  doses.of.interest = c(10, 24),
#'  dose1 = c(2, 4, 8, 10, 14, 0,  0,  0,  4,  8,  10, 12),
#'  dose2 = c(0, 0, 0, 0,  0,  12, 24, 36, 24, 24, 24, 24),
#'  n.pat = c(3, 3, 3, 9,  6,  3,  12, 3,  3,  3,  9,  3),
#'  n.dlt = c(0, 0, 0, 0,  2,  0,  1,  0,  0,  0,  1,  1),
#'  trial = c(1, 1, 1, 1,  1,  2,  2,  2,  3,  3,  3,  3),
#'  dose.ref1 = 14,
#'  dose.ref2 = 24,
#'  return.samples = TRUE
#')
#'
#'#summary
#'pred$summary
#'
#'#first eight samples from predictive DLT rate
#'pred$samples[, 1:8]
#'
#'}
#'@export

get_pred_tox <- function(
  doses.of.interest,
  dose1,
  dose2,
  dose.ref1,
  dose.ref2,
  n.pat,
  n.dlt,
  trial,
  prior.mu = list(mu_a1 =  c(logit(0.33), 2),
                  mu_b1 =  c(0,          1),
                  mu_a2 =  c(logit(0.33), 2),
                  mu_b2 =  c(0,          1),
                  mu_eta = c(0,          1.121)
  ),
  prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
                   tau_b1 =  c(log(0.125), log(2)/1.96),
                   tau_a2 =  c(log(0.25),  log(2)/1.96),
                   tau_b2 =  c(log(0.125), log(2)/1.96),
                   tau_eta = c(log(0.125), log(2)/1.96)
  ),
  saturating = FALSE,
  return.samples = FALSE,
  probs = c(0.025, 0.5, 0.975),
  dosing.intervals = c(0.16, 0.33),
  ...
){
  if(!is.logical(saturating)){
    stop("`saturating` must be logical.")
  }
  if(!is.numeric(doses.of.interest)){
    stop("`doses.of.interest` must be numeric.")
  }
  if(!length(dim(doses.of.interest))==2){
    if(length(doses.of.interest)==2){
      doses.of.interest <- rbind(doses.of.interest[1], doses.of.interest[2])
    }else{
      stop("`doses.of.interest` must be a matrix or a vector of length 2.")
    }
  }
  if(!dim(doses.of.interest)[1]==2){
    stop("`doses.of.interest` must be a matrix with 2 rows or a vector of length 2.")
  }
  doses.of.interest[1, which(is.na(doses.of.interest[1, ]))] <- rep(0, times = length(which(is.na(doses.of.interest[1, ]))))
  doses.of.interest[2, which(is.na(doses.of.interest[2, ]))] <- rep(0, times = length(which(is.na(doses.of.interest[2, ]))))

  if(!is.num(doses.of.interest, low=0, lB=T, uB=F)){
    stop("`doses.of.interest` must contain non-negative (or NA) entries.")
  }

  if(missing(dose1) & !missing(dose2)){
    dose1 <- rep(0, times = length(dose2))
    if(max(doses.of.interest[1,])==0 & missing(dose.ref1)){
      dose.ref1 <- 1
    }else if(missing(dose.ref1)){
      stop("`doses.of.interest` contain mono 1 doses, so `dose.ref1` must be given.")
    }
  }

  if(missing(dose2) & !missing(dose1)){
    dose2 <- rep(0, times = length(dose1))
    if(max(doses.of.interest[2,])==0 & missing(dose.ref2)){
      dose.ref2 <- 1
    }else if(missing(dose.ref2)){
      stop("`doses.of.interest` contain mono 2 doses, so `dose.ref2` must be given.")
    }
  }

  if(missing(dose1) & missing(dose2)){
    stop("Either `dose1` or `dose2` or both must be specified")
  }

  if(!is.num(dosing.intervals, low=0, up=1, lB=F, uB=F)){
    stop("`dosing.intervals` must contain non-NA numerical entries strictly between 0 and 1")
  }

  if(is.unsorted(dosing.intervals, strictly = T)){
    stop("`dosing.intervals` must have strictly increasing entries")
  }

  if(!length(dosing.intervals)>0){
    stop("`dosing.intervals` must at least have length 1.")
  }



  #extract samples from MAP prior
  smpl <- get_MAP(
    dose1 = dose1,
    dose2 = dose2,
    dose.ref1 = dose.ref1,
    dose.ref2 = dose.ref2,
    n.pat = n.pat,
    n.dlt = n.dlt,
    trial = trial,
    prior.mu = prior.mu,
    prior.tau = prior.tau,
    return.samples = T,
    probs = probs,
    saturating = saturating,
    ...
  )$samples

  #transform log(beta) to beta itself
  smpl[2,] <- exp(smpl[2,])
  smpl[4,] <- exp(smpl[4,])

  #pre-transform doses
  doset <- doses.of.interest
  logdoset <- doses.of.interest
  doset[1,] <- doset[1,]/dose.ref1
  logdoset[1,] <- log(doset[1,])
  doset[2,] <- doset[2,]/dose.ref2
  logdoset[2,] <- log(doset[2,])

  #matrix for samples from predictive tox
  nsmp <- length(smpl[1,])
  ndos <- length(doses.of.interest[1,])
  smpl_tox <- matrix(NA, nrow=ndos, ncol=nsmp,
                     dimnames=list("dose"=paste0(doses.of.interest[1,], "+",
                                                 doses.of.interest[2,]),
                                  "iteration" = 1:nsmp))

  #intermediate results for combi treatment
  t0 <- rep(0, nsmp)
  t1 <- rep(0, nsmp)
  t2 <- rep(0, nsmp)

  #compute samples from predictive tox for all doses of interest
  for(i in 1:ndos){
    if(doses.of.interest[1, i]<=0 & doses.of.interest[2, i]<=0){
      #no dose, fill in 0 tox
      smpl_tox[i, ] <- rep(0, times = nsmp)
    }else if(doses.of.interest[1,i]>0 & doses.of.interest[2,i]<=0){
      #treatment is mono1
      smpl_tox[i, ] <- inv_logit(smpl[1,] + smpl[2, ]*logdoset[1,i])
    }else if(doses.of.interest[2,i]>0 & doses.of.interest[1,i]<=0){
      #treatment is mono2
      smpl_tox[i, ] <- inv_logit(smpl[3,] + smpl[4, ]*logdoset[2,i])
    }else{
      #treatment is combi
      t1 <- inv_logit(smpl[1,] + smpl[2, ]*logdoset[1,i])
      t2 <- inv_logit(smpl[3,] + smpl[4, ]*logdoset[2,i])
      t0 <- t1 + t2 - t1*t2
      smpl_tox[i, ] <- inv_logit(logit(t0) + smpl[5,]*doset[1,i]*doset[2,i])
    }
  }

  #number of quantiles and intervals
  nquan <- length(probs)
  nint <- length(dosing.intervals)+1

  namesint <- rep("", times = nint)
  for(j in 1:nint){
    if(j==1){
      namesint[j] <- paste0("P([0,", dosing.intervals[j], "))")
    }else if(j==nint){
      namesint[j] <- paste0("P([", dosing.intervals[j-1], ",1])")
    }else{
      namesint[j] <- paste0("P([",dosing.intervals[j-1], "," , dosing.intervals[j], "))")
    }
  }

  #prepare summary matrix
  summ <- matrix(NA, nrow = ndos, ncol = 2 + nquan + nint,
                 dimnames = list("dose" = paste0(doses.of.interest[1,], "+",
                                                        doses.of.interest[2,]),
                                 "summary" = c("mean", "SD",
                                        paste0(round(100*probs, digits = 3), "%"),
                                        namesint)))
  for(i in 1:ndos){
    summ[i,1] <- mean(smpl_tox[i, ])
    summ[i,2] <- sd(smpl_tox[i, ])
    summ[i,3:(2+nquan)] <- quantile(smpl_tox[i,], probs = probs)
    summ[i, (2+nquan+1):(2+nquan+nint)] <- get_int_probs_BLRM(smpl_tox[i,], dosing.intervals)
  }

  if(!return.samples){
    return(summ)
  }else{
    return(list(
      summary = summ,
      samples = smpl_tox
    ))
  }

}




#-------------------------------------------------------------------------------
#Monotherapy variant
#-------------------------------------------------------------------------------
#'@rdname get_pred_tox
#'@param dose Numeric vector with non-negative entries that defines the monotherapy dose
#'levels administered to each cohort.
#'@param dose.ref Positive numerical, reference dose for monotherapy.
#'@param prior.mu.mono Optional named list that specifies the prior. Must have the
#' entries \code{mu_a} and \code{mu_b}, both of which must be vectors of length 2,
#' where the second entry is a positive number.\cr
#' The entries need to be:
#' \itemize{
#' \item{\code{...$mu_a}:\cr
#' Defaults to \code{c(logit(0.33), 2)}, provides mean and SD of the normal distribution used as
#' hypermean of \eqn{log(\alpha)}.}
#' \item{\code{...$mu_b}:\cr
#' Defaults to \code{c(0, 1)}, provides mean and SD of the normal distribution used as
#' hypermean of \eqn{log(\beta)}.}
#' }
#' @param prior.tau.mono Optional named list that specifies the prior. Must have the
#' entries \code{tau_a} and \code{tau_b}, both of which must be vectors of length 2,
#' where the second entry is a positive number.\cr
#' The entries need to be:
#' \itemize{
#' \item{\code{...$tau_a}:\cr
#' Defaults to \code{c(log(0.25), log(2)/1.96)}, provides mean and SD (on log-scale) of the log-normal distribution used as
#' between-trial heterogeneity of \eqn{log(\alpha)}.}
#' \item{\code{...$tau_b}:\cr
#' Defaults to \code{c(log(0.125), log(2)/1.96)}, provides mean and SD (on log-scale) of the log-normal distribution used as
#' between-trial heterogeneity of \eqn{log(\beta)}.}
#' }
#'@usage
#'get_pred_tox_mono(
#'  doses.of.interest,
#'  dose,
#'  dose.ref,
#'  n.pat,
#'  n.dlt,
#'  trial,
#'  prior.mu.mono = list(mu_a =  c(logit(0.33), 2),
#'                       mu_b =  c(0,          1)),
#'  prior.tau.mono = list(tau_a =  c(log(0.25),  log(2)/1.96),
#'                        tau_b =  c(log(0.125), log(2)/1.96)),
#'  return.samples = FALSE,
#'  probs = c(0.025, 0.5, 0.975),
#'  dosing.intervals = c(0.16, 0.33),
#'  ...
#')
#'@examples
#'\dontrun{
#'#derive the predictive distribution of the DLT rates of the doses 8, 12, and
#'#14 (monotherapy) from historical data of two previous trials
#'pred_mono <- get_pred_tox_mono(
#'  doses.of.interest = c(8, 12, 14),
#'  dose  = c(2, 4, 8, 10, 14, 10, 14),
#'  n.pat = c(3, 3, 3, 9,  6,  6,  12),
#'  n.dlt = c(0, 0, 0, 0,  2,  0,  1),
#'  trial = c(1, 1, 1, 1,  1,  2,  2),
#'  dose.ref = 14,
#'  return.samples = TRUE
#')
#'
#'#summary
#'pred_mono$summary
#'
#'#first eight samples from predictive DLT rates
#'pred_mono$samples[, 1:8]
#'
#'}
#'@export

get_pred_tox_mono <- function(
  doses.of.interest,
  dose,
  dose.ref,
  n.pat,
  n.dlt,
  trial,
  prior.mu.mono = list(mu_a =  c(logit(0.33), 2),
                       mu_b =  c(0,          1)),
  prior.tau.mono = list(tau_a =  c(log(0.25),  log(2)/1.96),
                        tau_b =  c(log(0.125), log(2)/1.96)),
  return.samples = FALSE,
  probs = c(0.025, 0.5, 0.975),
  dosing.intervals = c(0.16, 0.33),
  ...
){
  if(!is.numeric(doses.of.interest)){
    stop("`doses.of.interest` must be numeric.")
  }

  if(!length(dim(doses.of.interest))%in%c(0, 2)){
    stop("For `get_pred_tox_mono`, `doses.of.interest` must be a vector or a matrix either with one row or one column.")
  }else if(length(dim(doses.of.interest))==2){
    if(!dim(doses.of.interest)[1]==1 & !dim(doses.of.interest)[2]==1){
      stop("The argument `doses.of.interest` was given as a matrix. For `get_pred_tox_mono()`, it must in this case either have only one row or only one column. Note that vectors are also allowed.")
    }else if(dim(doses.of.interest)[1]==1){
      doses.of.interest <- doses.of.interest[1,]
    }else{
      doses.of.interest <- doses.of.interest[,1]
    }
  }

  if(!is.num(doses.of.interest, low=0, lB=T, uB=F)){
    stop("`doses.of.interest` must contain non-negative (or NA) entries.")
  }

  if(!is.num(dosing.intervals, low=0, up=1, lB=F, uB=F)){
    stop("`dosing.intervals` must contain non-NA numerical entries strictly between 0 and 1")
  }

  if(is.unsorted(dosing.intervals, strictly = T)){
    stop("`dosing.intervals` must have strictly increasing entries")
  }

  if(!length(dosing.intervals)>0){
    stop("`dosing.intervals` must at least have length 1.")
  }



  #extract samples from MAP prior
  smpl <- get_MAP_mono(
    dose = dose,
    dose.ref = dose.ref,
    n.pat = n.pat,
    n.dlt = n.dlt,
    trial = trial,
    prior.mu.mono = prior.mu.mono,
    prior.tau.mono = prior.tau.mono,
    return.samples = T,
    probs = probs,
    ...
  )$samples

  #transform log(beta) to beta itself
  smpl[2,] <- exp(smpl[2,])

  #pre-transform doses
  doset <- doses.of.interest/dose.ref
  logdoset <- log(doset)

  #matrix for samples from predictive tox
  nsmp <- length(smpl[1,])
  ndos <- length(doses.of.interest)
  smpl_tox <- matrix(NA, nrow=ndos, ncol=nsmp,
                     dimnames=list("dose"=doses.of.interest,
                                   "iteration" = 1:nsmp))


  #compute samples from predictive tox for all doses of interest
  for(i in 1:ndos){
    if(doses.of.interest[i]<=0){
      #no dose, fill in 0 tox
      smpl_tox[i, ] <- rep(0, times = nsmp)
    }else{
      #treatment is mono1
      smpl_tox[i, ] <- inv_logit(smpl[1,] + smpl[2, ]*logdoset[i])
    }
  }

  #number of quantiles and intervals
  nquan <- length(probs)
  nint <- length(dosing.intervals)+1

  namesint <- rep("", times = nint)
  for(j in 1:nint){
    if(j==1){
      namesint[j] <- paste0("P([0,", dosing.intervals[j], "))")
    }else if(j==nint){
      namesint[j] <- paste0("P([", dosing.intervals[j-1], ",1])")
    }else{
      namesint[j] <- paste0("P([",dosing.intervals[j-1], "," , dosing.intervals[j], "))")
    }
  }

  #prepare summary matrix
  summ <- matrix(NA, nrow = ndos, ncol = 2 + nquan + nint,
                 dimnames = list("dose" = doses.of.interest,
                                 "summary" = c("mean", "SD",
                                               paste0(round(100*probs, digits = 3), "%"),
                                               namesint)))
  for(i in 1:ndos){
    summ[i,1] <- mean(smpl_tox[i, ])
    summ[i,2] <- sd(smpl_tox[i, ])
    summ[i,3:(2+nquan)] <- quantile(smpl_tox[i,], probs = probs)
    summ[i, (2+nquan+1):(2+nquan+nint)] <- get_int_probs_BLRM(smpl_tox[i,], dosing.intervals)
  }

  if(!return.samples){
    return(summ)
  }else{
    return(list(
      summary = summ,
      samples = smpl_tox
    ))
  }

}

