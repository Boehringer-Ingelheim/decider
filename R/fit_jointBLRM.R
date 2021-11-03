

#------------------------------------------------------------------------------
#Function to compute a fitted stanmodel
#------------------------------------------------------------------------------

#' Fit the joint BLRM and extract a stanfit object
#' @description Call the \code{\link[rstan:sampling]{rstan::sampling}()} method for
#' the joint BLRM to extract a \code{\link[rstan:stanfit-class]{rstan::stanfit-class}} object.
#' Refer to \code{\link[OncoBLRM:scenario_jointBLRM]{scenario_jointBLRM}()} for a
#' documentation of the underlying model.
#' @param dose1 Numerical vector that lists the dose levels of the first compound applied to each
#' patient. The length must be equal to the number of cohorts. If the first compound was not
#' administered to a cohort, use 0 for this entry.
#' @param dose2 Numerical vector that lists the dose levels of the second compound applied to each
#' patient. The length must be equal to the number of cohorts. If the second compound was not
#' administered to a cohort, use 0 for this entry.
#' @param dose.ref1 Positive number, reference dose for compound 1.
#' @param dose.ref2 Positive number, reference dose for compound 2.
#' @param n.pat Numerical vector that lists the number of patients for each cohort.
#' @param n.dlt Numerical vector that lists the number of DLTs for each cohort.
#' @param trial Numerical vector that indicates to which trial a cohort is assigned to.
#' The entries can either be numericals or strings, where \code{1} and \code{"1"}
#' would be interpreted as the same trial. The entries of \code{trial} are internally
#' converted to numbers from 1 to the number of different studies if not already
#' given in this form.
#' @param MAP.prior Optional logical, defaults to \code{FALSE}. If \code{TRUE}, the returned stanfit contains an additional
#' trial number that provides the resulting meta-analytic predictive prior for a new trial.
#' @param prior.mu Optional list that specifies the prior.
#' Same format as documented in \code{\link[OncoBLRM:scenario_jointBLRM]{scenario_jointBLRM}()}.
#' @param prior.tau Optional list that specifies the prior.
#' Same format as documented in \code{\link[OncoBLRM:scenario_jointBLRM]{scenario_jointBLRM}()}.
#' @param saturating Optional logical, defaults to \code{FALSE}. If \code{TRUE}, the BLRM will be using a saturating interaction term as described in
#' \code{\link[OncoBayes2:blrm_formula_saturating]{OncoBayes2::blrm_formula_saturating}()}. Also refer to the Details section in the documentation
#' of \code{\link[OncoBLRM:scenario_jointBLRM]{scenario_jointBLRM}()}.
#' @param ... Additional optional arguments that are passed to \code{\link[rstan:sampling]{rstan::sampling}()}. These
#' include the settings for MCMC and the list of parameters which are returned in the stanfit object. Refer
#' to the \code{\link[rstan:rstan-package]{rstan-package}} for an overview.
#' @returns An object of class \code{\link[rstan:stanfit-class]{rstan::stanfit-class}} from the
#' \code{\link[rstan:rstan-package]{rstan-package}}.
#'
#' @details Refer to \code{\link[OncoBLRM:scenario_jointBLRM]{scenario_jointBLRM}()} for
#' a description of the underlying model.
#'
#' @seealso \code{\link[OncoBLRM:scenario_jointBLRM]{scenario_jointBLRM}()},
#' \code{\link[OncoBLRM:get_pred_tox]{get_pred_tox}()}, \code{\link[OncoBLRM:get_MAP]{get_MAP}()},
#'  \code{\link[rstan:sampling]{rstan::sampling}()},
#' \code{\link[rstan:stanfit-class]{rstan::stanfit-class}}, \code{\link[rstan:rstan-package]{rstan-package}}.
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. \url{https://mc-stan.org}.
#'
#' Neuenschwander, B., Branson, M., & Gsponer, T. (2008). Critical aspects of the Bayesian approach to phase I cancer trials.
#' Statistics in medicine, 27(13), 2420-2439, <doi:10.1002/sim.3230>.
#'
#' Neuenschwander, B., Matano, A., Tang, Z., Roychoudhury, S., Wandel, S., & Bailey, S. (2014). A Bayesian Industry Approach to Phase I Combination Trials in Oncology.
#' In: Zhao. W & Yang, H. (editors). Statistical methods in drug combination studies. Chapman and Hall/CRC, 95-135, <doi:10.1201/b17965>.
#'
#' Neuenschwander, B., Roychoudhury, S., & Schmidli, H. (2016). On the use of co-data in clinical trials.
#' Statistics in Biopharmaceutical Research, 8(3), 345-354, <doi:10.1080/19466315.2016.1174149>.
#'
#' @usage
#' fit_jointBLRM(
#'    dose1,
#'    dose2,
#'    dose.ref1,
#'    dose.ref2,
#'    n.pat,
#'    n.dlt,
#'    trial,
#'    MAP.prior = FALSE,
#'    saturating = FALSE,
#'    prior.mu = list(mu_a1 =  c(logit(0.33), 2),
#'                    mu_b1 =  c(0,          1),
#'                    mu_a2 =  c(logit(0.33), 2),
#'                    mu_b2 =  c(0,          1),
#'                    mu_eta = c(0,          1.121)),
#'    prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
#'                     tau_b1 =  c(log(0.125), log(2)/1.96),
#'                     tau_a2 =  c(log(0.25),  log(2)/1.96),
#'                     tau_b2 =  c(log(0.125), log(2)/1.96),
#'                     tau_eta = c(log(0.125), log(2)/1.96)),
#'    ...)
#' @examples
#' \dontrun{
#'   fit <- fit_jointBLRM(
#'               dose1 = c(1, 2, 4, 8, 12, 0,  0,  0,  2,  4,  6),
#'               dose2 = c(0, 0, 0, 0, 0,  10, 20, 30, 10, 10, 20),
#'               dose.ref1 = 12,
#'               dose.ref2 = 30,
#'               n.pat = c(3, 3, 3, 3, 9,  3,  3,  6,  3,  3,  3),
#'               n.dlt = c(0, 0, 0, 0, 2,  0,  0,  1,  0,  0,  1),
#'               trial = c(1, 1, 1, 1, 1,  2,  2,  2,  3,  3,  3),
#'
#'               iter = 10000,
#'               chains = 4,
#'               control = list(adapt_delta = 0.95),
#'               pars = c("log_ab")
#'  )
#'
#'  #print posterior summary
#'  print(fit)
#'
#'  #posterior density plots for parameters of trial 3.
#'  plot <- rstan::stan_plot(fit, pars=paste0("log_ab[3,", 1:5, "]"),
#'                           show_density=TRUE)
#'  print(plot)
#'
#'  #extract samples in matrix format
#'  samples <- as.matrix(fit, pars="log_ab")
#'  head(samples)
#' }
#' @export
fit_jointBLRM <- function(
  dose1,
  dose2,
  dose.ref1,
  dose.ref2,
  n.pat,
  n.dlt,
  trial,
  MAP.prior = FALSE,
  saturating = FALSE,
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
  ...
){
  nobs <- length(dose1)
  if(! (nobs==length(dose2)
        & nobs==length(n.pat)
        & nobs==length(n.dlt)
        & nobs==length(trial)
  )
  ){
    stop("dose1, dose2, n.pat, n.dlt, and trial must have the same length.")
  }

  if(!is.logical(MAP.prior)){
    stop("`MAP.prior` must be logical.")
  }
  if(!is.logical(saturating)){
    stop("`saturating` must be logical.")
  }
  if(!is.num(dose1, low=0, up= Inf, len=nobs, lB=T, uB=F)){
    stop("`dose1` must be a numerical with positive entries.")
  }
  if(!is.num(dose2, low=0, up= Inf, len=nobs, lB=T, uB=F)){
    stop("`dose2` must be a numerical with positive entries.")
  }
  if(!is.num(dose.ref1, low=0, up= Inf, len=1, lB=T, uB=F)){
    stop("`dose.ref1` must be a positive single number.")
  }
  if(!is.num(dose.ref2, low=0, up= Inf, len=1, lB=T, uB=F)){
    stop("`dose.ref2` must be a positive single number.")
  }

  if(!is.wholenumber(n.pat)){
    stop("`n.pat` must be an vector with integer entries.")
  }

  if(!is.wholenumber(n.dlt)){
    stop("`n.dlt` must be an vector with integer entries.")
  }

  if(!all(n.pat>=0)){
    stop("`n.pat` must have non-negative entries.")
  }
  if(!(all(n.dlt>=0) & all(n.dlt <=n.pat))){
    stop("`n.dlt` must have non-negative entries that are smaller than the number\n",
         "of patients of the corresponding cohorts.")
  }

  if(!is.character(trial) & !is.num(trial,uB=F, lB=F)){
    stop("`trial` must contain either character or numeric entries that indicate\n",
         "the corresponding trial of the cohorts in the data.")
  }

  nstd <- length(unique(trial))
  stdnames <- unique(trial)
  for(i in 1:nstd){
    trial[which(trial==stdnames[i])] <- rep(i, times=length(which(trial==stdnames[i])))
  }

  #new short prior checks
  if(!is.list(prior.mu)){
    stop("`prior.mu` must be a named list.")
  }
  if(is.null(names(prior.mu))){
    stop("`prior.mu` must be a named list.")
  }
  names(prior.mu) <- tolower(names(prior.mu))
  test.prior.mu(prior.mu)

  if(!is.list(prior.tau)){
    stop("`prior.tau` must be a named list.")
  }
  if(is.null(names(prior.tau))){
    stop("prior.tau must be a named list.")
  }
  names(prior.tau) <- tolower(names(prior.tau))
  test.prior.tau(prior.tau)


  return(sampling(
    object = stanmodels$jointBLRM,
    data = list(
      n_obs = nobs,
      n_studies = nstd,
      dose_1 = as.array(dose1),
      dose_2 = as.array(dose2),
      dose_c = c(dose.ref1, dose.ref2),
      n = as.array(n.pat),
      r = as.array(n.dlt),
      s = as.array(trial),
      doMAP = ifelse(MAP.prior, 1, 0),
      saturating = ifelse(saturating, 1, 0),
      mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
                  prior.mu$mu_a2[1], prior.mu$mu_b2[1], prior.mu$mu_eta[1]),
      sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
                prior.mu$mu_a2[2], prior.mu$mu_b2[2], prior.mu$mu_eta[2]),
      mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
                   prior.tau$tau_a2[1], prior.tau$tau_b2[1], prior.tau$tau_eta[1]),
      sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
                 prior.tau$tau_a2[2], prior.tau$tau_b2[2], prior.tau$tau_eta[2])
    ),
    ...))
}







