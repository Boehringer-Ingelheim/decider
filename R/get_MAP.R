#' Compute MAP priors for BLRMs
#' @name get_MAP
#' @rdname get_MAP
#' @aliases get_MAP_mono get_MAP
#' @description Fit data from one or more monotherapy or combination therapy trials
#' to a hierarchical BLRM and extract estimates of the resulting meta-analytic
#' predictive prior for a new monotherapy or combination therapy trial.
#'
#' The function \code{get_MAP_mono()} supports only monotherapy data (from one or more trials)
#' as input and will produce the MAP prior for a two-parameter monotherapy BLRM.\cr
#' For the function \code{get_MAP()}, data from one or more monotherapy and/or
#' combination therapy trials can be given, which is then used to derive the MAP
#' prior for a five-parameter combination therapy BLRM. Note that this also contains
#' the parameter for the implicit model for the monotherapy DLT rates, so that one
#' can use \code{get_MAP()} in principle also to derive MAP priors for monotherapy
#' from both monotherapy and combination therapy data.
#' @param dose1 Numerical vector of non-negative entries that provide the
#' administered dose level of compound 1 for each cohort. To include monotherapy
#' dose levels in the data given to \code{get_MAP()}, one can use \code{0}
#' or \code{NA} to indicate that compound 1 was not administered to a cohort.
#' Note that cohorts for which both \code{dose1} and \code{dose2} are either
#' \code{0} or \code{NA} will be removed from the data.
#' @param dose2 Numerical vector of non-negative entries that provide the
#' administered dose level of compound 2 for each cohort. To include monotherapy
#' dose levels in the data given to \code{get_MAP()}, one can use \code{0}
#' or \code{NA} to indicate that compound 2 was not administered to a cohort.
#' Note that cohorts for which both \code{dose1} and \code{dose2} are either
#' \code{0} or \code{NA} will be removed from the data.
#' @param dose.ref1 Positive number, reference dose for compound 1.
#' @param dose.ref2 Positive number, reference dose for compound 2.
#' @param n.pat Numerical vector that lists the number of patients for each cohort.
#' @param n.dlt Numerical vector that lists the number of DLTs for each cohort.
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
#' @param return.samples Optional logical, defaults to \code{FALSE}. If \code{TRUE},
#' the function will return a matrix with samples from the MAP prior additionally
#' to a prior summary.
#' @param probs Optional numerical, defaults to \code{c(0.025, 0.5, 0.975)}.
#' Indicates which quantiles should be included in the posterior summaries.
#' @param ... Optional additional arguments that are passed to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()}. These
#' include the settings for MCMC and the random seed. Refer
#' to the \code{\link[rstan:rstan]{rstan-package}} for an overview. Note that
#' the function will automatically use \code{iter=10000} and a \code{control} argument
#' containing \code{adapt_delta=0.9} to ensure a relatively good fit if the
#' arguments \code{iter} and \code{control} for \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()}
#' are not provided manually in \code{...}.
#' @returns A table containing prior summaries of the MAP priors for either the parameters
#' \eqn{log(\alpha)} and \eqn{log(\beta)} (function \code{get_MAP_mono()}) or of the
#' parameters \eqn{log(\alpha_1)}, \eqn{log(\beta_1)}, \eqn{log(\alpha_2)}, \eqn{log(\beta_2)},
#' and \eqn{\eta} (function \code{get_MAP()}). If \code{return.samples} is \code{TRUE},
#' the output is instead a list with entries \code{...$summary} and \code{...$samples},
#' where the first entry contains the prior summary and the second entry a table with
#' the drawn samples for each parameter.
#' @details
#' The model for computing meta-analytic predictive prior is specified according to
#' Neuenschwander et al. (2014, 2016). A model description of the hierarchical BLRM
#' used by the function to obtain the predictive distribution can be found in the documentation
#' of \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()}.
#'
#' The underlying model for monotherapy MAP priors is the two-parameter monotherapy variant of the
#' more general hierarchical model that is described in the documentation of
#' \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()}.
#' As the function \code{get_MAP_mono()} only considers monotherapy, only the parameters \eqn{log(\alpha)}
#' and \eqn{log(\beta)} are needed to model the DLT rate of each dose level.
#' Besides this, the same hierarchical model structure as described in
#' \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} is used.
#' @seealso  \code{\link[decider:get_pred_tox]{get_pred_tox}()},
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
#'
#' @usage
#' get_MAP(
#'   dose1,
#'   dose2,
#'   dose.ref1,
#'   dose.ref2,
#'   n.pat,
#'   n.dlt,
#'   trial,
#'   prior.mu = list(mu_a1 =  c(logit(0.33), 2),
#'                   mu_b1 =  c(0,           1),
#'                   mu_a2 =  c(logit(0.33), 2),
#'                   mu_b2 =  c(0,           1),
#'                   mu_eta = c(0,           1.121)
#'   ),
#'   prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
#'                    tau_b1 =  c(log(0.125), log(2)/1.96),
#'                    tau_a2 =  c(log(0.25),  log(2)/1.96),
#'                    tau_b2 =  c(log(0.125), log(2)/1.96),
#'                    tau_eta = c(log(0.125), log(2)/1.96)
#'   ),
#'   saturating = FALSE,
#'   return.samples = FALSE,
#'   probs = c(0.025, 0.5, 0.975),
#'   ...
#' )
#' @examples
#' \dontrun{
#' MAPprior_combi <- get_MAP(
#'   dose1 = c(1, 2, 3, 4,   0, 0,   2, 2),
#'   dose2 = c(0, 0, 0, 0,   5, 8,   5, 8),
#'   n.pat = c(3, 3, 3, 6,   3, 9,   3, 6),
#'   n.dlt = c(0, 0, 0, 1,   0, 0,   0, 0),
#'   trial = c(1, 1, 1, 1,   2, 2,   3, 3),
#'   dose.ref1 = 5,
#'   dose.ref2 = 8,
#'   return.samples = TRUE
#' )
#'
#' MAPprior_combi$summary
#' MAPprior_combi$samples[, 1:8]
#' }
#' @export

get_MAP <- function(
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
  probs = c(0.025,  0.5, 0.975),
  ...
){
  nobs <- length(dose1)
  if(!nobs>=1){
    stop("`dose1` must have at least one entry.")
  }
  if(! (nobs==length(n.pat)
        & nobs==length(n.dlt)
        & nobs==length(trial)
        & nobs==length(dose2)
  )
  ){
    stop("`dose1`, `dose2`, `n.pat`, `n.dlt`, and `trial` must have the same length.")
  }

  if(!is.logical(return.samples)){
    stop("`return.samples` must be logical.")
  }
  if(!is.logical(saturating)){
    stop("`saturating` must be logical.")
  }

  dose1[which(is.na(dose1))] <- rep(0, times= length(which(is.na(dose1))))
  dose2[which(is.na(dose2))] <- rep(0, times= length(which(is.na(dose2))))

  if(!is.num(dose1, low=0, up= Inf, len=nobs, lB=T, uB=F)){
    stop("`dose1` must be a numerical with positive entries.")
  }
  if(!is.num(dose.ref1, low=0, up= Inf, len=1, lB=F, uB=F)){
    stop("`dose.ref1` must be a positive single number.")
  }
  if(!is.num(dose2, low=0, up= Inf, len=nobs, lB=T, uB=F)){
    stop("`dose2` must be a numerical with positive entries.")
  }
  if(!is.num(dose.ref2, low=0, up= Inf, len=1, lB=F, uB=F)){
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


  #remove observations which do not contribute.
  idx_noninf_obs <- which((dose1==0 & dose2==0))

  if(length(idx_noninf_obs)>0){
    if(length(idx_noninf_obs)<nobs){
      message("Note: Detected and removed cohorts for which both doses are 0 or NA.")
      allid <- 1:nobs
      keep <- allid[which(!allid%in%idx_noninf_obs)]
      dose1 <- dose1[keep]
      dose2 <- dose2[keep]
      n.pat <- n.pat[keep]
      n.dlt <- n.dlt[keep]
      trial <- trial[keep]
    }else if (length(idx_noninf_obs)==nobs){
      message("Note: None of the given cohorts were administered one of the\n",
              "dose levels. Sampling without data from prior.")
      dose1 <- c(dose.ref1)
      dose2 <- c(0)
      n.pat <- c(0)
      n.dlt <- c(0)
      trial <- c(1)

    }
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
    stop("`prior.tau` must be a named list.")
  }
  names(prior.tau) <- tolower(names(prior.tau))
  test.prior.tau(prior.tau)

  if(!is.num(probs, low=0, up=1)){
    stop("`probs` must be a numeric vector with entries between 0 and 1.")
  }
  if(!length(probs)>0){
    stop("`probs` must be at least of length 1.")
  }

  dots <- list(...)
  if("refresh"%in% names(dots)){
    if("iter"%in%names(dots) & "control"%in%names(dots)){
    smpl <- as.matrix(sampling(
      object = stanmodels$jointBLRM,
      data = list(
        n_obs = length(dose1),
        n_studies = nstd,
        dose_1 = as.array(dose1),
        dose_2 = as.array(dose2),
        dose_c = c(dose.ref1, dose.ref2),
        n = as.array(n.pat),
        r = as.array(n.dlt),
        s = as.array(trial),
        doMAP = 1,
        saturating = ifelse(saturating, 1, 0),
        mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
                    prior.mu$mu_a2[1], prior.mu$mu_b2[1],prior.mu$mu_eta[1]),
        sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
                  prior.mu$mu_a2[2], prior.mu$mu_b2[2],prior.mu$mu_eta[2]),
        mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
                     prior.tau$tau_a2[1], prior.tau$tau_b2[1],prior.tau$tau_eta[1]),
        sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
                   prior.tau$tau_a2[2], prior.tau$tau_b2[2],prior.tau$tau_eta[2])
      ),
      ...),
      pars = paste0("log_ab[",nstd+1, ",", 1:5 ,"]"))

  }else if((!"iter"%in%names(dots)) & "control"%in%names(dots)){

    smpl <- as.matrix(sampling(
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
        doMAP = 1,
        saturating = ifelse(saturating, 1, 0),
        mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
                    prior.mu$mu_a2[1], prior.mu$mu_b2[1],prior.mu$mu_eta[1]),
        sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
                  prior.mu$mu_a2[2], prior.mu$mu_b2[2],prior.mu$mu_eta[2]),
        mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
                     prior.tau$tau_a2[1], prior.tau$tau_b2[1],prior.tau$tau_eta[1]),
        sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
                   prior.tau$tau_a2[2], prior.tau$tau_b2[2],prior.tau$tau_eta[2])
      ),
      iter = 10000,
      ...),
      pars = paste0("log_ab[",nstd+1, ",", 1:5 ,"]"))


  }else if((!"control"%in%names(dots)) & "iter"%in%names(dots)){


    smpl <- as.matrix(sampling(
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
        doMAP = 1,
        saturating = ifelse(saturating, 1, 0),
        mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
                    prior.mu$mu_a2[1], prior.mu$mu_b2[1],prior.mu$mu_eta[1]),
        sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
                  prior.mu$mu_a2[2], prior.mu$mu_b2[2],prior.mu$mu_eta[2]),
        mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
                     prior.tau$tau_a2[1], prior.tau$tau_b2[1],prior.tau$tau_eta[1]),
        sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
                   prior.tau$tau_a2[2], prior.tau$tau_b2[2],prior.tau$tau_eta[2])
      ),
      control = list(adapt_delta=0.9),
      ...),
      pars = paste0("log_ab[",nstd+1, ",", 1:5 ,"]"))


  }else {

    smpl <- as.matrix(sampling(
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
        doMAP = 1,
        saturating = ifelse(saturating, 1, 0),
        mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
                    prior.mu$mu_a2[1], prior.mu$mu_b2[1],prior.mu$mu_eta[1]),
        sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
                  prior.mu$mu_a2[2], prior.mu$mu_b2[2],prior.mu$mu_eta[2]),
        mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
                     prior.tau$tau_a2[1], prior.tau$tau_b2[1],prior.tau$tau_eta[1]),
        sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
                   prior.tau$tau_a2[2], prior.tau$tau_b2[2],prior.tau$tau_eta[2])
      ),
      iter = 10000,
      control = list(adapt_delta=0.9),
      ...),
      pars = paste0("log_ab[",nstd+1, ",", 1:5 ,"]"))

  }
  }else{
    #refresh not in ...

    if("iter"%in%names(dots) & "control"%in%names(dots)){
      smpl <- as.matrix(sampling(
        object = stanmodels$jointBLRM,
        data = list(
          n_obs = length(dose1),
          n_studies = nstd,
          dose_1 = as.array(dose1),
          dose_2 = as.array(dose2),
          dose_c = c(dose.ref1, dose.ref2),
          n = as.array(n.pat),
          r = as.array(n.dlt),
          s = as.array(trial),
          doMAP = 1,
          saturating = ifelse(saturating, 1, 0),
          mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
                      prior.mu$mu_a2[1], prior.mu$mu_b2[1],prior.mu$mu_eta[1]),
          sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
                    prior.mu$mu_a2[2], prior.mu$mu_b2[2],prior.mu$mu_eta[2]),
          mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
                       prior.tau$tau_a2[1], prior.tau$tau_b2[1],prior.tau$tau_eta[1]),
          sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
                     prior.tau$tau_a2[2], prior.tau$tau_b2[2],prior.tau$tau_eta[2])
        ),
        refresh = 0,
        ...),
        pars = paste0("log_ab[",nstd+1, ",", 1:5 ,"]"))

    }else if((!"iter"%in%names(dots)) & "control"%in%names(dots)){

      smpl <- as.matrix(sampling(
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
          doMAP = 1,
          saturating = ifelse(saturating, 1, 0),
          mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
                      prior.mu$mu_a2[1], prior.mu$mu_b2[1],prior.mu$mu_eta[1]),
          sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
                    prior.mu$mu_a2[2], prior.mu$mu_b2[2],prior.mu$mu_eta[2]),
          mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
                       prior.tau$tau_a2[1], prior.tau$tau_b2[1],prior.tau$tau_eta[1]),
          sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
                     prior.tau$tau_a2[2], prior.tau$tau_b2[2],prior.tau$tau_eta[2])
        ),
        iter = 10000,
        refresh = 0,
        ...),
        pars = paste0("log_ab[",nstd+1, ",", 1:5 ,"]"))


    }else if((!"control"%in%names(dots)) & "iter"%in%names(dots)){


      smpl <- as.matrix(sampling(
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
          doMAP = 1,
          saturating = ifelse(saturating, 1, 0),
          mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
                      prior.mu$mu_a2[1], prior.mu$mu_b2[1],prior.mu$mu_eta[1]),
          sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
                    prior.mu$mu_a2[2], prior.mu$mu_b2[2],prior.mu$mu_eta[2]),
          mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
                       prior.tau$tau_a2[1], prior.tau$tau_b2[1],prior.tau$tau_eta[1]),
          sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
                     prior.tau$tau_a2[2], prior.tau$tau_b2[2],prior.tau$tau_eta[2])
        ),
        control = list(adapt_delta=0.9),
        refresh = 0,
        ...),
        pars = paste0("log_ab[",nstd+1, ",", 1:5 ,"]"))


    }else {

      smpl <- as.matrix(sampling(
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
          doMAP = 1,
          saturating = ifelse(saturating, 1, 0),
          mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
                      prior.mu$mu_a2[1], prior.mu$mu_b2[1],prior.mu$mu_eta[1]),
          sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
                    prior.mu$mu_a2[2], prior.mu$mu_b2[2],prior.mu$mu_eta[2]),
          mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
                       prior.tau$tau_a2[1], prior.tau$tau_b2[1],prior.tau$tau_eta[1]),
          sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
                     prior.tau$tau_a2[2], prior.tau$tau_b2[2],prior.tau$tau_eta[2])
        ),
        iter = 10000,
        control = list(adapt_delta=0.9),
        refresh = 0,
        ...),
        pars = paste0("log_ab[",nstd+1, ",", 1:5 ,"]"))

    }

  }
  #extract samples per parameter
  la1 <- smpl[, paste0("log_ab[",nstd+1, ",", 1 ,"]")]
  lb1 <- smpl[, paste0("log_ab[",nstd+1, ",", 2 ,"]")]
  la2 <- smpl[, paste0("log_ab[",nstd+1, ",", 3 ,"]")]
  lb2 <- smpl[, paste0("log_ab[",nstd+1, ",", 4 ,"]")]
  eta <- smpl[, paste0("log_ab[",nstd+1, ",", 5 ,"]")]


  n.quan <- length(probs)
  summ <- matrix(nrow = 5, ncol = 2 + n.quan)
  rownames(summ) <- c("log(alpha_1)", "log(beta_1)","log(alpha_2)", "log(beta_2)", "eta")
  colnames(summ) <- c("mean", "SD", paste0(round(probs*100, digits=4), "%"))
  summ[1, 1] <- mean(la1)
  summ[2, 1] <- mean(lb1)
  summ[3, 1] <- mean(la2)
  summ[4, 1] <- mean(lb2)
  summ[5, 1] <- mean(eta)
  summ[1, 2] <- sd(la1)
  summ[2, 2] <- sd(lb1)
  summ[3, 2] <- sd(la2)
  summ[4, 2] <- sd(lb2)
  summ[5, 2] <- sd(eta)
  summ[1, 3:(n.quan+2)] <- quantile(la1, probs=probs)
  summ[2, 3:(n.quan+2)] <- quantile(lb1, probs=probs)
  summ[3, 3:(n.quan+2)] <- quantile(la2, probs=probs)
  summ[4, 3:(n.quan+2)] <- quantile(lb2, probs=probs)
  summ[5, 3:(n.quan+2)] <- quantile(eta, probs=probs)

  if(!return.samples){
    return(summ)
  }else{
    smpl_ret <- t(smpl)
    rownames(smpl_ret) <- c("log(alpha_1)", "log(beta_1)","log(alpha_2)", "log(beta_2)", "eta")
    return(list(
      summary = summ,
      samples = smpl_ret
    ))
  }
}






#-------------------------------------------------------------------------------
#Function: get_MAP
#-------------------------------------------------------------------------------

#' @rdname get_MAP
#' @param dose Numerical vector that lists the dose levels of a monotherapy compound applied to each
#' patient. The length must be equal to the number of cohorts in the data.
#' @param dose.ref Positive number, reference dose for the monotherapy compound.
#' @param prior.mu.mono Optional named list that specifies the prior. Must have the
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
#' @examples
#' \dontrun{
#' MAPprior <- get_MAP_mono(
#'    dose =  c(1, 2, 3, 4),
#'    n.pat = c(3, 3, 3, 6),
#'    n.dlt = c(0, 0, 0, 1),
#'    trial = c(1, 2, 1, 2),
#'    dose.ref = 5,
#'    return.samples = TRUE
#' )
#'
#' MAPprior$summary
#' MAPprior$samples[, 1:8]
#' }
#' @usage
#' get_MAP_mono(
#'    dose,
#'    dose.ref,
#'    n.pat,
#'    n.dlt,
#'    trial,
#'    prior.mu.mono = list(mu_a =  c(logit(0.33), 2),
#'                         mu_b =  c(0,          1)),
#'    prior.tau.mono = list(tau_a =  c(log(0.25),  log(2)/1.96),
#'                          tau_b =  c(log(0.125), log(2)/1.96)),
#'    return.samples = FALSE,
#'    probs = c(0.025, 0.5, 0.975),
#'    ...
#' )
#' @export
get_MAP_mono <- function(
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
  probs = c(0.025,  0.5, 0.975),
  ...
){
  nobs <- length(dose)
  if(!nobs>1){
    stop("dose must have at least one entry.")
  }
  if(! (nobs==length(n.pat)
        & nobs==length(n.dlt)
        & nobs==length(trial)
  )
  ){
    stop("`dose`, `n.pat`, `n.dlt`, and `trial` must have the same length.")
  }

  if(!is.logical(return.samples)){
    stop("`return.samples` must be logical.")
  }
  if(!is.num(dose, low=0, up= Inf, len=nobs, lB=F, uB=F)){
    stop("`dose` must be a numerical with positive entries.")
  }
  if(!is.num(dose.ref, low=0, up= Inf, len=1, lB=F, uB=F)){
    stop("`dose.ref` must be a positive single number.")
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
  if(!is.list(prior.mu.mono)){
    stop("`prior.mu.mono` must be a named list.")
  }
  if(is.null(names(prior.mu.mono))){
    stop("`prior.mu.mono` must be a named list.")
  }
  names(prior.mu.mono) <- tolower(names(prior.mu.mono))
  test.prior.mu.mono(prior.mu.mono)

  if(!is.list(prior.tau.mono)){
    stop("`prior.tau.mono` must be a named list.")
  }
  if(is.null(names(prior.tau.mono))){
    stop("`prior.tau.mono` must be a named list.")
  }
  names(prior.tau.mono) <- tolower(names(prior.tau.mono))
  test.prior.tau.mono(prior.tau.mono)

  if(!is.num(probs, low=0, up=1)){
    stop("`probs` must be a numeric vector with entries between 0 and 1.")
  }
  if(!length(probs)>0){
    stop("`probs` must be at least of length 1.")
  }

  dots <- list(...)
  if("refresh"%in%names(dots)){
    if("iter"%in%names(dots) & "control"%in%names(dots)){
      smpl <- as.matrix(sampling(
        object = stanmodels$jointBLRM,
        data = list(
          n_obs = nobs,
          n_studies = nstd,
          dose_1 = as.array(dose),
          dose_2 = as.array(rep(0, times=nobs)),
          dose_c = c(dose.ref, 1),
          n = as.array(n.pat),
          r = as.array(n.dlt),
          s = as.array(trial),
          doMAP = 1,
          saturating = 0,
          mean_mu = c(prior.mu.mono$mu_a[1], prior.mu.mono$mu_b[1],
                      -0.8, 0, 0),
          sd_mu = c(prior.mu.mono$mu_a[2], prior.mu.mono$mu_b[2],
                    2, 1, 1),
          mean_tau = c(prior.tau.mono$tau_a[1], prior.tau.mono$tau_b[1],
                       log(0.25), log(0.125), log(0.125)),
          sd_tau = c(prior.tau.mono$tau_a[2], prior.tau.mono$tau_b[2],
                     log(2)/1.96, log(2)/1.96, log(2)/1.96)
        ),
        ...),
        pars = paste0("log_ab[",nstd+1, ",", 1:2 ,"]"))
    }else if((!"iter"%in%names(dots)) & "control"%in%names(dots)){
      smpl <- as.matrix(sampling(
        object = stanmodels$jointBLRM,
        data = list(
          n_obs = nobs,
          n_studies = nstd,
          dose_1 = as.array(dose),
          dose_2 = as.array(rep(0, times=nobs)),
          dose_c = c(dose.ref, 1),
          n = as.array(n.pat),
          r = as.array(n.dlt),
          s = as.array(trial),
          doMAP = 1,
          saturating = 0,
          mean_mu = c(prior.mu.mono$mu_a[1], prior.mu.mono$mu_b[1],
                      -0.8, 0, 0),
          sd_mu = c(prior.mu.mono$mu_a[2], prior.mu.mono$mu_b[2],
                    2, 1, 1),
          mean_tau = c(prior.tau.mono$tau_a[1], prior.tau.mono$tau_b[1],
                       log(0.25), log(0.125), log(0.125)),
          sd_tau = c(prior.tau.mono$tau_a[2], prior.tau.mono$tau_b[2],
                     log(2)/1.96, log(2)/1.96, log(2)/1.96)
        ),
        iter = 10000,
        ...),
        pars = paste0("log_ab[",nstd+1, ",", 1:2 ,"]"))
    }else if((!"control"%in%names(dots)) & "iter"%in%names(dots)){
      smpl <- as.matrix(sampling(
        object = stanmodels$jointBLRM,
        data = list(
          n_obs = nobs,
          n_studies = nstd,
          dose_1 = as.array(dose),
          dose_2 = as.array(rep(0, times=nobs)),
          dose_c = c(dose.ref, 1),
          n = as.array(n.pat),
          r = as.array(n.dlt),
          s = as.array(trial),
          doMAP = 1,
          saturating = 0,
          mean_mu = c(prior.mu.mono$mu_a[1], prior.mu.mono$mu_b[1],
                      -0.8, 0, 0),
          sd_mu = c(prior.mu.mono$mu_a[2], prior.mu.mono$mu_b[2],
                    2, 1, 1),
          mean_tau = c(prior.tau.mono$tau_a[1], prior.tau.mono$tau_b[1],
                       log(0.25), log(0.125), log(0.125)),
          sd_tau = c(prior.tau.mono$tau_a[2], prior.tau.mono$tau_b[2],
                     log(2)/1.96, log(2)/1.96, log(2)/1.96)
        ),
        control = list(adapt_delta=0.9),
        ...),
        pars = paste0("log_ab[",nstd+1, ",", 1:2 ,"]"))
    }else {
      smpl <- as.matrix(sampling(
        object = stanmodels$jointBLRM,
        data = list(
          n_obs = nobs,
          n_studies = nstd,
          dose_1 = as.array(dose),
          dose_2 = as.array(rep(0, times=nobs)),
          dose_c = c(dose.ref, 1),
          n = as.array(n.pat),
          r = as.array(n.dlt),
          s = as.array(trial),
          doMAP = 1,
          saturating = 0,
          mean_mu = c(prior.mu.mono$mu_a[1], prior.mu.mono$mu_b[1],
                      -0.8, 0, 0),
          sd_mu = c(prior.mu.mono$mu_a[2], prior.mu.mono$mu_b[2],
                    2, 1, 1),
          mean_tau = c(prior.tau.mono$tau_a[1], prior.tau.mono$tau_b[1],
                       log(0.25), log(0.125), log(0.125)),
          sd_tau = c(prior.tau.mono$tau_a[2], prior.tau.mono$tau_b[2],
                     log(2)/1.96, log(2)/1.96, log(2)/1.96)
        ),
        control = list(adapt_delta=0.9),
        iter = 10000,
        ...),
        pars = paste0("log_ab[",nstd+1, ",", 1:2 ,"]"))
    }
  }else{
    #set refresh to 0 manually
    if("iter"%in%names(dots) & "control"%in%names(dots)){
      smpl <- as.matrix(sampling(
        object = stanmodels$jointBLRM,
        data = list(
          n_obs = nobs,
          n_studies = nstd,
          dose_1 = as.array(dose),
          dose_2 = as.array(rep(0, times=nobs)),
          dose_c = c(dose.ref, 1),
          n = as.array(n.pat),
          r = as.array(n.dlt),
          s = as.array(trial),
          doMAP = 1,
          saturating = 0,
          mean_mu = c(prior.mu.mono$mu_a[1], prior.mu.mono$mu_b[1],
                      -0.8, 0, 0),
          sd_mu = c(prior.mu.mono$mu_a[2], prior.mu.mono$mu_b[2],
                    2, 1, 1),
          mean_tau = c(prior.tau.mono$tau_a[1], prior.tau.mono$tau_b[1],
                       log(0.25), log(0.125), log(0.125)),
          sd_tau = c(prior.tau.mono$tau_a[2], prior.tau.mono$tau_b[2],
                     log(2)/1.96, log(2)/1.96, log(2)/1.96)
        ),
        refresh = 0,
        ...),
        pars = paste0("log_ab[",nstd+1, ",", 1:2 ,"]"))
    }else if((!"iter"%in%names(dots)) & "control"%in%names(dots)){
      smpl <- as.matrix(sampling(
        object = stanmodels$jointBLRM,
        data = list(
          n_obs = nobs,
          n_studies = nstd,
          dose_1 = as.array(dose),
          dose_2 = as.array(rep(0, times=nobs)),
          dose_c = c(dose.ref, 1),
          n = as.array(n.pat),
          r = as.array(n.dlt),
          s = as.array(trial),
          doMAP = 1,
          saturating = 0,
          mean_mu = c(prior.mu.mono$mu_a[1], prior.mu.mono$mu_b[1],
                      -0.8, 0, 0),
          sd_mu = c(prior.mu.mono$mu_a[2], prior.mu.mono$mu_b[2],
                    2, 1, 1),
          mean_tau = c(prior.tau.mono$tau_a[1], prior.tau.mono$tau_b[1],
                       log(0.25), log(0.125), log(0.125)),
          sd_tau = c(prior.tau.mono$tau_a[2], prior.tau.mono$tau_b[2],
                     log(2)/1.96, log(2)/1.96, log(2)/1.96)
        ),
        iter = 10000,
        refresh = 0,
        ...),
        pars = paste0("log_ab[",nstd+1, ",", 1:2 ,"]"))
    }else if((!"control"%in%names(dots)) & "iter"%in%names(dots)){
      smpl <- as.matrix(sampling(
        object = stanmodels$jointBLRM,
        data = list(
          n_obs = nobs,
          n_studies = nstd,
          dose_1 = as.array(dose),
          dose_2 = as.array(rep(0, times=nobs)),
          dose_c = c(dose.ref, 1),
          n = as.array(n.pat),
          r = as.array(n.dlt),
          s = as.array(trial),
          doMAP = 1,
          saturating = 0,
          mean_mu = c(prior.mu.mono$mu_a[1], prior.mu.mono$mu_b[1],
                      -0.8, 0, 0),
          sd_mu = c(prior.mu.mono$mu_a[2], prior.mu.mono$mu_b[2],
                    2, 1, 1),
          mean_tau = c(prior.tau.mono$tau_a[1], prior.tau.mono$tau_b[1],
                       log(0.25), log(0.125), log(0.125)),
          sd_tau = c(prior.tau.mono$tau_a[2], prior.tau.mono$tau_b[2],
                     log(2)/1.96, log(2)/1.96, log(2)/1.96)
        ),
        control = list(adapt_delta=0.9),
        refresh = 0,
        ...),
        pars = paste0("log_ab[",nstd+1, ",", 1:2 ,"]"))
    }else {
      smpl <- as.matrix(sampling(
        object = stanmodels$jointBLRM,
        data = list(
          n_obs = nobs,
          n_studies = nstd,
          dose_1 = as.array(dose),
          dose_2 = as.array(rep(0, times=nobs)),
          dose_c = c(dose.ref, 1),
          n = as.array(n.pat),
          r = as.array(n.dlt),
          s = as.array(trial),
          doMAP = 1,
          saturating = 0,
          mean_mu = c(prior.mu.mono$mu_a[1], prior.mu.mono$mu_b[1],
                      -0.8, 0, 0),
          sd_mu = c(prior.mu.mono$mu_a[2], prior.mu.mono$mu_b[2],
                    2, 1, 1),
          mean_tau = c(prior.tau.mono$tau_a[1], prior.tau.mono$tau_b[1],
                       log(0.25), log(0.125), log(0.125)),
          sd_tau = c(prior.tau.mono$tau_a[2], prior.tau.mono$tau_b[2],
                     log(2)/1.96, log(2)/1.96, log(2)/1.96)
        ),
        control = list(adapt_delta=0.9),
        iter = 10000,
        refresh = 0,
        ...),
        pars = paste0("log_ab[",nstd+1, ",", 1:2 ,"]"))
    }
  }

  #extract samples per parameter
  log_alpha <- smpl[, paste0("log_ab[",nstd+1, ",", 1 ,"]")]
  log_beta <- smpl[, paste0("log_ab[",nstd+1, ",", 2 ,"]")]

  n.quan <- length(probs)
  summ <- matrix(nrow = 2, ncol = 2 + n.quan)
  rownames(summ) <- c("log(alpha)", "log(beta)")
  colnames(summ) <- c("mean", "SD", paste0(round(probs*100, digits=4), "%"))
  summ[1, 1] <- mean(log_alpha)
  summ[2, 1] <- mean(log_beta)
  summ[1, 2] <- sd(log_alpha)
  summ[2, 2] <- sd(log_beta)
  summ[1, 3:(n.quan+2)] <- quantile(log_alpha, probs=probs)
  summ[2, 3:(n.quan+2)] <- quantile(log_beta, probs=probs)

  if(!return.samples){
    return(summ)
  }else{
    smpl_ret <- t(smpl)
    rownames(smpl_ret) <- c("log(alpha)", "log(beta)")
    return(list(
      summary = summ,
      samples = smpl_ret
    ))
  }
}


