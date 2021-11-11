#'Binomial-beta-mixture model for single dose levels
#'@name binomial_betamix
#'@rdname binomial_betamix
#'@aliases sim_binomial_betamix
#'@description Evaluates a Bayesian model with binomial likelihood and beta-mixture as prior
#'for binary DLT data at a single dose level, where the DLT rate of this dose
#' is the parameter of interest. The mixture prior can be given manually
#'or, alternatively, an (e.g.) MCMC sample from a predictive prior (e.g. MAP prior) can
#'be used to derive the beta mixture that the model employs as prior. Function
#'\code{binomial_betamix()} derives the beta mixture (if a sample is given),
#'adds optional additional prior components for robustification, and determines
#'the resulting posterior based on given data.  Function \code{sim_binomial_betamix()} allows to
#'estimate the probability of confirmation if the binomial-beta-mixture model
#'is used to perform a dose confirmation (in terms of the safety of the DLT rate)
#'for a single dose level, where confirmation is defined via the EWOC criterion from
#'Babb et al. (1998) after a fixed number of patients. The probability of confirmation
#'is estimated under the assumption
#'of an underlying true DLT rate of the dose of interest and a fixed number of new patients
#'treated at it. The model and function are defined for such dose confirmations of single
#'dose levels, which could e.g. be performed to
#'confirm the safety of a dose of interest in a new population or e.g. the largest dose of interest, if
#'multiple dose levels shall be tested in the new population. Refer to the corresponding
#'vignette on single-dose confirmation for further detail.
#'@param n.pat Positive whole number, gives the number of patients in the dose confirmation cohort.
#'@param n.dlt Positive whole number, gives the number of patients with DLTs in the dose confirmation cohort.
#'@param sample.pred.tox Optional numerical vector (or matrix with either one row or one column)
#'that contains samples from the predictive distribution of the DLT rate from
#'the dose level of interest. The sample is approximated with a mixture of beta
#'distributions via the function \code{\link[RBesT:mixfit]{RBesT::mixfit}()}.
#'A sample from the predictive distribution can e.g. be obtained via
#'\code{\link[OncoBLRM:get_pred_tox]{get_pred_tox}()}. Note that either \code{samples.pred.tox}
#'of \code{robust.comp} must be specified.
#'@param Nc Optional positive whole number, defaults to 3. Is passed to \code{\link[RBesT:mixfit]{RBesT::mixfit}()}.
#'Indicates the number of mixture components to be used during approximation with a
#'beta mixture.
#'@param robust.comp Optional numerical vector of length 3 with positive entries,
#'or a list with multiple such vectors. Indicates the robustifying mixture components
#'to be used during the analysis. When the value is a vector of length 3, it provides
#'the weight \eqn{w} (between 0 and 1) and shape parameters \eqn{a>0} and \eqn{b>0} for a beta distribution
#'\eqn{\mathrm{Beta}(a, b)} with mixture weight \eqn{w}, i.e., is of the form
#'\code{robust.comp = c(w, a, b)}. Otherwise, \code{robust.comp} must be a list with
#'one or more such vectors of lenght 3 as entries, where each list entry is interpreted
#'in the same way as before, i.e., as \code{c(w, a, b)} with weight \code{w} and
#'shape parameters \code{a} and \code{b}. Note that the sum of the given weights
#'must not exceed 1, but can be equal to 1 in case \code{sample.pred.tox} is \code{NULL}.
#'Note that the mixture weights of the mixture that approximates \code{sample.pred.tox}
#'will be rescaled accordingly to ensure that the robustifying compounds have the
#'mixture weights given in \code{robust.comp}.
#'
#'If \code{NULL} is given as value of \code{robust.comp} or if the parameter is left out, the function
#'will not perform robustification and only use the given \code{sample.pred.tox} to determine
#'a beta mixture as prior. Note that either \code{samples.pred.tox}
#'of \code{robust.comp} must be specified.
#'@param probs Optional numerical vector with values between 0 and 1. Provides the
#'quantiles to be included in the summary. Defaults to \code{c(0.025, 0.5, 0.975)}.
#'@param dosing.intervals Numerical vector with two ascending entries between 0 and 1.
#'Indicates the boundaries of the target dosing interval.
#'Defaults to \code{c(0.16, 0.33)}. Note:  function \code{sim_binomial_betamix()}
#'implicitly assumes that the overdosing interval consists of the DLT rates greater than
#'\code{dosing.intervals[2]}.
#'@param return.prior Optional logical, indicates whether the prior mixture
#'shall be contained in the output. Defaults to \code{FALSE}. The prior is
#'returned as a list of its components, where each component is described
#'by a vector that contains the weight and two shape parameters of the beta distribution.
#'This can be used to perform further evaluations on the prior, e.g. by using the
#'functions provided by the \code{\link[RBesT:RBesT-package]{RBesT-package}}.
#'@param constrain_gt1 Logical, given to \code{\link[RBesT:mixfit]{RBesT::mixfit}()}.
#'Default is \code{FALSE}. If the value is \code{TRUE}, \code{RBesT::mixfit()} will
#'constrain the parameters of the approximating Beta mixture to be greater or equal to 1.
#'Note that this is mostly required to ensure that the effective sample size as assessed
#'by the ELIR method (Neuenschwander et al., 2020) is well-defined, which allows
#'to calculate the effective sample size using \code{\link[RBesT:ess]{RBesT::ess}()}.
#'@param ... Optional further parameters given to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()}.
#'@returns For function \code{binomial_betamix()}: List with the entries
#'\itemize{
#'\item{\code{...$summary}:\cr
#'Matrix that contains mean, SD, quantiles and interval probabilities of the posterior DLT rate.}
#'\item{\code{...$posterior}:\cr
#'List that describes the posterior distribution analytically as a mixture of beta distributions.
#'This is possible as the binomial-beta-mixture model is conjugated. Each list entry describes
#'a compound of the posterior distribution by giving the weight and the two shape parameters
#'of the beta distribution.
#'This can be used to perform further evaluations on the posterior, e.g. by using the
#'functions provided by the \code{\link[RBesT:RBesT-package]{RBesT-package}}.}
#'}
#'In case \code{return.prior} is set to \code{TRUE}, the following additional list
#'entry will be contained in the output:
#'\itemize{
#'\item{\code{...$prior}:\cr
#'List that describes the prior distribution as a mixture of beta distributions.
#'Each list entry describes
#'a compound of the prior distribution by giving the weight and the two shape parameters
#'of the beta distribution.
#'}}
#'
#'For function \code{sim_binomial_betamix()}: Single numerical value, which is the simulated probability of
#'confirming the dose based on the assumed DLT rate and number of patients.
#'@usage
#'binomial_betamix(
#'  n.pat,
#'  n.dlt,
#'  sample.pred.tox=NULL,
#'  Nc = 3,
#'  robust.comp = NULL,
#'  probs = c(0.025, 0.5, 0.975),
#'  dosing.intervals = c(0.16, 0.33),
#'  return.prior = FALSE,
#'  constrain_gt1 = TRUE,
#'  ...
#')
#'@details
#'The binomial-beta-mixture model assumes the likelihood
#'\deqn{
#' r \vert \pi_d \ \sim \ \mathrm{Binomial}(n, \pi_d)
#'}
#'where \eqn{n} denotes the number of patients treated at the dose \eqn{d}, while
#'\eqn{r} denotes the observed number of patients with DLT (which is assumed to be
#'a realization of an underlying random variable), and \eqn{\pi_d} denotes the
#'model parameter that encodes the DLT rate of the dose \eqn{d}.
#'
#'The model is completed by assuming a prior distribution defined by a beta mixture, i.e.
#'is given by a probability density function of the form
#'\deqn{
#'p(\pi_d) = \sum\limits_{i=1}^{N_c} w_i \cdot p_i(\pi_d),
#'}
#'where \eqn{N_c} is
#'the number of mixture components, \eqn{0<w_i\leq 1} are mixture weights
#'(prior probabilities of the mixture components) with \eqn{sum_i w_i = 1},
#'while \eqn{p_i(\cdot)} is the density function of the \eqn{i^{th}} mixture component,
#'and \eqn{p(\cdot)} the density of the prior.
#'TO obtain a beta mixture, we assume that
#'\deqn{
#'p_i(x) \propto x^{a-1}\cdot (1-x)^{b-1}
#'}
#'is the density of a beta distribution \eqn{\mathrm{Beta}(a, b)} with shape parameters
#'\eqn{a>0} and \eqn{b>0}. Using standard calculations, one verifies readily that
#'the resulting model is conjugated and that the posterior \eqn{p(\pi_d \vert r)}
#'can again be written as a mixture of beta distributions. The function will compute both the
#'analytical expression of the posterior as a beta mixture, but can also draw MCMC
#'samples from the posterior.
#'
#'The function \code{binomial_betamix()} does not only compute the posterior of
#'this model, but also allows to determine a beta mixture based on an MCMC sample
#'from a predictive distribution of the DLT rate. That is, one can e.g. use
#'meta-analytic predictive distributions of the DLT rate as obtained from a joint BLRM
#'and historical data to construct a beta mixture which is then used as the prior
#'in the binomial-beta-mixture model. Additionally, one may provide the weight
#'and shape parameters of additional mixture components which are added to the
#'approximation for robustification against prior-data conflicts (see e.g. Schmidli et al., 2014,
#'or Neuenschwander et al., 2020).
#'
#'The approximation of the MCMC sample with a mixture of beta distributions (with
#'a prespecified number of components) is carried out via the function
#' \code{\link[RBesT:mixfit]{RBesT::mixfit}()} from the \code{\link[RBesT:RBesT-package]{RBesT-package}}.
#' The function uses the so-called Expectation Maximization algorithm for this,
#' refer to the documentation of it for further detail.
#'
#'The function \code{binomial_betamix()} is designed to be used in combination with \code{\link[OncoBLRM:get_pred_tox]{get_pred_tox}()},
#'which allows to derive predictive distributions of the DLT rate for individual dose levels. The function
#'constructs an MCMC sample from this distribution, which can then be approximated with a mixture
#'of beta distributions by the function \code{binomial_betamix()}. A detailed description
#'of how these functions are intended to be combined can be found in a vignette of the
#'package \code{\link[OncoBLRM:OncoBLRM-package]{OncoBLRM-package}}.
#'
#'Function \code{sim_binomial_betamix()} assumes that \code{n.pat} patients are
#'treated at a new dose with a fixed DLT rate \code{p.dlt}, with the aim of
#'confirming the safety of this dose. Confirmation is in this setting defined by
#'an overdosing probability below \code{confirmation.bound}, where the overdosing
#'probability is the posterior probability
#'\deqn{
#'\Pr(\pi_d\geq \delta)
#'}
#' based on the binomia-beta-mixture model and for a fixed boundary \eqn{\delta}
#' which is given in \code{dosing.intervals[2]} (upper boundary of the target dosing interval).
#'
#'Using this, it can be estimated how likely confirmation via simulation. This
#'simulation is performed by function \code{sim_binomial_betamix()}.
#'
#' @seealso \code{\link[OncoBLRM:get_pred_tox]{get_pred_tox}()}, \code{\link[OncoBLRM:get_MAP]{get_MAP}()},
#' \code{\link[OncoBLRM:fit_jointBLRM]{fit_jointBLRM}()},
#' \code{\link[RBesT:mixfit]{RBesT::mixfit}()}, \code{\link[RBesT:RBesT-package]{RBesT-package}}.
#'
#' #' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. \url{https://mc-stan.org}.
#'
#' Schmidli, H., Gsteiger, S., Roychoudhury, S., O'Hagan, A., Spiegelhalter, D., & Neuenschwander B. (2014).
#' Robust meta-analytic-predictive priors in clinical trials with historical control information.
#'
#' Neuenschwander, B., Roychoudhury, S., & Schmidli, H. (2016). On the use of co-data in clinical trials.
#' Statistics in Biopharmaceutical Research, 8(3), 345-354, <doi:10.1080/19466315.2016.1174149>.
#'
#' Babb, J., Rogatko, A., & Zacks, S. (1998). Cancer phase I clinical trials: Efficient dose escalation with overdose control.
#' Statistics in medicine 17(10), 1103-1120.
#'
#'
#'
#'@export
binomial_betamix <- function(
  n.pat,
  n.dlt,
  sample.pred.tox=NULL,
  Nc = 3,
  robust.comp = NULL,
  probs = c(0.025, 0.5, 0.975),
  dosing.intervals = c(0.16, 0.33),
  return.prior = FALSE,
  constrain_gt1 = TRUE,
  ...
){
  if(!is.logical(return.prior)){
    stop("`return.prior` must be of type logical.")
  }
  if(missing(n.pat)){
    n.pat <- 0
  }
  if(missing(n.dlt)){
    if(n.pat==0){
      n.dlt <-0
    }else{
      stop("When `n.pat` is positive, the argument `n.dlt` must be given.")
    }
  }
  if(!is.wholenumber(n.pat)){
    stop("`n.pat` must be a whole number.")
  }
  if(!n.pat>=0){
    stop("`n.pat` must be positive.")
  }
  if(!is.num(probs, low=0, up=1)){
    stop("`probs` must be a numeric vector with entries between 0 and 1.")
  }
  if(!length(probs)>0){
    stop("`probs` must be at least of length 1.")
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
  if(!is.wholenumber(n.dlt)){
    stop("`n.dlt` must be a whole number.")
  }
  if(!n.dlt>=0){
    stop("`n.dlt` must be positive.")
  }
  if(!n.dlt <=n.pat){
    stop("`n.dlt` must be smaller or equal to n.pat")
  }
  if(is.null(sample.pred.tox) & is.null(robust.comp)){
    stop("Either `sample.pred.tox` or `robust.comp` or both must be given.")
  }

  if(!is.wholenumber(Nc)){
    stop("`Nc` must be a whole number")
  }



  if(!is.null(sample.pred.tox)){
    dims <- dim(sample.pred.tox)
    if(!length(dims)%in%c(0, 2)){
      stop("`sample.pred.tox` must be a vector, or a matrix with either 1 row or with 1 column")
    }else if(length(dims)==2){
      if(!dims[1]==1 & !dims[2]==1){
        stop("`sample.pred.tox` must be a vector, or a matrix with either 1 row or with 1 column")
      }else if(dims[1]==1){
        smp <- sample.pred.tox[1, ]
      }else{
        smp <- sample.pred.tox[, 1]
      }
      #head(smp)
    }else{
      smp <- sample.pred.tox
    }
    nsmp <- length(smp)
  }

  if(!is.null(robust.comp)){
    if(!is.num(robust.comp, low = 0, lB=F, len=3) & !is.list(robust.comp)){
      stop("`robust.comp` must be numeric of length 3 with positive entries or a list of such vectors.")
    }
    if(is.num(robust.comp, low=0, lB=F, len=3)){
      Ncrob <- 1
      if(!is.num(robust.comp[1], low =0, up=1, lB=T, uB=T, len = 1)){
        stop("When `robust.comp` is not a list, `robust.comp[1]` is the mixture weight which must be a number between 0 and 1.")
      }
      robconfig <- matrix(NA, nrow = 3, ncol = 1)
      robconfig[1:3, 1] <-robust.comp
    }else{
      Ncrob <- length(robust.comp)
      robconfig <- matrix(NA, nrow = 3, ncol = Ncrob)
      for(comp_curr in 1:Ncrob){
        if(!is.num(robust.comp[[comp_curr]], low=0, lB=F, len=3)){
          stop("When `robust.comp` is a list, each entry must be a numeric with length 3 and positive entries.")
        }
        if(!robust.comp[[comp_curr]][1]<1){
          stop("When `robust.comp` is a list, entry `robust.comp[[i]][1]` must be a number between 0 and 1 which indicates the mixture weight.")
        }
        robconfig[1:3, comp_curr] <-robust.comp[[comp_curr]]
      }

    }
    if(!sum(robconfig[1, ])<=1){
      stop("The sum of the given mixture weights from `robust.comp` exceeds 1.")
    }else if(sum(robconfig[1, ])==1 & !is.null(sample.pred.tox)){
      stop("The sum of the given mixture weights from `robust.comp` is 1. This is not allowed when `sample.pred.tox` is given.")
    }


  }else{
    robust.weight <- 0
  }

  #derive beta mixture from sample and potentially robustify,
  #further, write mixture comps in matrix as input for stan
  if(!is.null(sample.pred.tox)){
    bmix1 <- mixfit(smp, type = "beta", Nc=Nc, constrain_gt1=constrain_gt1)
    if(!is.null(robust.comp)){
      nmix <- Nc+Ncrob
    }else{
      nmix <- Nc
    }
    bmix <- matrix(NA,  nrow = 3, ncol = nmix,
                   dimnames = list("parameter" = c("w", "a", "b"), "component" = 1:nmix))
    bmix[, 1:Nc] <- bmix1[, 1:Nc]
    if(!is.null(robust.comp)){
      bmix[1, 1:Nc] <- bmix[1, 1:Nc]*(1-sum(robconfig[1, ]))
      for(i in 1:Ncrob){
        bmix[1:3, Nc+i] <- robconfig[, i]
      }
      bmix[1, 1:nmix] <- bmix[1, 1:nmix]/sum(bmix[1, 1:nmix])
    }
  }else{
    #only fill in robust component
    nmix <- Ncrob
    bmix <- matrix(NA,  nrow = 3, ncol = nmix,
                   dimnames = list("parameter" = c("w", "a", "b"), "component" = 1:nmix))
    bmix[1:3, 1:nmix] <- robconfig[1:3, 1:nmix]
    bmix[1, 1:nmix] <-bmix[1, 1:nmix]/sum(bmix[1, 1:nmix])
  }

  dots <- list(...)
  if("iter" %in% names(dots) & "refresh" %in% names(dots)){
    smp_post <- as.matrix(sampling(
      object = stanmodels$betamix_binomial,
      data = list(
        nobs = 1,
        n = as.array(n.pat),
        r = as.array(n.dlt),
        num_mix_comps = nmix,
        mix_probs = as.array(bmix[1, ]),
        mix_a = as.array(bmix[2, ]),
        mix_b = as.array(bmix[3, ])
      ),
      pars = "p",
      ...
    ), pars = "p")
  }else if(!"iter" %in% names(dots) & "refresh" %in% names(dots)){
    smp_post <- as.matrix(sampling(
      object = stanmodels$betamix_binomial,
      data = list(
        nobs = 1,
        n = as.array(n.pat),
        r = as.array(n.dlt),
        num_mix_comps = nmix,
        mix_probs = as.array(bmix[1, ]),
        mix_a = as.array(bmix[2, ]),
        mix_b = as.array(bmix[3, ])
      ),
      pars = "p",
      iter = 5000,
      ...
    ), pars = "p")
  }else if("iter" %in% names(dots) & !"refresh" %in% names(dots)){
    smp_post <- as.matrix(sampling(
      object = stanmodels$betamix_binomial,
      data = list(
        nobs = 1,
        n = as.array(n.pat),
        r = as.array(n.dlt),
        num_mix_comps = nmix,
        mix_probs = as.array(bmix[1, ]),
        mix_a = as.array(bmix[2, ]),
        mix_b = as.array(bmix[3, ])
      ),
      pars = "p",
      refresh = 0,
      ...
    ), pars = "p")
  }else if(!"iter" %in% names(dots) & !"refresh" %in% names(dots)){
    smp_post <- as.matrix(sampling(
      object = stanmodels$betamix_binomial,
      data = list(
        nobs = 1,
        n = as.array(n.pat),
        r = as.array(n.dlt),
        num_mix_comps = nmix,
        mix_probs = as.array(bmix[1, ]),
        mix_a = as.array(bmix[2, ]),
        mix_b = as.array(bmix[3, ])
      ),
      pars = "p",
      iter = 5000,
      refresh = 0,
      ...
    ), pars = "p")
  }

  #prepare posterior summary from samples
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
  summ <- matrix(NA, nrow = 1, ncol = 2 + nquan + nint)
  colnames(summ) <- c("mean", "SD", paste0(round(100*probs, digits = 4), "%"), namesint)
  rownames(summ) <- "Dose"
  summ[1, 1] <- mean(smp_post[,1])
  summ[1, 2] <- sd(smp_post[,1])
  summ[1, 3:(2+nquan)] <- quantile(smp_post[,1], probs = probs)
  summ[1, (2+nquan+1):(2+nquan+nint)] <- get_int_probs_BLRM(smpl=smp_post[,1], ints = dosing.intervals)


  #also calculate posterior analytically
  post_param <- matrix(NA,  nrow = 3, ncol = nmix,
         dimnames = list("parameter" = c("w", "a", "b"), "component" = 1:nmix))
  #post_param <- matrix(NA, nrow = 3, ncol = nmix)
  #rownames(post_param) <- c("weight","a", "b")
  #colnames(post_param) <- paste0("comp. ", 1:nmix)
  c_i <- rep(0, times = nmix)

  for(k in 1:nmix){
    a_curr <- bmix[2, k]
    b_curr <- bmix[3, k]
    post_param[2, k] <-  a_curr+ n.dlt
    post_param[3, k] <-  b_curr+ n.pat - n.dlt
    c_i[k] <- choose(n.pat, n.dlt) * (gamma(a_curr + b_curr) / (gamma(a_curr)*gamma(b_curr) ))*
      ((gamma(a_curr + n.dlt)*gamma(b_curr + n.pat -n.dlt) )  /
         gamma(a_curr + b_curr + n.pat))
  }
  c_all <- sum(c_i*bmix[1,])
  for(k in 1:nmix){
    post_param[1, k] <-  (bmix[1,k] * c_i[k])/c_all
  }

  if(return.prior){
    return(list("summary" = summ,
                "posterior" = lapply(1:nmix, function(i) return(post_param[,i])),
                "prior" = lapply(1:nmix, function(i) return(bmix[,i]))))
  }else{
    return(list("summary" = summ,
                "posterior" = lapply(1:nmix, function(i) return(post_param[,i]))))
  }
}


#'@rdname binomial_betamix
#'@param p.dlt Non-negative numerical between 0 and 1, provides the assumed true
#'DLT rate of the dose of interest.
#'@param n.studies Positive integer, number of simulated studies.
#'@param confirmation.bound Maximum allowed probability of having a DLT rate within
#'the overdosing interval (referred to as overdosing probability) based on
#'the posterior. The dose is considered to be confirmed if the overdosing probability
#'is below \code{confirmation.bound}, otherwise the dose is considered to be
#'non-confirmed. Note that the last dosing interval as specified by
#'\code{dosing.intervals} is assumed to be the overdosing interval. By default,
#'this means that the overdosing interval is between 0.33 and 1.
#'@usage
#'sim_binomial_betamix(
#'  n.pat,
#'  p.dlt,
#'  n.studies = 1,
#'  confirmation.bound = 0.25,
#'  sample.pred.tox=NULL,
#'  Nc = 3,
#'  robust.comp = NULL,
#'  dosing.intervals = c(0.16, 0.33),
#'  ...
#')
#'@examples
#'\dontrun{
#'
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
#'#evaluate binomial-beta-mixture model using a mixture approximation of
#'#the predictive DLT rate and Beta(1, 1) with weight 0.1 as robustifying
#'#component. Assumes 2 out of 3 treated new patients suffered DLT at
#'#the dose 10+24.
#'post <- binomial_betamix(
#'  n.pat = 3,
#'  n.dlt = 2,
#'  sample.pred.tox = pred$samples,
#'  Nc = 3,
#'  robust.comp = c(0.1, 1, 1),
#'  return.prior = TRUE
#')
#'#summary
#'post$summary
#'
#'#posterior determined analytically as beta mixture
#'post$posterior
#'
#'#beta-mixture prior defined by the sample and robustification.
#'post$prior
#'
#'#Simulate probability of confirmation (i.e., EWOC criterion satisfied) when
#'#using the binomial-beta model as specified above under the assumption that
#'#6 patients are treated at 10+24 and that the dose has a true DLT rate of 0.5
#'sim_binomial_betamix(
#'  n.pat = 6,
#'  p.dlt = 0.5,
#'  n.studies = 1000,
#'  confirmation.bound = 0.25,
#'  sample.pred.tox = pred$samples,
#'  Nc = 3,
#'  robust.comp = c(0.1, 1, 1)
#')
#'
#'#equivalent: use the previously derived prior again
#'#(avoids that the beta-approximation needs to be computed again)
#'sim_binomial_betamix(
#'  n.pat = 6,
#'  p.dlt = 0.5,
#'  n.studies = 1000,
#'  confirmation.bound = 0.25,
#'  robust.comp = post$prior
#')
#'}
#'@export

sim_binomial_betamix <- function(
  n.pat,
  p.dlt,
  n.studies = 1,
  confirmation.bound = 0.25,
  sample.pred.tox=NULL,
  Nc = 3,
  robust.comp = NULL,
  dosing.intervals = c(0.16, 0.33),
  ...
) {
  if(!is.num(p.dlt, low=0, up=1, len=1)){
    stop("`p.dlt` must be a real number between 0 and 1")
  }
  if(!is.num(confirmation.bound, low=0, up=1, len=1)){
    stop("`confirmation.bound` must be a real number between 0 and 1")
  }

  if(!is.num(n.studies, low=0, lB=F, uB=F, len=1)){
    stop("`n.studies` must be an positive integer.")
  }
  if(!is.wholenumber(n.studies)){
    stop("`n.studies` must be an positive integer.")
  }


  n.dlt <- rbinom(n.studies, size = n.pat, prob = p.dlt)
  n.sims <- rep(0, times = n.pat + 1)
  res <- rep(0, times = n.pat + 1)
  for (i in 0:n.pat) {
    n.sims[i + 1] <- sum(ifelse(n.dlt == i, 1, 0))
    if (n.sims[i + 1] > 0) {
      curr <- binomial_betamix(
        n.pat = n.pat,
        n.dlt = i,
        sample.pred.tox = sample.pred.tox,
        Nc = Nc,
        robust.comp = robust.comp,
        dosing.intervals = dosing.intervals,
        ...
      )
      p_over <- curr$summary[1, length(curr$summary[1,])]
      if (p_over >= confirmation.bound) {
        res[i + 1] <- 0
      } else{
        res[i + 1] <- 1
      }
    }
  }
  n_conf <- sum(res*n.sims)

  return(n_conf/n.studies)
}
