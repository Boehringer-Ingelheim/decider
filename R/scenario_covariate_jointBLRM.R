#'@title Evaluate data from oncology dose finding using a joint BLRM with covariate
#'@description Evaluates data scenarios consisting of observations of one or more
#'monotherapy or two-drug combination therapy dose-finding trials including a
#'binary covariate and computes
#'posterior toxicities for a trial of interest and a set of doses of interest.
#'The function supports one-sided and two-sided covariates. This can be
#'controlled separately for each of thw two compounds.
#'
#'If multiple scenarios need to be evaluated, consider using the function
#'\code{\link[decider:scenario_list_covariate_jointBLRM]{scenario_list_covariate_jointBLRM}()}
#'instead, which is a parallelized wrapper that processes a list of data scenarios
#'within the same setting via \code{scenario_covariate_jointBLRM()}.
#'
#'A description of the underlying model and methods are given in the section Details.
#'@usage
#'scenario_covariate_jointBLRM(
#'    data=NULL,
#'    historical.data=NULL,
#'    doses.of.interest,
#'    dose.ref1,
#'    dose.ref2,
#'    trials.of.interest,
#'    types.of.interest=NULL,
#'    trials.of.interest.covars=NULL,
#'    esc.rule=c("ewoc", "loss", "dynamic.loss"),
#'    dosing.intervals = c(0.16, 0.33, 0.6),
#'    ewoc.threshold = 0.25,
#'    loss.weights = c(1, 0, 1, 2),
#'    dynamic.weights = rbind(c(0.32, 0, 0.32, 0.36),
#'                            c(0.29, 0, 0.31, 0.4),
#'                            c(0.27, 0, 0.33, 0.4),
#'                            c(0.2,  0, 0.3,  0.5)
#'    ),
#'    prior.mu = list(mu_a1 =  c(logit(0.33), 2),
#'                    mu_b1 =  c(0,           1),
#'                    mu_a2 =  c(logit(0.33), 2),
#'                    mu_b2 =  c(0,           1),
#'                    mu_eta = c(0,           1.121)
#'    ),
#'    prior.mu.covar = list(mu_g1 = c(0,  1),
#'                          mu_g2 = c(0,  1)
#'    ),
#'    prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
#'                     tau_b1 =  c(log(0.125), log(2)/1.96),
#'                     tau_a2 =  c(log(0.25),  log(2)/1.96),
#'                     tau_b2 =  c(log(0.125), log(2)/1.96),
#'                     tau_eta = c(log(0.125), log(2)/1.96)
#'    ),
#'    prior.tau.covar = list(tau_g1 = c(log(0.125), log(2)/1.96),
#'                           tau_g2 = c(log(0.125), log(2)/1.96)
#'    ),
#'    two_sided1 = TRUE,
#'    two_sided2 = TRUE,
#'    saturating = FALSE,
#'    probs = c(0.025, 0.5, 0.975),
#'    iter = 26000,
#'    warmup = 1000,
#'    refresh = 0,
#'    adapt_delta = 0.8,
#'    max_treedepth = 15,
#'    chains = 4,
#'    seed=sample.int(.Machine$integer.max, 1),
#'    path = NULL,
#'    file.name = NULL,
#'    plot.decisions = FALSE,
#'    plot.combi.heatmap = TRUE,
#'    plot.int.probs.loss = FALSE,
#'    plot.return = FALSE,
#'    plot.file.format = "pdf",
#'    plot.width,
#'    plot.height,
#'    plot.unit,
#'    output.scen.config = FALSE
#')
#'@param data List that contains the data scenario to be evaluated. Can be \code{NULL} if the prior
#'shall be computed. The data list should have the following named entries, all of which need to
#'be vectors of the same length (length should be the number of cohorts in the data):
#'\itemize{
#'\item{\code{data$dose1}\cr Numeric vector, each entry must be non-negative.
#'Entry \eqn{i} should provide the dose level of compound 1 administered to cohort \eqn{i} in the data.
#'Use \code{0} or \code{NA} to state that compound 1 was not used during treatment of a cohort.
#'Note: For each cohort, either \code{data$dose1} or \code{data$dose2} must be positive.}
#'\item{\code{data$dose2}\cr Numeric vector, each entry must be non-negative.
#'Entry \eqn{i} should provide the dose level of compound 2 administered to cohort \eqn{i} in the data.
#'Use \code{0} or \code{NA} to state that compound 2 was not used during treatment of a cohort.
#'Note: For each cohort, either \code{data$dose1} or \code{data$dose2} must be positive.}
#'\item{\code{data$n.pat}\cr Numeric vector, each entry must be a non-negative integer.
#'Entry \eqn{i} should provide the number of patients in cohort \eqn{i} in the data.
#'\code{0} is interpreted as the cohort not having been treated yet. If \code{data} contains solely
#'cohorts with 0 patients, the function will sample from the prior.}
#'\item{\code{data$n.dlt}\cr Numeric vector, each entry must be a non-negative integer.
#'Entry \eqn{i} should provide the number of DLTs in cohort \eqn{i} in the data. In particular, the
#'values need to be smaller or equal to the patient number of the corresponding cohort.}
#'\item{\code{data$trial}\cr Numeric or character vector. The entries should be trial names, i.e. indicators for the trial
#'to which the cohort belongs. These can either be numbers or strings (both will be converted to numbers internally).
#'Note: As mixed vectors of numbers and strings will be converted to strings, entries such as \code{1} and \code{"1"}
#'will be interpreted as the same trial.}
#'\item{\code{data$covar}\cr Numeric vector with 0-1 entries indicating the value of binary covariate for the
#'cohorts included in \code{data}. That is, when entry
#'\eqn{i} is 0, the value of the binary covariate is assumed to be 0 for the patients in cohort \eqn{i}, and similarly
#'for entry 1.
#'}}
#'@param historical.data Optional named list, must have the same structure as \code{data}. It is equivalent to include observations
#'as \code{data} or \code{historical.data}. Trial names across data and historical data must be consistent,
#'in the sense that observations with the same entries in \code{data$trial}, respectively \code{historical.data$trial}
#'are interpreted to belong to the same trial.
#'@param trials.of.interest Optional vector of numerical or character trial numbers/names, for which the posterior is to be computed.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param doses.of.interest Numeric matrix with two rows and non-negative entries. Each column gives a dose combination
#'of interest.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param types.of.interest Optional character vector with one entry for each entry of \code{trials.of.interest} which specifies the
#'trial type of the corresponding trial of interest. Supported trial types are \code{"mono1"}, \code{"mono2"},  \code{"combi"},
#'and \code{"all"}.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param trials.of.interest.covars Optional numerical vector with either length 1 or the same length as \code{trials.of.interest}.
#'Indicates (if applicable) the value of the covariate that is of interest for one or all trials of interest. Entries can be either
#'\code{0}, \code{1}, or \code{NA},
#'where \code{1} represents a value of 1 for the binary covariate, \code{0} indicates a value of 0, and \code{NA}
#'indicates that both values are of interest for a given trial.
#'Can be used to set a fixed covariate value for a trial of interest, which will cause the function to only compute, return, or plot
#'the results using this value for some trial of interest. If the parameter \code{trials.of.interest.covars}
#'is not given, the function will include the results for both values of the covariate for each trial
#'(equivalent to setting \code{trials.of.interest.covars = NA}).
#'If \code{trials.of.interest.covars} has length 1, the function will use the specified option for all trials of interest as given in
#'\code{trials.of.interest}, otherwise
#'it is assumed that \code{trials.of.interest.covars[i]} provides the option for the trial specified in \code{trials.of.interest[i]}.
#'@param dose.ref1 Numeric, must be positive. Reference dose for compound 1.
#'@param dose.ref2 Numeric, must be positive. Reference dose for compound 2.
#'@param esc.rule  Optional character. Can be either \code{"ewoc"}, \code{"loss"}, \code{"dynamic"} or \code{"dynamic.loss"}, where the
#'latter two are treated synonymously.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param dosing.intervals Optional numeric with 1, 2, or 3 ascending positive entries. Must have three entries when the \code{esc.rule} is
#'set to \code{"loss"}, \code{"dynamic.loss"}, or \code{"dynamic"}, otherwise (i.e. when \code{esc.rule} is \code{"ewoc"}) lengths 1, 2, or
#'3 are permitted.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param ewoc.threshold Optional numeric between 0 and 1. Overdosing thresholds for EWOC plots. Defaults to 0.25.
#'@param loss.weights Optional numerical vector with four entries (which can be arbitrary numbers), the default is \code{c(1,0,1,2)}.
#'Specifies the weights used for loss escalation.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param dynamic.weights Optional numerical matrix with four rows and four columns, and arbitrary numbers as entries.  Specifies the interval weights/penalties
#'that are used for dynamic loss escalation.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param prior.mu Optional list that gives the prior distribution for the hyper means \eqn{\mu}.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param prior.tau Optional list that gives the prior distribution for the between-trial heterogeneities (hyper SD) \eqn{\tau}.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param prior.mu.covar Optional named list that gives the prior distribution for the hyper-means of the additional parameters
#'included in the joint BLRM to realize the binary covariate. Also refer to the Details section for the notation used in the following.
#'The argument\code{prior.mu.covar} must be a list with named entries \code{mu_g1} and \code{mu_g2}.
#'Both must have length 2, and provide the mean and SD of the hyper-means for the parameters corresponding to the covariate
#'in compound 1 and 2. More precisely:
#'\itemize{
#'\item{\code{prior.mu.covar$mu_g1}\cr Numeric with length two, defaults to \code{c(0, 1)}. Specifies mean and SD of the hypermean \eqn{\mu_6} of the parameter \eqn{\gamma_1} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'
#'\item{\code{prior.mu.covar$mu_g2}\cr Numeric with length two, defaults to \code{c(0, 1)}. Specifies mean and SD of the hypermean \eqn{\mu_7} of the parameter \eqn{\gamma_2} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'}
#'@param prior.tau.covar Optional named list that gives the prior distribution for the between-trial heterogeneities of the additional parameters
#'included in the joint BLRM to realize the binary covariate. Also refer to the Details section for the notation used in the following.
#'The argument\code{prior.tau.covar} must be a list with named entries \code{tau_g1} and \code{tau_g2}.
#'Both must have length 2, and provide the mean and SD of the hyper-means for the parameters corresponding to the covariate
#'in compound 1 and 2. More precisely:
#'\itemize{
#'\item{\code{prior.tau.covar$tau_g1}\cr Numeric with length two, defaults to \code{c(0, 1)}. Specifies mean and SD of the between-trial
#'heterogeneity \eqn{\tau_{\gamma_1}} of the parameter \eqn{\gamma_1} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'
#'\item{\code{prior.tau.covar$tau_g2}\cr Numeric with length two, defaults to \code{c(log(0.125), log(2)/1.96)}. Specifies mean and SD of the between-trial
#'heterogeneity \eqn{\tau_{\gamma_2}} of the parameter \eqn{\gamma_2} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'}
#'@param two_sided1 Optional logical, defaults to \code{TRUE}. Indicates whether the covariate is assumed to have a two-sided effect
#'on the DLT rate of compound 1. If \code{FALSE}, the function will assume a one-sided effect of the covariate on the DLT rate for compound 1.
#'See the section Details below for a formal description.
#'@param two_sided2 Optional logical. Optional logical, defaults to \code{TRUE}. Indicates whether the covariate is assumed to have a two-sided effect
#'on the DLT rate of compound 2. If \code{FALSE}, the function will assume a one-sided effect of the covariate on the DLT rate for compound 2.
#'See the section Details below for a formal description.
#'@param saturating Optional logical that activates the use of a saturating interaction term (instead of linear), defaults to \code{FALSE}.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param probs Optional numeric with arbitrary entries between 0 and 1. Provides the levels for the quantiles displayed in
#'the output. Defaults to \code{c(0.025, 0.5, 0.975)}.
#'@param path Optional character that specified the path to save the resulting output files. If \code{NULL} (the default), no
#'output is written (but still returned to R). Otherwise, it is checked whether \code{path} specifies a directory, and, if yes,
#'all output is saved there.
#'@param file.name Optional name for the output file. If \code{NULL} or missing, no output is saved (and also if no valid path is given). Results are only returned to R in this case. Note
#'that plots will not be returned to R unless the additional argument \code{plot.return} is specified.
#'@param plot.decisions Optional logical, defaults to \code{FALSE}. If \code{TRUE}, plots of escalation decisions according to the
#'specified escalation rule are created.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param plot.combi.heatmap Optional logical, defaults to \code{TRUE}. If the value is \code{TRUE}, combination therapy plots
#'are created as heatmaps instead of bar plots. This affects all escalation rules.
#'@param plot.int.probs.loss Optional logical, defaults to \code{FALSE}. Only has an effect if \code{esc.rule} is either
#'\code{"loss"} or \code{"dynamic.loss"}. In this case, if the value is \code{TRUE}, additional plots will
#'be created that display the interval probabilities to complement the (always created) plots that display the resulting expected loss.
#'@param plot.return Optional logical, defaults to \code{FALSE}. If set to \code{TRUE}, the functions return the created plots
#'to R in the result list.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param plot.file.format Optional character, defaults to \code{"pdf"}. Can either be \code{"pdf"}, \code{"jpeg"} or \code{"png"}.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param plot.unit Optional character string, can be "in", "cm", or "mm".
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#' @param plot.width Optional numerical value or vector, can have length 1 or 3 and must have positive entries.
#' Provides the width of the output plots measured in the unit given in \code{plot.unit}.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#' @param plot.height Optional numerical value or vector, can have length 1 or 3 and must have positive entries.
#' Provides the height of the output plots measured in the unit given in \code{plot.unit}.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param iter Optional integer, number of total MCMC iterations per chain. Defaults to 26000. Note: Number of warmup iterations is counted towards \code{iter}, i.e.
#'of the \code{iter} many iterations, the first \code{warmup} many samples are not saved.
#'@param warmup Optional integer, number of warmup iterations discarded from total MCMC iterations per chain. Defaults to \code{1000}.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param chains Optional integer. Number of Markov chains constructed by Stan. Defaults to 4.
#'@param refresh Optional integer. Given to Stan's \code{refresh} argument for \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()}, defaults to \code{0}.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param adapt_delta Optional numeric between 0.6 and 1, default is 0.8. Given to Stan's \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()} method in the \code{control} argument of Stan.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param max_treedepth Optional integer, defaults to 15.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param seed Optional positive integer that specifies the seed to be used for the simulation.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@param output.scen.config Optional logical, defaults to \code{FALSE}.
#'See \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details.
#'@details
#'The basic joint BLRM is defined according to (Neuenschwander et al., 2014 and 2016). Refer to the section Details in the documentation of
#'\code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()} for details of the general model definition.
#'
#'In the following, we will describe how the function includes a binary covariate in the basic joint BLRM.
#'Recall the parameter vector used for the joint BLRM without covariates for trial \eqn{j}, i.e.
#'\deqn{\theta_j = (log(\alpha_{1j}), log(\beta_{1j}), log(\alpha_{2j}), log(\beta_{2j}), \eta_j).}
#'
#'To model the effect of a binary covariate \eqn{c\in\{0,1\}} on the DLT rate \eqn{\pi_j(d_i)} of dose \eqn{d_i} of
#'compound \eqn{i} in monotherapy in trial \eqn{j}, we will assume that there exist two additional parameters:
#'\eqn{\gamma_{1j}} and \eqn{\gamma_{2j}}.
#'Here, we will assume further that the new parameters can take arbitrary values, i.e. \eqn{-\infty < \gamma_{ij} < \infty}
#'and will later assume a normal distribution as prior.
#'Each of the new parameters will be assumed to shift the intercept (\eqn{log(\alpha_{ij})}) of the logistic model for
#'patients with \eqn{c=1}, i.e. the covariate is assumed to act on the intercept compared to patients with \eqn{c=0}.
#'
#'The function \code{scenario_covariate_jointBLRM()} implements two-sided and one-sided action of the covariate. More precisely:
#'\itemize{
#'\item{Two-sided:\cr Two-sided effect of the covariate is (here) interpreted in the sense that patients with \code{c=1} could
#'have a larger, lower, or identical probability of experiencing a DLT at some dose compared to patients with \code{c=0}. It is assumed
#'that we do not know which of these cases will apply, so that the covariate could increase or decrease the DLT rate of the model.}
#'\item{One-sided:\cr One-sided effect of the covariate is (here) interpreted in the sense that patients with \code{c=1} are assumed
#'to have a larger probability of experiencing a DLT compared to patients with \code{c=0}, i.e., that the model for the DLT rate
#'is enforcing the DLT rate for patients with \eqn{c=1} to be larger than for patients with \eqn{c=0}. For obvious reasons,
#'one-sided covariate action should only be applied when there is strong evidence that the patients with \eqn{c=1} must have a larger
#'DLT rate than patients with \eqn{c=0}.}}
#'
#'## One-sided one two-sided model for covariates
#'Formally, speaking, one-sided and two-sided covariates are implemented as follows in the joint BLRM with covariates.
#'Writing \eqn{\pi_j(d_i)} for the DLT rate of dose \eqn{d_i} of compound \eqn{i} in trial \eqn{j},
#'the BLRM with two-sided covariates assumes:
#'\deqn{
#'logit(\pi_j(d_i)) = log(\alpha_i) + \beta \cdot  log\left(\frac{d_i}{d^*_i}\right) + c\cdot \gamma_{ij}.
#'}
#'In other words, the DLT rate for patients with \eqn{c=0} is assumed to follow the usual logistic model, while for
#'patients with \eqn{c=1} the intercept is shifted by \eqn{\gamma_{ij}}. In particular, the intercept of the logistic
#'model could become larger or smaller, depending on the value of \eqn{\gamma_{ij}}.
#'
#'The one-sided model assumes instead
#'\deqn{
#'logit(\pi_j(d_i)) = log(\alpha_i) + \beta \cdot  log\left(\frac{d_i}{d^*_i}\right) + c\cdot exp(\gamma_{ij}).
#'}
#'That is, the shift of the intercept is assumed to be a log-normally distributed variable, or, more specifically,
#'we shift by \eqn{exp(\gamma_{ij})} for a normally distributed \eqn{\gamma_{ij}}. In particular, the logit of the DLT rate,
#'and therefore the DLT rate itself, is guaranteed to increase for \eqn{c=1} compared to \eqn{c=0}.
#'
#'## Prior specification for joint BLRM with binary covariate
#'Note that we can assume \eqn{\gamma_{ij}} to be normally distributed both for the one-sided and two-sided models.
#'This allows to include it in the hierarchical random effects model of the joint BLRM in the same fashion as the
#'other model parameters. More specifically, we consider for trial \eqn{j} the seven-dimensional parameter vector
#'\deqn{\theta_j = \theta_j = (log(\alpha_{1j}), log(\beta_{1j}), log(\alpha_{2j}), log(\beta_{2j}), \eta_j, \gamma_{1j}, \gamma_{2j}).}
#'
#'For the random effects distribution, we can proceed analogously to the case of the joint BLRM without covariates.
#'That is, for the general hierarchical prior assume that
#'\deqn{\theta_j|\mu,\Sigma \sim Normal_7(\mu, \Sigma)}
#'for a shared hyper mean vector \eqn{\mu} and a shared hyper covariance matrix \eqn{\Sigma}, with entries
#'\deqn{\mu=(\mu_1, \mu_2, \mu_3, \mu_4, \mu_5, \mu_6, \mu_7)}
#'and
#'\deqn{\Sigma = (\Sigma_{kl})}
#'for \eqn{k=1,...,7}, \eqn{l=1,...,7}. The entries \eqn{\Sigma_{kl}} are defined as
#'\deqn{\Sigma_{kk}=\tau_k\cdot\tau_k}
#'on the diagonal, and as
#'\deqn{\Sigma_{kl}=\rho_{kl}\cdot\tau_k\cdot\tau_l}
#'for \eqn{k} not equal to \eqn{l}. Here, \eqn{\rho_{kl}} are the correlation coefficients of parameters \eqn{k} and \eqn{l},
#'and \eqn{\tau_k} the standard deviations of their respective parameters. As before, we will assume all parameters to be uncorrelated,
#'except for the intercept and log-slope of the model for same compound.
#'That is, \eqn{\rho_{kl}=0} for all \eqn{kl}, except for \eqn{\rho_{12}} and \eqn{\rho_{34}}. The latter two are again assumed to be
#'uniformly distributed on \eqn{[-1, 1]}.
#'
#'The model is completed in the same way as the joint BLRM without covariate, by assuming
#'\deqn{\mu_k \sim Normal(m_{\mu_k}, s_{\mu_k}^2)}
#'as hyper-prior for \eqn{\mu_k} respectively
#'\deqn{\tau_k \sim logNormal(m_{\tau_k}, s_{\tau_k}^2)}
#'for \eqn{\tau_k}.
#'
#'Note that this model allows to use the same prior distributions as the joint BLRM without covariate for all parameters
#'\eqn{\mu_k} and \eqn{\tau_k}, except for the newly added \eqn{\mu_6}, \eqn{\mu_7}, \eqn{\tau_6}, \eqn{\tau_7}, which
#'correspond to the hyper-mean and between-trial heterogeneity for the parameters \eqn{\gamma_{ij}}. As before, we will use
#'the parameter name as subscript of \eqn{\mu_k} and \eqn{\tau_k} to indicate the corresponding parameter.
#'
#'For the function \code{scenario_covariate_jointBLRM()} and \code{sim_covariate_jointBLRM()}, the priors
#'for \eqn{\mu_6}, \eqn{\mu_7}, \eqn{\tau_6}, \eqn{\tau_7} are controlled in the arguments \code{prior.mu.covar}
#'and \code{prior.tau.covar}, while the priors for all remaining \eqn{\mu_k} and \eqn{\tau_k} are specified
#'in the arguments \code{prior.mu} and \code{prior.tau}, which use the same format and names as the joint BLRM
#'functions without covariate.
#'
#'The entries of \code{prior.mu.covar} correspond to the hyper-means for \eqn{\gamma_{ij}} as specified in the following table.
#'
#'| Entry            | \code{mu_g1}               | \code{mu_g2}               |
#'| :--------------- | :------------------:       | :------------------:       |
#'| Parameter        | \eqn{\gamma_{1j}}          | \eqn{\gamma_{2j}}          |
#'| Hyper Mean       | \eqn{\mu_6=\mu_{\gamma_1}} | \eqn{\mu_7=\mu_{\gamma_2}} |
#'
#'Similarly, the entries of \code{prior.tau.covar} follow the conventions stated below.
#'
#'| Entry            | \code{tau_g1}                | \code{tau_g2}                |
#'| :--------------- | :------------------:         | :------------------:         |
#'| Parameter        | \eqn{\gamma_{1j}}            | \eqn{\gamma_{2j}}            |
#'| Hyper Mean       | \eqn{\tau_6=\tau_{\gamma_1}} | \eqn{\tau_7=\tau_{\gamma_2}} |
#'
#'Note further that \code{prior.mu.covar} can therefore be used to control the size of the assumed difference
#'across patients with \eqn{c=0} and \eqn{c=1}. Assuming a two-sided covariate, setting e.g. \code{prior.mu.covar$mu_g1[1]=0}
#'as the mean, one essentially states that it is not known whether the covariate increases the DLT rate, while
#'\code{prior.mu.covar$mu_g1[1]>0} assumes that \eqn{c=1} is more likely to increase the DLT rate.
#'Further, when increasing the SD by setting an
#'higher value for \code{prior.mu.covar$mu_g1[2]}, one essentially assumes larger uncertainty with
#'respect to the potential influence of the covariate, allowing for greater differences across trials.
#'Setting \code{prior.mu.covar$mu_g1[2]} to a value similar to \code{prior.mu$mu_a1[2]} would
#'allow that the influence of the covariate potentially dominates the influence of the (shared) intercept
#'\eqn{log(\alpha_{ij})}, thereby reducing the degree of information sharing across cohorts with
#'\eqn{c=1} and \eqn{c=0} substantially. Hence, using the default value of \code{prior.mu$mu_a1[2]=2},
#'one will usually not want to choose \code{prior.mu.covar$mu_g1[2]} larger than about 1 to 1.5. If one
#'expects even larger uncertainty across \eqn{c=1} and \eqn{c=0}, it may be advantageous to use a different
#'model which includes this uncertainty explicitly (e.g. exchangeale-non-exchangeable-type random effects models as discussed in e.g.
#'in (Schmidli et al., 2015), (Neuenschwander et al., 2016)).
#'
#'The default
#'values for \code{prior.mu.covar} are a relatively conservative setting and assume that there is little
#'concrete knowledge about how the covariate influences the DLT rates, which is often a reasonable expectation
#'in early trials. The default values should therefore be applicable in a wide range of settings, but can be
#'adjusted and optimized in situations where there is concrete information about the influence of the covariate.
#'@returns
#'List. The output list will have at least one entry for each trial of interest that provides a summary of the posterior toxicities.
#'The summary will contain the posterior DLT rates using the covariate of interest of the corresponding trial as provided
#'in the argument \code{trials.of.interest.covars}. If the latter is not specified or set to \code{NA} for a trial, the returned list
#'will contain two entries for this trial, one for each value of the covariate.
#'If \code{output.scen.config}
#'is \code{TRUE}, additional entries that give the input data are included.
#'
#'If \code{plot.return} and \code{plot.decisions} are both \code{TRUE}, an additional entry is created that holds a list of all
#'output plots as \code{\link[ggplot2:ggplot]{ggplot2::ggplot}} objects.
#'
#'More precisely, the following list entries are always generated:
#'\itemize{
#'\item{\code{$trial-[...]_covar-[0/1]}\cr
#'Here, \code{[...]} denotes the given trial name in \code{trials.of.interest}, while
#'\code{[0/1]} is the value of the binary covariate used for the posterior. The entry
#'gives a matrix that lists summary statistics and interval probabilities based on the
#'posterior of the DLT rate of each of the
#'doses of interest for a trial. Additionally, if loss escalation is active,
#'the expected loss of each dose is listed.}
#'}
#'If additionally \code{output.scen.config} is active, there will be the following
#'additional entries:
#'\itemize{
#'\item{\code{$data}\cr
#'Input data used to fit the joint BLRM. Includes all cohorts from both \code{data}
#'and \code{historical.data} in a merged data matrix.}
#'
#'\item{\code{$prior}\cr
#'Contains the specified (hyper-)prior distribution used by the joint BLRM.}
#'
#'\item{\code{$configuration}\cr
#'Contains the remaining configurations, e.g. seed and escalation rule.}
#'
#'\item{\code{$'Stan options'}\cr
#'Contains the arguments given to Stan, e.g. number of MCMC iterations and chains.}
#'}
#'
#'If the additional plot list is generated, there will be an entry:
#'\itemize{
#'\item{\code{$plots}\cr
#'  List that contains the output plots as \code{\link[ggplot2:ggplot]{ggplot2::ggplot}} objects.
#'  There will be one entry for each trial, namely
#'  \itemize{
#'    \item{\code{$plots$trial-[...]_covar-[0/1]}\cr
#'    Here, \code{[...]} denotes the trial name of one of the trials of interest, while
#'    \code{[0/1]} is the value of the binary covariate used for the posterior. The entry for a trial
#'    is either a single \code{\link[ggplot2:ggplot]{ggplot2::ggplot}} object or a list of such objects.
#'    The number of entries is determined by the remaining specifications for plots. E.g., if loss escalation is performed,
#'    there can be either just a plot of the expected loss for each dose, or additionally a second
#'    plot with the usual interval probabilities. Similarly, if the trial of interest is of type \code{"all"},
#'    there will be up to 3 plots, one for each trial type for which dose levels are included (mono 1, mono 2, combination).}
#'  }
#'}
#'}
#'@examples
#' \dontrun{
#' result <- scenario_covariate_jointBLRM(
#'                    data=list(dose1 = c(1, 2, 4, 6, 8, 0,  0,  0,  1,  2),
#'                              dose2 = c(0, 0, 0, 0, 0, 10, 20, 30, 10, 10),
#'                              n.pat = c(3, 3, 3, 3, 3, 3,  6,  9,  3,  3),
#'                              n.dlt = c(0, 0, 0, 0, 1, 0,  0,  1,  0,  0),
#'                              trial = c(1, 1, 1, 1, 1, 2,  2,  3,  3,  3),
#'                              covar = c(0, 0, 0, 0, 0, 1,  1,  1,  1,  1)
#'                              ),
#'                    trials.of.interest =        c(1,       3),
#'                    types.of.interest =         c("mono1", "combi"),
#'                    trials.of.interest.covars = c(0,       1),
#'                    doses.of.interest = rbind(
#'                      c(1, 2, 4, 6, 8, 12,  rep(c(1, 2, 4, 6, 8, 12),
#'                                                times=3)),
#'                      c(0, 0, 0, 0, 0, 0,   rep(c(10, 20, 30),
#'                                                each = 6 ))),
#'                    dose.ref1 = 12,
#'                    dose.ref2 = 30,
#'                    esc.rule = "dynamic.loss",
#'                    prior.mu = list(mu_a1 =  c(logit(0.33), 2),
#'                                    mu_b1 =  c(0,          1),
#'                                    mu_a2 =  c(logit(0.33), 2),
#'                                    mu_b2 =  c(0,          1),
#'                                    mu_eta = c(0,          1.121)),
#'                    prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
#'                                     tau_b1 =  c(log(0.125), log(2)/1.96),
#'                                     tau_a2 =  c(log(0.25),  log(2)/1.96),
#'                                     tau_b2 =  c(log(0.125), log(2)/1.96),
#'                                     tau_eta = c(log(0.125), log(2)/1.96)),
#'                    path = getwd(),
#'                    file.name = NULL,
#'                    iter=10000,
#'                    chains=4
#'                    )
#'}
#'@references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. \url{https://mc-stan.org}
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
#' Schmidli, H., Gsteiger, S., Roychoudhury, S., O'Hagan, A., Spiegelhalter, D., Neuenschwander, B. (2014). Robust Meta-Analytic-Predictive Priors in Clinical Trials with Historical Control Information.
#' Biometrics, 70(4), 1023-1032 <doi: 10.1111/biom.12242>.
#'
#' Babb, J., Rogatko, A., & Zacks, S. (1998). Cancer phase I clinical trials: Efficient dose escalation with overdose control.
#' Statistics in medicine 17(10), 1103-1120.
#'
#' Zhou, H.,  Yuan, Y., & Nie, L. (2018). Accuracy, safety, and reliability of novel phase I designs.
#' Clinical Cancer Research, 24(21), 5483-5484 <doi: 10.1158/1078-0432.ccr-18-0168>.
#'
#' @seealso \code{\link[decider:scenario_list_covariate_jointBLRM]{scenario_list_covariate_jointBLRM}()},
#' \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()},
#' \code{\link[decider:scenario_list_jointBLRM]{scenario_list_jointBLRM}()},
#' \code{\link[decider:sim_jointBLRM]{sim_jointBLRM}()}, \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()},
#' \code{\link[rstan:rstan]{rstan-package}},
#' \code{\link[ggplot2:ggplot]{ggplot2::ggplot}}, \code{\link[ggplot2:ggplot2-package]{ggplot2-package}}.
#'
#'@md
#'@export scenario_covariate_jointBLRM

scenario_covariate_jointBLRM <- function(
                               data=NULL,
                               historical.data=NULL,
                               doses.of.interest,
                               dose.ref1,
                               dose.ref2,
                               trials.of.interest,
                               types.of.interest=NULL,
                               trials.of.interest.covars=NULL,                         #NEW
                               esc.rule=c("ewoc", "loss", "dynamic.loss"),
                               dosing.intervals = c(0.16, 0.33, 0.6),
                               ewoc.threshold = 0.25,
                               loss.weights = c(1, 0, 1, 2),
                               dynamic.weights = rbind(c(0.32, 0, 0.32, 0.36),
                                                       c(0.29, 0, 0.31, 0.4),
                                                       c(0.27, 0, 0.33, 0.4),
                                                       c(0.2,  0, 0.3,  0.5)
                                                       ),


                               prior.mu = list(mu_a1 =  c(logit(0.33), 2),
                                               mu_b1 =  c(0,          1),
                                               mu_a2 =  c(logit(0.33), 2),
                                               mu_b2 =  c(0,          1),
                                               mu_eta = c(0,          1.121)
                               ),
                               prior.mu.covar = list(mu_g1 = c(0,  1),   #NEW
                                                     mu_g2 = c(0,  1)    #NEW
                               ),
                               prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
                                                tau_b1 =  c(log(0.125), log(2)/1.96),
                                                tau_a2 =  c(log(0.25),  log(2)/1.96),
                                                tau_b2 =  c(log(0.125), log(2)/1.96),
                                                tau_eta = c(log(0.125), log(2)/1.96)
                               ),
                               prior.tau.covar = list(tau_g1 = c(log(0.125), log(2)/1.96),   #NEW
                                                      tau_g2 = c(log(0.125), log(2)/1.96)    #NEW
                               ),

                               two_sided1 = TRUE,                               #New
                               two_sided2 = TRUE,                               #New
                               saturating = FALSE,

                               probs = c(0.025, 0.5, 0.975),
                               iter = 26000,
                               warmup = 1000,
                               refresh = 0,
                               adapt_delta = 0.8,
                               max_treedepth = 15,
                               chains = 4,
                               seed=sample.int(.Machine$integer.max, 1),

                               path = NULL,
                               file.name = NULL,
                               plot.decisions = FALSE,
                               plot.combi.heatmap = TRUE,
                               plot.int.probs.loss = FALSE,
                               plot.return = FALSE,
                               plot.file.format = "pdf",
                               plot.width,
                               plot.height,
                               plot.unit,
                               output.scen.config = FALSE
                               )
{



  #--------------------------------------------------------------------------------------------------
  #Input checks and handling of special cases
  #--------------------------------------------------------------------------------------------------
  if(is.list(data)){
    names_dat <- names(data)
    if(!all(c("dose1", "dose2", "n.pat", "n.dlt", "trial", "covar")%in%names_dat)){

      stop("`data` must be a named list with entries `dose1`, `dose2`, `n.pat`, `n.dlt`, `trial`, `covar`.")
    }
    nobs <- length(data$dose1)
    if(!nobs>0){

      stop("if `data` is a list, at least one observation must be included. Otherwise, set `data=NULL`.")
    }
    if(!is.numeric(data$dose1)){

      stop("`data$dose1` must be numeric.")
    }
    if(!is.numeric(data$dose2)|!(length(data$dose2)==nobs)){

      stop("`data$dose2` must be numeric with the same length as `data$dose1`.")
    }

    if(!is.numeric(data$n.pat)|!(length(data$n.pat)==nobs)){

      stop("`data$n.pat` must be numeric with the same length as `data$dose1`.")
    }
    if(any(is.na(data$n.pat))){
      stop("`data$n.pat` may not be NA.")
    }

    if(!is.numeric(data$n.dlt)|!(length(data$n.dlt)==nobs)){

      stop("`data$n.dlt` must be numeric with the same length as `data$dose1`.")
    }
    if(any(is.na(data$n.dlt))){
      stop("`data$n.dlt` may not be NA.")
    }
    if(!is.numeric(data$trial) & !is.character(data$trial)){

      stop("`data$trial` must either be numeric or of type character.")
    }
    if(!(length(data$trial)==nobs)|any(is.na(data$trial))){

      stop("`data$trial` must be of the same length as `data$dose1` and no trial can be NA.")
    }

    if(!is.numeric(data$covar)){

      stop("`data$covar` must be numeric.")
    }
    if(!(length(data$covar)==nobs)|any(is.na(data$covar))){

      stop("`data$covar` must be of the same length as `data$dose1` and must not contain NA.")
    }
    if(!all(data$covar==1 | data$covar==0)){
      stop("Entries of `data$covar` must be either 0 or 1.")
    }

    data$dose1[which(is.na(data$dose1))] <- rep(0, times = length(which(is.na(data$dose1))))
    data$dose2[which(is.na(data$dose2))] <- rep(0, times = length(which(is.na(data$dose2))))
    data$n.dlt[which(is.na(data$n.dlt))] <- rep(0, times = length(which(is.na(data$n.dlt))))
    data$n.pat[which(is.na(data$n.pat))] <- rep(0, times = length(which(is.na(data$n.pat))))

    if(!all(data$dose1>=0) | !all(data$dose2>=0)){

      stop("negative doses detected in `data`.")
    }
    if(!all(data$n.pat>=0) | !all(data$n.dlt>=0)){

      stop("negative `n.pat` or `n.dlt` detected in `data`.")
    }
    if(!all(floor(data$n.pat)==data$n.pat) | !all(floor(data$n.dlt)==data$n.dlt)){

      stop("non-integer `n.pat` or `n.dlt` detected in `data`.")
    }
    if(!all(data$n.dlt<=data$n.pat)){

      stop("observation with `n.pat`<`n.dlt` detected in `data`.")
    }
  }

  if(!is.wholenumber(seed)){

    stop("`seed` must be a whole number.")
  }
  set.seed(seed)
  seed.internal <- sample.int(.Machine$integer.max, 1)
  set.seed(seed.internal)
  if(!is.list(historical.data)&!is.null(historical.data)){

    stop("`historical.data` must be either NULL or a list.")
  }
  if(is.list(historical.data)){
    names_dat <- names(historical.data)
    if(!all(c("dose1", "dose2", "n.pat", "n.dlt", "trial", "covar")%in%names_dat)){

      stop("`historical.data` must be a named list with entries `dose1`, `dose2`, `n.pat`, `n.dlt`, `trial`, `covar`.")
    }
    nobs <- length(historical.data$dose1)
    if(!nobs>0){

      stop("If `historical.data` is a list, at least one observation must be included.\n Otherwise, set `historical.data=NULL`.")
    }
    if(!is.numeric(historical.data$dose1)){

      stop("`historical.data$dose1` must be numeric.")
    }
    if(!is.numeric(historical.data$dose2)|!(length(historical.data$dose2)==nobs)){

      stop("`historical.data$doses` must be numeric with the same length as `historical.data$dose1`.")
    }

    if(!is.numeric(historical.data$n.pat)|!(length(historical.data$n.pat)==nobs)){

      stop("`historical.data$n.pat` must be numeric with the same length as `historical.data$dose1`.")
    }

    if(any(is.na(historical.data$n.pat))){
      stop("`historical.data$n.pat` may not be NA.")
    }


    if(!is.numeric(historical.data$n.dlt)|!(length(historical.data$n.dlt)==nobs)){

      stop("`historical.data$n.dlt` must be numeric with the same length as `historical.data$dose1`.")
    }

    if(any(is.na(historical.data$n.dlt))){
      stop("`historical.data$n.dlt` may not be NA.")
    }

    if(!is.numeric(historical.data$trial) & !is.character(historical.data$trial)){

      stop("`historical.data$trial` must either be numeric or of type character.")
    }
    if(!(length(historical.data$trial)==nobs)){

      stop("`historical.data$trial` must be of the same length as `historical.data$dose1`.")
    }
    if(any(is.na(historical.data$trial))){

      stop("`historical.data$trial` cannot contain NA.")
    }

    if(!is.numeric(historical.data$covar)){

      stop("`historical.data$covar` must be numeric.")
    }
    if(!(length(historical.data$covar)==nobs)|any(is.na(historical.data$covar))){

      stop("`historical.data$covar` must be of the same length as `historical.data$dose1` and must not contain NA.")
    }
    if(!all(historical.data$covar==1 | historical.data$covar==0)){
      stop("Entries of `historical.data$covar` must be either 0 or 1.")
    }

    historical.data$dose1[which(is.na(historical.data$dose1))] <- rep(0, times = length(which(is.na(historical.data$dose1))))
    historical.data$dose2[which(is.na(historical.data$dose2))] <- rep(0, times = length(which(is.na(historical.data$dose2))))
    historical.data$n.dlt[which(is.na(historical.data$n.dlt))] <- rep(0, times = length(which(is.na(historical.data$n.dlt))))
    historical.data$n.pat[which(is.na(historical.data$n.pat))] <- rep(0, times = length(which(is.na(historical.data$n.pat))))

    if(!all(historical.data$dose1>=0) | !all(historical.data$dose2>=0)){

      stop("negative doses detected in `historical.data`.")
    }
    if(!all(historical.data$n.pat>=0) | !all(historical.data$n.dlt>=0)){

      stop("negative `n.pat` or `n.dlt` detected in `historical.data`.")
    }
    if(!all(floor(historical.data$n.pat)==historical.data$n.pat) | !all(floor(historical.data$n.dlt)==historical.data$n.dlt)){

      stop("non-integer `n.pat` or `n.dlt`  detected in `historical.data`.")
    }
    if(!all(historical.data$n.dlt<=historical.data$n.pat)){

      stop("observation with `n.dlt`>`n.pat` detected in `historical.data`.")
    }



  }

  #reset data to combination of hist. data and data
  if(is.null(data) & is.null(historical.data)){
    message("Note: `data` and `historical.data` is NULL, only the prior will be computed.")
    data <- list(
      dose1 = as.array(c(1)),
      dose2 = as.array(c(0)),
      n.pat = as.array(c(0)),
      n.dlt = as.array(c(0)),
      trial = as.array(c(1)),
      covar = as.array(c(0))
    )
  }else if(is.null(data) & !is.null(historical.data)){
    data <- historical.data
  }else if(!is.null(data) & !is.null(historical.data)){
    #combine data and historical data
    new_dose1 <- c(historical.data$dose1, data$dose1)
    new_dose2 <- c(historical.data$dose2, data$dose2)
    new_n.pat <- c(historical.data$n.pat, data$n.pat)
    new_n.dlt <- c(historical.data$n.dlt, data$n.dlt)
    new_trial <- c(historical.data$trial, data$trial)
    new_covar <- c(historical.data$covar, data$covar)
    data <- list(
      dose1=as.array(new_dose1),
      dose2=as.array(new_dose2),
      n.pat=as.array(new_n.pat),
      n.dlt=as.array(new_n.dlt),
      trial=as.array(new_trial),
      covar=as.array(new_covar)
    )
  }


  idx_noninf_obs <- which((data$dose1==0 &
                           data$dose2==0))
  if(!length(idx_noninf_obs)==0){
    if(length(idx_noninf_obs)==length(data$dose1)){
      message("Note: Neither `data` nor `historical.data` contain cohorts to which a dose\n",
              " was administered. Non-informative data is used instead to prevent errors.")
      data <- list(
        dose1 = as.array(c(1)),
        dose2 = as.array(c(0)),
        n.pat = as.array(c(0)),
        n.dlt = as.array(c(0)),
        trial = as.array(c(1)),
        covar = as.array(c(0))
      )
    }else{
      message("Note: Detected and removed cohorts in `data` or `historical.data`\n",
              "for which both doses were 0 or NA.")
      data <- remove.noninf.obs(data, idx_noninf_obs)
    }
  }

  #in the remaining case, data is not NULL and historical data is,
  #so the data argument can stay as it is.
  nobs <- length(data$dose1)


  if(missing(doses.of.interest)){

    stop("`doses.of.interest` are missing.")
  }
  if(!is.numeric(doses.of.interest)){

    stop("`doses.of.interest` must be numeric.")
  }
  if(is.null(dim(doses.of.interest))){

    stop("`doses.of.interest` must be a matrix with two rows and at least one column.")
  }
  if(!dim(doses.of.interest)[1]==2){

    stop("`doses.of.interest` must be a matrix with two rows and at least one column.")
  }
  doses.of.interest[1, which(is.na(doses.of.interest[1,]))] <- rep(0, length(which(is.na(doses.of.interest[1,]))))
  doses.of.interest[2, which(is.na(doses.of.interest[2,]))] <- rep(0, length(which(is.na(doses.of.interest[2,]))))
  if(!all(doses.of.interest[1, ]>=0)){

    stop("`doses.of.interest` must contain non-negative entries.")
  }
  if(!all(doses.of.interest[2, ]>=0)){

    stop("`doses.of.interest` must contain non-negative entries.")
  }
  if(any(doses.of.interest[1,]==0 & doses.of.interest[2, ]==0)){

    stop("dose combination of 0/NA dose with another 0/NA dose detected in `doses.of.interest`.\nCombination of zero and zero (or NA) is not a valid dose.")
  }

  #reorder trial levels (in case user provided weird trial names for scenario_jointBLRM)
  fstudy <- factor(data$trial)
  lvstudy <- levels(fstudy)
  nstd <- length(lvstudy)
  study_proc <- rep(0, nobs)
  #for translating back in the end
  new_studies <- rep(0, nstd)
  names(new_studies) <- lvstudy
  for(level in 1:nstd){
    study_proc[which(data$trial==lvstudy[level])] <- level
    new_studies[paste0(lvstudy[level])] <- level
  }

  if(missing(trials.of.interest)){
    message("Note: no `trials.of.interest` detected, so all trials in data are used.")
    trials.of.interest <- lvstudy
  }

  if(is.null(trials.of.interest)){
    message("Note: no `trials.of.interest` detected, so all trials in data are used.")
    trials.of.interest <- lvstudy
  }

  std_interest_proc <- trials.of.interest
  is_MAP <- trials.of.interest
  n_std_interest <- length(trials.of.interest)
  MAP.prior <- FALSE
  for(i in 1:n_std_interest){
    if(paste0(trials.of.interest[i])%in%lvstudy){
      std_interest_proc[i] <- new_studies[paste0(trials.of.interest[i])]
      is_MAP[i] <- 0
    }else{
      #study not in data: use MAP prior instead
      MAP.prior <- TRUE
      is_MAP[i] <- 1
      std_interest_proc[i] <- nstd+1
    }
  }

  if(!is.null(types.of.interest)){
    if(!length(types.of.interest)==n_std_interest){

      stop("`types.of.interest` must contain one entry for each trial of interest.")
    }
    if(!all(tolower(types.of.interest)%in%c("mono1", "mono2", "combi", "all"))){

      stop("`types.of.interest` can only contain entries \"mono1\", \"mono2\", \"combi\", or \"all\".")
    }

    if(any(types.of.interest%in%c("mono1")) &
       length(
          which(doses.of.interest[1,]>0 & doses.of.interest[2,]==0)
          )==0){

      stop("`types.of.interest` is not allowed to contain \"mono1\" as the specified `doses.of.interest`\ndo not contain mono 1 doses.")
    }


    if(any(types.of.interest%in%c("mono2")) &
       length(
         which(doses.of.interest[2,]>0 & doses.of.interest[1,]==0)
       )==0){

      stop("`types.of.interest` is not allowed to contain \"mono2\" as the specified `doses.of.interest`\ndo not contain mono 2 doses.")
    }


    if(any(types.of.interest%in%c("combi")) &
       length(
         which(doses.of.interest[1,]>0 & doses.of.interest[2,]>0)
       )==0){

      stop("`types.of.interest` is not allowed to contain \"combi\" as the specified `doses.of.interest` \ndo not contain combination doses.")
    }

    # for(i in 1:n_std_interest){
    #   if(paste0(trials.of.interest[i])%in%lvstudy){
    #     type_curr <- types.of.interest[i]
    #     if(type_curr)
    #   }
    # }

  }else{

    #message(paste0("Note: trial types of the trials.of.interest are inferred from the data."))
    types.of.interest <- rep("all", times = n_std_interest)
    #try to infer type from data or set to "all" otherwise
    for(i in 1:n_std_interest){
      if(paste0(trials.of.interest[i])%in%lvstudy){
        d1_std <- data$dose1[which(paste0(data$trial)==paste0(trials.of.interest[i]))]
        d2_std <- data$dose2[which(paste0(data$trial)==paste0(trials.of.interest[i]))]
        if(all(d1_std>0 & d2_std==0)){
          types.of.interest[i] <- "mono1"
        }else if(all(d2_std>0 & d1_std==0)){
          types.of.interest[i] <- "mono2"
        }else if(all(d2_std>0 & d1_std>0)){
          types.of.interest[i] <- "combi"
        }else{
          message(paste0("Note: trial ", trials.of.interest[i], " contains observations from more than one\n",
                         "      trial type (mono 1, mono 2, combination)."))
        }
      }
    }

  }

  if(!is.null(trials.of.interest.covars)){
    if(all(!is.numeric(trials.of.interest.covars) & !is.na(trials.of.interest.covars))){
      stop("`trials.of.interest.covars` must be numeric.")
    }
    if(!length(trials.of.interest.covars)== n_std_interest & !length(trials.of.interest.covars)==1){
      stop("'trials.of.interest.covars' must have either lenght 1 or provide a value for each entry of trials.of.interest.covars")
    }
    if(length(trials.of.interest.covars)==1 & n_std_interest>1){
      trials.of.interest.covars <- rep(trials.of.interest.covars, times = n_std_interest)
    }
    if(!all(is.na(trials.of.interest.covars)|
            trials.of.interest.covars==1 |
            trials.of.interest.covars==0)){
      stop("Each entry of `trial.of.interest.covars` must be either 0, 1, or NA.")
    }

    #replace NA entries with "2" (internal code for "both")
    trials.of.interest.covars[which(is.na(trials.of.interest.covars))] <- rep(2, times = length(which(is.na(trials.of.interest.covars))))

  }else{
    #use 2 as entry, displays results with and without cov
    trials.of.interest.covars <- rep(2, times = n_std_interest)

    #we do not try to automatically detect covariates, as there might be
    #situations where patients with and without covariate are in the same study
  }


  #input data and trials are processed now.
  #continue with further input and consistency checks

  if(missing(dose.ref1)){

    stop("`dose.ref1` missing.")
  }
  if(missing(dose.ref2)){

    stop("`dose.ref2` missing.")
  }
  if(!(is.numeric(dose.ref1) & is.numeric(dose.ref2))){

    stop("`dose.ref1` and `dose.ref2` must be numeric.")
  }
  if((!dose.ref1>0) | (!dose.ref2>0)){

    stop("`dose.ref1` and `dose.ref2` must be positive.")
  }


  if(!is.character(esc.rule)){


    stop("`esc.rule` must be a character type with entries in \"ewoc\", \"loss\", or \"dynamic.loss\".")
  }
  esc.rule<-tolower(esc.rule[1])
  if(!tolower(esc.rule)%in%c("ewoc", "loss", "dynamic.loss", "dynamic")){

    stop("`esc.rule` must be a character type with entries in \"ewoc\", \"loss\", \"dynamic\" or \"dynamic.loss\".")
  }
  if(esc.rule%in%c("ewoc")){
    if(!is.numeric(ewoc.threshold)){

      stop("`ewoc.threshold` must be numeric between 0 and 1.")
    }
    if(!(0<ewoc.threshold & ewoc.threshold<=1)){

      stop("`ewoc.threshold` must be numeric between 0 and 1.")
    }
    if(!is.numeric(dosing.intervals)){

      stop("`dosing.intervals` must be numeric.")
    }
    if(!length(dosing.intervals)%in%c(1, 2, 3)){

      stop("`dosing.intervals` must be of length 1,  2, or 3")
    }
    underint.deact <- FALSE
    if(length(dosing.intervals)==1){
      dosing.intervals <- c(0, dosing.intervals[1])
      underint.deact <- TRUE
    }
    if(!(0<=dosing.intervals[1]
        & dosing.intervals[1]<dosing.intervals[2]
        & dosing.intervals[2]<=1)){

      stop(paste0("`dosing.intervals[1]` must be greater or equal to 0, and smaller than \n",
           "`dosing.intervals[2]`. Further, `dosing.intervals[2]` must be smaller or equal to 1"))
    }
    dosing.intervals.internal <- dosing.intervals[1:2]
  }
  if(esc.rule%in%c("loss")){
    if(!is.numeric(loss.weights)){

      stop("`loss.weights` must be numeric.")
    }
    if(!(length(loss.weights)==4 & all(!is.na(loss.weights))) ){

      stop("`loss.weights` must have four entries, which are not allowed to be NA.")
    }
    if(!is.numeric(dosing.intervals)){

      stop("`dosing.intervals` must be numeric.")
    }
    if(!length(dosing.intervals)==3){

      stop("`dosing.intervals` must be of length 3 for loss escalation.")
    }
    if(!(0<=dosing.intervals[1]
         & dosing.intervals[1]<dosing.intervals[2]
         & dosing.intervals[2]<dosing.intervals[3]
         & dosing.intervals[3]<=1)){

      stop(paste0("`dosing.intervals[1]` must be greater or equal to 0, and smaller than \n",
                  "`dosing.intervals[2]`.Further, `dosing.intervals[2]` must be smaller than \n",
                  "`dosing.intervals[3]`. Further, `dosing.intervals[3]` must be smaller or equal to 1"))
    }
    dosing.intervals.internal <- dosing.intervals[1:3]
  }
  if(esc.rule%in%c("dynamic.loss", "dynamic")){
    esc.rule <- "dynamic.loss"
    if(!is.numeric(dynamic.weights)){

      stop("`dynamic.weights` must be numeric.")
    }
    if(is.null(dim(dynamic.weights)) ){

      stop("`dynamic.weights` must be a 4x4 matrix.")
    }
    if(!dim(dynamic.weights)[1]==4 ){

      stop("`dynamic.weights` must be a 4x4 matrix.")
    }

    if(!dim(dynamic.weights)[2]==4 ){

      stop("`dynamic.weights` must be a 4x4 matrix.")
    }
    if(any(is.na(dynamic.weights))){

      stop("`dynamic.weights` must not contain NA.")
    }

    norm_l1 <- sum(abs(dynamic.weights[1,]))
    if(!sum(abs(dynamic.weights[2,]))==norm_l1){
      stop("`dynamic.weights` must be normalized: the sum of the abolute values of\n",
           "must be the same for each row (failed in row 2).")
    }

    if(!sum(abs(dynamic.weights[3,]))==norm_l1){
      stop("`dynamic.weights` must be normalized: the sum of the abolute values of\n",
           "must be the same for each row (failed in row 3).")
    }

    if(!sum(abs(dynamic.weights[4,]))==norm_l1){
      stop("`dynamic.weights` must be normalized: the sum of the abolute values of\n",
           "must be the same for each row (failed in row 4).")
    }
    if(!is.numeric(dosing.intervals)){

      stop("`dosing.intervals` must be numeric.")
    }
    if(!length(dosing.intervals)==3){

      stop("`dosing.intervals` must be of length 3 for loss escalation.")
    }
    if(!(0<=dosing.intervals[1]
         & dosing.intervals[1]<dosing.intervals[2]
         & dosing.intervals[2]<dosing.intervals[3]
         & dosing.intervals[3]<=1)){

      stop(paste0("`dosing.intervals[1]` must be greater or equal to 0, and smaller than \n",
                  "`dosing.intervals[2]`.Further, `dosing.intervals[2]` must be smaller than \n",
                  "`dosing.intervals[3]`. Further, `dosing.intervals[3]` must be smaller or equal to 1"))
    }
    dosing.intervals.internal <- dosing.intervals[1:3]
  }

  if(!is.null(path) & !is.character(path)){
    stop("`path` must be NULL or a character type.")
  }

  if(!is.null(file.name) & !is.character(file.name)){

    stop("`file.name` must be NULL or a character type.")
  }
  if(!is.null(file.name) & !is.null(path)){
    if(!dir.exists(file.path(path))){
      message("Note: `path` does not specify a writable directory.\nNo output will be written.")
    }
  }

  if(!is.logical(saturating)){

    stop("`saturating` must be logical.")
  }

  if(!is.logical(plot.decisions)){

    stop("`plot.decisions` must be logical.")
  }
  if(!is.logical(plot.combi.heatmap)){

    stop("`plot.combi.heatmap` must be logical.")
  }

  if(!is.logical(plot.int.probs.loss)){

    stop("`plot.int.probs.loss` must be logical.")
  }

  if(!is.logical(plot.return)){

    stop("`plot.return` must be logical.")
  }
  if(!is.character(plot.file.format)){

    stop("`plot.file.format` must be a character type.")
  }
  if(!length(plot.file.format)==1){

    stop("`plot.file.format` must be a of length 1.")
  }
  if(!tolower(plot.file.format)%in%c("pdf", "png", "jpeg")){

    stop("`plot.file.format` must be either \"pdf\" or \"png\" or \"jpeg\".")
  }
  plot.file.format <- tolower(plot.file.format)

  if(!missing(plot.unit)){
  if(!is.character(plot.unit) | !length(plot.unit)==1){

    stop("`plot.unit` must be of character type and have length 1.")
  }

  if(!tolower(plot.unit)%in%c("mm", "cm", "in")){

    stop("`plot.unit` must be one of \"mm\", \"cm\", or \"in\".")
  }
  plot.unit <- tolower(plot.unit)
  }

  if(!missing(plot.height)){
  if(!is.numeric(plot.height) | any(is.na(plot.height))){

    stop("`plot.height` must be a non-NA numeric.")
  }
  if((!length(plot.height)%in%c(1, 3)) | any(plot.height<=0)){

    stop("`plot.height` must be a single positive number or a vector of length 3.")
  }
  #plot.height <- floor(plot.height)
  }

  if(!missing(plot.width)){
  if(!is.numeric(plot.width) | any(is.na(plot.width))){

    stop("`plot.width` must be a non-NA numeric.")
  }
  if((!length(plot.width)%in%c(1, 3))  | any(plot.width<=0)){

    stop("`plot.width` must be a single positive number or a vector of length 3.")
  }
  }

  height.not.given <- missing(plot.height)
  unit.not.given <- missing(plot.unit)
  width.not.given <- missing(plot.width)
  plot.size.args.not.given <- c(height.not.given, width.not.given, unit.not.given)
  if(any(plot.size.args.not.given) & !all(plot.size.args.not.given)){


    message("Note: At least one of the arguments `plot.height`, `plot.width`, `plot.unit` was specified,")
    message("      but not all of them. Please be aware that the size of output plots is only")
    message("      affected if all three arguments are given.")
  }

  if(!is.numeric(probs)){

    stop("`probs` must be numeric.")
  }
  if(!length(probs)>=1){

    stop("`probs` must have at least one entry.")
  }
  if(!all(probs<=1 & probs>=0)){

    stop("`probs` must have entries between 0 and 1.")
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


  if(!is.list(prior.mu.covar)){
    stop("`prior.mu.covar` must be a named list.")
  }
  if(is.null(names(prior.mu.covar))){
    stop("`prior.mu.covar` must be a named list.")
  }
  names(prior.mu.covar) <- tolower(names(prior.mu.covar))
  test.prior.mu.covar(prior.mu.covar)

  if(!is.list(prior.tau.covar)){
    stop("`prior.tau.covar` must be a named list.")
  }
  if(is.null(names(prior.tau.covar))){
    stop("`prior.tau.covar` must be a named list.")
  }
  names(prior.tau.covar) <- tolower(names(prior.tau.covar))
  test.prior.tau.covar(prior.tau.covar)

  prior.mu[["mu_c1"]] <- prior.mu.covar$mu_g1
  prior.mu[["mu_c2"]] <- prior.mu.covar$mu_g2

  prior.tau[["tau_c1"]] <- prior.tau.covar$tau_g1
  prior.tau[["tau_c2"]] <- prior.tau.covar$tau_g2

  # #OLD: prior checks, kept for reference until next version.
  # if(!is.list(prior.mu)){
  #
  #   stop("prior.mu must be a list.")
  # }
  #
  # if(is.null(names(prior.mu))){
  #   stop("prior.mu must be a named list.")
  # }
  #
  # names(prior.mu) <- tolower(names(prior.mu))
  # names_pmu <- names(prior.mu)
  #
  # if(!is.list(prior.tau)){
  #
  #   stop("prior.tau must be a list.")
  # }
  #
  # if(is.null(names(prior.tau))){
  #   stop("prior.tau must be a named list.")
  # }
  # names(prior.tau) <- tolower(names(prior.tau))
  # names_ptau <- names(prior.tau)
  #
  # if(!(all(c("mu_a1", "mu_a2", "mu_b1", "mu_b2", "mu_eta")%in%names_pmu))){
  #
  #   stop("prior.mu must have entries \"mu_a1\",\"mu_b1\",\"mu_a2\",\"mu_b2\", and \"mu_eta\".")
  # }
  # if(!is.numeric(prior.mu$mu_a1)){
  #
  #   stop("prior.mu$mu_a1 must be numeric.")
  # }
  # if(any(is.na(prior.mu$mu_a1))){
  #
  #   stop("prior.mu$mu_a1 must not be NA.")
  # }
  # if(!length(prior.mu$mu_a1)==2){
  #
  #   stop("prior.mu$mu_a1 must have length 2.")
  # }
  # if(!prior.mu$mu_a1[2]>0){
  #
  #   stop("prior.mu$mu_a1[2] is the SD and cannot be negative.")
  # }
  # if(!is.numeric(prior.mu$mu_b1)){
  #
  #   stop("prior.mu$mu_b1 must be numeric.")
  # }
  # if(any(is.na(prior.mu$mu_b1))){
  #
  #   stop("prior.mu$mu_b1 must not be NA.")
  # }
  # if(!length(prior.mu$mu_b1)==2){
  #
  #   stop("prior.mu$mu_b1 must have length 2.")
  # }
  # if(!prior.mu$mu_b1[2]>0){
  #
  #   stop("prior.mu$mu_b1[2] is the SD and cannot be negative.")
  # }
  # if(!is.numeric(prior.mu$mu_a2)){
  #
  #   stop("prior.mu$mu_a2 must be numeric.")
  # }
  # if(any(is.na(prior.mu$mu_a2))){
  #
  #   stop("prior.mu$mu_a2 must not be NA.")
  # }
  # if(!length(prior.mu$mu_a2)==2){
  #
  #   stop("prior.mu$mu_a2 must have length 2.")
  # }
  # if(!prior.mu$mu_a2[2]>0){
  #
  #   stop("prior.mu$mu_a2[2] is the SD and cannot be negative.")
  # }
  #
  # if(!is.numeric(prior.mu$mu_b2)){
  #
  #   stop("prior.mu$mu_b2 must be numeric.")
  # }
  # if(any(is.na(prior.mu$mu_b2))){
  #
  #   stop("prior.mu$mu_b2 must not be NA.")
  # }
  # if(!length(prior.mu$mu_b2)==2){
  #
  #   stop("prior.mu$mu_b2 must have length 2.")
  # }
  # if(!prior.mu$mu_b2[2]>0){
  #
  #   stop("prior.mu$mu_b2[2] is the SD and cannot be negative.")
  # }
  # if(!is.numeric(prior.mu$mu_eta)){
  #
  #   stop("prior.mu$mu_eta must be numeric.")
  # }
  # if(any(is.na(prior.mu$mu_eta))){
  #
  #   stop("prior.mu$mu_eta must not be NA.")
  # }
  # if(!length(prior.mu$mu_eta)==2){
  #
  #   stop("prior.mu$mu_eta must have length 2.")
  # }
  # if(!prior.mu$mu_eta[2]>0){
  #
  #   stop("prior.mu$mu_eta[2] is the SD and cannot be negative.")
  # }
  #
  # if(!(all(c("tau_a1", "tau_a2", "tau_b1", "tau_b2", "tau_eta")%in%names_ptau))){
  #
  #   stop("prior.tau must have entries \"tau_a1\",\"tau_b1\",\"tau_a2\",\"tau_b2\", and \"tau_eta\".")
  # }
  # if(!is.numeric(prior.tau$tau_a1)){
  #
  #   stop("prior.tau$tau_a1 must be numeric.")
  # }
  # if(any(is.na(prior.tau$tau_a1))){
  #
  #   stop("prior.tau$tau_a1 must not be NA.")
  # }
  # if(!length(prior.tau$tau_a1)==2){
  #
  #   stop("prior.tau$tau_a1 must have length 2.")
  # }
  # if(!prior.tau$tau_a1[2]>0){
  #
  #   stop("prior.tau$tau_a1[2] is the SD and cannot be negative.")
  # }
  # if(!is.numeric(prior.tau$tau_b1)){
  #
  #   stop("prior.tau$tau_b1 must be numeric.")
  # }
  # if(any(is.na(prior.tau$tau_b1))){
  #
  #   stop("prior.tau$tau_b1 must not be NA.")
  # }
  # if(!length(prior.tau$tau_b1)==2){
  #
  #   stop("prior.tau$tau_b1 must have length 2.")
  # }
  # if(!prior.tau$tau_b1[2]>0){
  #
  #   stop("prior.tau$tau_b1[2] is the SD and cannot be negative.")
  # }
  # if(!is.numeric(prior.tau$tau_a2)){
  #
  #   stop("prior.tau$tau_a2 must be numeric.")
  # }
  # if(any(is.na(prior.tau$tau_a2))){
  #
  #   stop("prior.tau$tau_a2 must not be NA.")
  # }
  # if(!length(prior.tau$tau_a2)==2){
  #
  #   stop("prior.tau$tau_a2 must have length 2.")
  # }
  # if(!prior.tau$tau_a2[2]>0){
  #
  #   stop("prior.tau$tau_a2[2] is the SD and cannot be negative.")
  # }
  #
  # if(!is.numeric(prior.tau$tau_b2)){
  #
  #   stop("prior.tau$tau_b2 must be numeric.")
  # }
  # if(any(is.na(prior.tau$tau_b2))){
  #
  #   stop("prior.tau$tau_b2 must not be NA.")
  # }
  # if(!length(prior.tau$tau_b2)==2){
  #
  #   stop("prior.tau$tau_b2 must have length 2.")
  # }
  # if(!prior.tau$tau_b2[2]>0){
  #
  #   stop("prior.tau$tau_b2[2] is the SD and cannot be negative.")
  # }
  # if(!is.numeric(prior.tau$tau_eta)){
  #
  #   stop("prior.tau$tau_eta must be numeric.")
  # }
  # if(any(is.na(prior.tau$tau_eta))){
  #
  #   stop("prior.tau$tau_eta must not be NA.")
  # }
  # if(!length(prior.tau$tau_eta)==2){
  #
  #   stop("prior.tau$tau_eta must have length 2.")
  # }
  # if(!prior.tau$tau_eta[2]>0){
  #
  #   stop("prior.tau$tau_eta[2] is the SD and cannot be negative.")
  # }



  if(!is.numeric(iter)){

    stop("`iter` must be numeric.")
  }
  if(is.na(iter)){

    stop("`iter` must not be NA.")
  }
  if(!is.numeric(warmup)){

    stop("`warmup` must be numeric.")
  }
  if(is.na(warmup)){

    stop("`warmup` must not be NA.")
  }

  if(!is.numeric(refresh)){

    stop("`refresh` must be numeric.")
  }
  if(is.na(refresh)){

    stop("`refresh` must not be NA.")
  }

  if(!is.numeric(chains)){

    stop("`chains` must be numeric.")
  }
  if(is.na(chains)){

    stop("`chains` must not be NA.")
  }

  if(!is.numeric(adapt_delta)){

    stop("`adapt_delta` must be numeric.")
  }
  if(is.na(adapt_delta)){

    stop("`adapt_delta` must not be NA.")
  }

  if(!is.numeric(max_treedepth)){

    stop("`max_treedepth` must be numeric.")
  }
  if(is.na(max_treedepth)){

    stop("`max_treedepth` must not be NA.")
  }

  if(!length(max_treedepth)==1){

    stop("`max_treedepth` must have length 1.")
  }

  if(!length(iter)==1){

    stop("`iter` must have length 1.")
  }
  if(!length(warmup)==1){

    stop("`warmup` must have length 1.")
  }
  if(!length(adapt_delta)==1){

    stop("`adapt_delta` must have length 1.")
  }
  if(!length(chains)==1){

    stop("`chains` must have length 1.")
  }
  if(!length(refresh)==1){

    stop("`refresh` must have length 1.")
  }

  if(!floor(refresh)==refresh){

    stop("`refresh` must be an integer.")
  }

  if(!floor(chains)==chains){

    stop("`chains` must be an integer.")
  }
  if(!floor(iter)==iter){

    stop("`iter` must be an integer.")
  }
  if(!floor(warmup)==warmup){

    stop("`warmup` must be an integer.")
  }
  if(!floor(max_treedepth)==max_treedepth){

    stop("`max_treedepth` must be an integer.")
  }

  if(!(0.6<=adapt_delta & adapt_delta <1)){

    stop("It is not permitted to set `adapt_delta` below 0.6, or larger or equal to 1.")
  }
  if(!(refresh>=0)){

    stop("`refresh` must be larger or equal to 0.")
  }

  if(!(iter>=2000)){

    stop("`iter` must be larger or equal to 2000 (although we recommend to use even more iterations).")
  }
  #else if(iter<4000){
  #  message(paste0("Note: It is recommended to use at least 4000 iterations per chain, otherwise the\n",
  #                 "      posterior may be unreliable."))
  #}

  if(!(chains>=1)){

    stop("`chains` must be larger or equal to 1.")
  }

  if(!(warmup>=1000)){

    stop("`warmup` must be larger or equal to 1000.")
  }
  if(!(warmup<=iter-1000)){


    stop("`warmup` must be smaller or equal to `iter-1000`, to ensure that at least 1000 samples \nare saved per chain.")
  }
  # if(!(iter-warmup)>=2000){
  #   warning_msg <- paste0("Note: based on your specifications, only ", iter-warmup, " samples are drawn per chain.\n",
  #                         "      It is recommended to ensure that `iter` and `warmup` differ by at least 2000.")
  #   message(warning_msg)
  # }
  if(!max_treedepth>=10){

    stop("`max_treedepth` should be larger or equal to 10.")
  }

  if(!is.logical(two_sided1)){
    stop("`two_sided1` must be logical.")
  }

  if(!is.logical(two_sided2)){
    stop("`two_sided2` must be logical.")
  }


  # type.mono1.active <- length(which(doses.of.interest[1,]>0&doses.of.interest[2,]==0))>0
  # type.mono2.active <- length(which(doses.of.interest[2,]>0&doses.of.interest[1,]==0))>0
  # type.combi.active <- length(which(doses.of.interest[1,]>0&doses.of.interest[2,]>0))>0

  # ldi <- length(doses.of.interest[1,])
  # if(esc.rule=="dynamic.loss"){
  #   doses.of.interest2 <- matrix(NA, nrow=2, ncol=(ldi+1))
  #   doses.of.interest2[, 1:ldi] <- doses.of.interest
  #   doses.of.interest2[1, ldi+1] <- dref1
  #   doses.of.interest2[2, ldi+1] <- dref2
  # }else{
  #   doses.of.interest2 <- doses.of.interest
  # }
  #--------------------------------------------------------------------------------------------------
  #Sampling from posterior and extracting results
  #--------------------------------------------------------------------------------------------------

  res_raw_all <- post_tox_covariate_jointBLRM(
    study.interest = std_interest_proc,
    type.interest = types.of.interest,
    cov.interest = trials.of.interest.covars,
    dose1.interest=doses.of.interest[1,],
    dose2.interest=doses.of.interest[2,],
    dose1 = data$dose1,
    dose2 = data$dose2,
    dose.ref1 = dose.ref1,
    dose.ref2 = dose.ref2,
    n.pat = data$n.pat,
    n.dlt = data$n.dlt,
    n.study = study_proc,
    cov = data$covar,
    MAP.prior = MAP.prior,
    dosing.intervals = dosing.intervals.internal,
    probs = probs,
    prior.mu = prior.mu,
    prior.tau = prior.tau,
    two_sided1 = two_sided1,
    two_sided2 = two_sided2,
    saturating = saturating,
    iter = iter,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    warmup = warmup,
    refresh = refresh,
    chains = chains,
    seed=seed.internal
  )

  res_raw <- res_raw_all[["postsummaries"]]
  ref.probs.all <- res_raw_all[["refprobs"]]

  summary_list <- list()
  plot_list <- list()

  for(std in 1:n_std_interest){
    std_curr <- trials.of.interest[std]
    isMAPcurr <- is_MAP[std]
    type <- types.of.interest[std]
    covcurr <- trials.of.interest.covars[std]

    if(covcurr==2){
      cov_inds <- c(0, 1)
    }else{
      cov_inds <- covcurr
    }

    #loop over artificial vector of cov of interest
    #for current study to generate two plots
    for(covar_proc in cov_inds){

      str_curr <- paste0("trial-", std_curr, "_covar-", covar_proc)
      plotstr_curr <- paste0("trial-", std_curr, "_covar-", covar_proc)
      if(isMAPcurr==1){
        #str_curr <- paste0("Trial ", std_curr, ": MAP toxicities")
        std_curr_proc <- nstd+1
      }else{
        #str_curr <- paste0("Trial ", std_curr, ": posterior toxicities")
        std_curr_proc <- new_studies[trials.of.interest[std]]
      }
      summ.current.std_src <- res_raw[[paste0(std)]]
      #row indices of summary entries for covariate of current value
      rows_curr <- which(summ.current.std_src[, length(summ.current.std_src[1, ])]==covar_proc)
      #dims <- dim(summ.current.std_src[rows_curr, ])

      if(esc.rule=="ewoc"){
        if(!underint.deact){
          #summ.current.std_src <- res_raw[[paste0(std)]]
          #row indices of summary entries for covariate of current value
          #rows_curr <- which(summ.current.std_src[, length(summ.current.std_src[1, ])]==covar_proc)
          summ.current.std <- summ.current.std_src[rows_curr, ]
        }else{
          #summ.current.std_src <- res_raw[[paste0(std)]]
          #row indices of summary entries for covariate of current value
          #rows_curr <- which(summ.current.std_src[, length(summ.current.std_src[1, ])]==covar_proc)
          dim_curr <- dim(summ.current.std_src[rows_curr, ])

          summ.current.std <- matrix(NA, nrow=dim_curr[1], ncol = dim_curr[2]-1)
          rownames(summ.current.std) <- rownames(res_raw[[paste0(std)]])
          colnames(summ.current.std) <- c(colnames(res_raw[[paste0(std)]][, 1:(dim_curr[2]-4)]),
                                          colnames(res_raw[[paste0(std)]][, (dim_curr[2]-2):(dim_curr[2])]))

          summ.current.std[, 1:(dim_curr[2]-4)] <- res_raw[[paste0(std)]][, 1:(dim_curr[2]-4)]
          summ.current.std[, (dim_curr[2]-3):(dim_curr[2]-1)] <- res_raw[[paste0(std)]][, (dim_curr[2]-2):(dim_curr[2])]
        }
        #summ_plot <- summ.current.std
      }else{
        #summ.current.std_src <- res_raw[[paste0(std)]]
        #row indices of summary entries for covariate of current value
        #rows_curr <- which(summ.current.std_src[, length(summ.current.std_src[1, ])]==covar_proc)
        dims <- dim(summ.current.std_src[rows_curr, ])
        #perform loss or dynamic loss escalation
        #dims <- dim(res_raw[[paste0(std)]])
        summ.current.std <- matrix(NA, nrow=dims[1], ncol=(dims[2]+1))
        rownames(summ.current.std)<-rownames(res_raw[[paste0(std)]])
        if(esc.rule=="loss"){
          losstype <- "Exp.Loss"
        }else{
          losstype <- "Exp.Loss"
        }
        colnames(summ.current.std)<-c(colnames(res_raw[[paste0(std)]])[1:dims[2]-1], losstype, colnames(res_raw[[paste0(std)]])[dims[2]])
        summ.current.std[, 1:dims[2]-1] <- res_raw[[paste0(std)]][rows_curr, 1:dims[2]-1]
        summ.current.std[, dims[2]+1] <- res_raw[[paste0(std)]][rows_curr, dims[2]]

        if(esc.rule=="loss"){
          for(d in 1:dims[1]){
            intprobs_curr <- res_raw[[paste0(std)]][rows_curr[d], (2+length(probs)+1):(2+length(probs)+length(dosing.intervals)+1)]
            summ.current.std[d, dims[2]] <- sum(intprobs_curr*loss.weights)
          }
        }else{
          #calculate dynamic loss weights
          ref.p.curr <- ref.probs.all[[paste0(std)]]
          if(covar_proc==0){
            dlwm1 <- ref.p.curr[1, 1]*dynamic.weights[1, ] + ref.p.curr[1, 2]*dynamic.weights[2, ]+
              ref.p.curr[1, 3]*dynamic.weights[3, ] + ref.p.curr[1, 4]*dynamic.weights[4, ]
            dlwm2 <- ref.p.curr[2, 1]*dynamic.weights[1, ] + ref.p.curr[2, 2]*dynamic.weights[2, ]+
              ref.p.curr[2, 3]*dynamic.weights[3, ] + ref.p.curr[2, 4]*dynamic.weights[4, ]
            dlwc <- ref.p.curr[3, 1]*dynamic.weights[1, ] + ref.p.curr[3, 2]*dynamic.weights[2, ]+
              ref.p.curr[3, 3]*dynamic.weights[3, ] + ref.p.curr[3, 4]*dynamic.weights[4, ]
          }else{
            dlwm1 <- ref.p.curr[4, 1]*dynamic.weights[1, ] + ref.p.curr[4, 2]*dynamic.weights[2, ]+
              ref.p.curr[4, 3]*dynamic.weights[3, ] + ref.p.curr[4, 4]*dynamic.weights[4, ]
            dlwm2 <- ref.p.curr[5, 1]*dynamic.weights[1, ] + ref.p.curr[5, 2]*dynamic.weights[2, ]+
              ref.p.curr[5, 3]*dynamic.weights[3, ] + ref.p.curr[5, 4]*dynamic.weights[4, ]
            dlwc <- ref.p.curr[6, 1]*dynamic.weights[1, ] + ref.p.curr[6, 2]*dynamic.weights[2, ]+
              ref.p.curr[6, 3]*dynamic.weights[3, ] + ref.p.curr[6, 4]*dynamic.weights[4, ]
          }
          doses_all <- rownames(res_raw[[paste0(std)]][rows_curr, ])

          dose_str <- strsplit(doses_all, split="+", fixed = TRUE)
          for(d in 1:dims[1]){
            dose1_curr_num <- as.numeric(dose_str[[d]][1])
            dose2_curr_num <- as.numeric(dose_str[[d]][2])

            #select which type the dose is and compute dynamic loss accordingly
            if(dose1_curr_num>0&dose2_curr_num==0){
              dlw_curr <- dlwm1
            }else if(dose2_curr_num>0&dose1_curr_num==0){
              dlw_curr <- dlwm2
            }else{
              dlw_curr <- dlwc
            }
            intprobs_curr <- res_raw[[paste0(std)]][rows_curr[d], (2+length(probs)+1):(2+length(probs)+length(dosing.intervals)+1)]
            summ.current.std[d, dims[2]] <- sum(intprobs_curr*dlw_curr)
          }
        }


        #summ_plot <- summ.current.std[, c(1:(length(summ.current.std)-1))]
      }

      summary_list[[str_curr]] <- summ.current.std

      if(plot.decisions){
        if(is.null(file.name)){
          file.int.name <- NULL
        }else{
          file.int.name <- paste0(file.name, "_trial-", std_curr, "_cov-", covar_proc)
        }

        if(covar_proc==1){
          inds_reftox_curr <- 4:6
        }else{
          inds_reftox_curr <- 1:3
        }

        if(missing(plot.width)|missing(plot.height)|missing(plot.unit)){

          plot_list[[plotstr_curr]] <- plot_decisions_jointBLRM_int(
            summary=res_raw[[paste0(std)]][rows_curr, ],
            probs=probs,
            type=type,
            dosing.intervals=dosing.intervals.internal,
            esc.rule=esc.rule,
            ref.probs.all=ref.probs.all[[paste0(std)]][inds_reftox_curr, ],
            ewoc.threshold=ewoc.threshold,
            loss.weights=loss.weights,
            dloss.weights= dynamic.weights,
            combi.heatmap = plot.combi.heatmap,
            p.int.probs = plot.int.probs.loss,
            file.name=file.int.name,
            path=path,
            file.format=plot.file.format,
            underint.deact = underint.deact)
        }else{
          if(length(plot.width)==1){
            plot.width <- c(plot.width, plot.width, plot.width)
          }

          if(length(plot.height)==1){
            plot.height <- c(plot.height, plot.height, plot.height)
          }

          plot_list[[plotstr_curr]] <- plot_decisions_jointBLRM_int(
            summary=res_raw[[paste0(std)]][rows_curr, ],
            probs=probs,
            type=type,
            dosing.intervals=dosing.intervals.internal,
            esc.rule=esc.rule,
            ref.probs.all=ref.probs.all[[paste0(std)]][inds_reftox_curr, ],
            ewoc.threshold=ewoc.threshold,
            loss.weights=loss.weights,
            dloss.weights= dynamic.weights,
            combi.heatmap = plot.combi.heatmap,
            p.int.probs = plot.int.probs.loss,
            file.name=file.int.name,
            path=path,
            file.format=plot.file.format,
            width=plot.width,
            height=plot.height,
            unit=plot.unit,
            underint.deact = underint.deact)
        }

      }# end of "if(plot.decisions)"

    } # end of loop over covariates

  }   #end of loop over studies

  if(output.scen.config){
    summary_list[["data"]] <- data_matrix_covariate_jointBLRM(data)
    summary_list[["prior"]] <- prior_mat_out_cov(m=prior.mu, t=prior.tau)
    if(esc.rule%in%c("ewoc")){
      outconf <- matrix("-", nrow=5, ncol=1)
      rownames(outconf) <- c("seed","saturating", "two_sided1","two_sided2", "ewoc.threshold")
      colnames(outconf) <- c("-")
      outconf[1, 1] <- seed
      outconf[2, 1] <- saturating
      outconf[3, 1] <- two_sided1
      outconf[4, 1] <- two_sided2
      outconf[5, 1] <- ewoc.threshold
      summary_list[["configuration"]] <- outconf
    }else if(esc.rule%in%c("loss")){
      outconf <- matrix("-", nrow=5, ncol=4)
      rownames(outconf) <- c("seed","saturating", "two_sided1","two_sided2", "loss.weights")
      colnames(outconf) <- c("-", "-", "-", "-")
      outconf[1, 1] <- seed
      outconf[2, 1] <- saturating
      outconf[3, 1] <- two_sided1
      outconf[4, 1] <- two_sided2
      outconf[5, ] <- loss.weights
      summary_list[["configuration"]] <- outconf
    }else{
      outconf <- matrix("-", nrow=8, ncol=4)
      rownames(outconf) <- c("seed","saturating", "two_sided1","two_sided2", "dynamic.weights[1,]", "dynamic.weights[2,]"
                             , "dynamic.weights[3,]", "dynamic.weights[4,]")
      colnames(outconf) <- c("-", "-", "-", "-")
      outconf[1, 1] <- seed
      outconf[2, 1] <- saturating
      outconf[3, 1] <- two_sided1
      outconf[4, 1] <- two_sided2
      outconf[5, ] <- dynamic.weights[1,]
      outconf[6, ] <- dynamic.weights[2,]
      outconf[7, ] <- dynamic.weights[3,]
      outconf[8, ] <- dynamic.weights[4,]
      summary_list[["configuration"]] <- outconf
    }

    #write the input parameters and simulation options into the output list.
    options.stan <- as.matrix(c(chains, iter, warmup, adapt_delta, max_treedepth))
    rownames(options.stan) <- c('chains','iter', 'warmup','adapt_delta', 'max_treedepth')
    colnames(options.stan) <- "Value"
    summary_list[['Stan options']] <-options.stan
  }

  if(!is.null(file.name) &!is.null(path)){
    if(dir.exists(file.path(path))){
      write.xlsx(summary_list,
                 file=paste0(file.path(path), "/", file.name, ".xlsx"),
                 rowNames=TRUE, colNames=TRUE, overwrite = TRUE)
    }
  }

  if(plot.decisions & plot.return){
    summary_list[["plots"]] <- plot_list
  }

  return(summary_list)
}
#end of function


#'Process a list of data scenarios including a binary covariate in parallel
#'@rdname scenario_list_covariate_jointBLRM
#'@description Wrapper for calling \code{\link[decider:scenario_covariate_jointBLRM]{scenario_covariate_jointBLRM}()}
#'for a list of data scenarios, and processing them in parallel. This function is mostly included
#'for convenience, e.g. to evaluate a number of scenarios that are used for prior adjustment
#'in a quick manner.
#'
#'Be aware that input checks are performed by \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()}.
#'If this results in an error for one or more scenarios, the error is returned in the output list.
#'@param data.list List of hypothetical data scenarios. Each entry is passed to
#'\code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()}.
#'@param n.cores Optional integer, defaults to \code{1}. Indicates the number of workers
#'to be used during parallel processing. A value of \code{1} causes the scenarios to
#'be processed sequentially.
#'@param file.names Optional character vector, defaults to \code{NULL}. The argument
#'has only an effect when results are saved to disk, i.e. when an additional argument
#'\code{path} is supplied that states the output location.
#'If a \code{path} is supplied and \code{file.names} is \code{NULL}, the \code{file.name} of
#'the output for scenario \eqn{i} will be \code{"scenario_i"}.
#'Otherwise, a vector with the same length as the \code{data.list} can be supplied as \code{file.names},
#'which is used as a file name for each scenario.
#'@param ... Arguments that are passed to \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()},
#'e.g. priors, number of MCMC iterations, and number of chains. Can contain any arguments of
#'\code{scenario_jointBLRM()}, except for the argument \code{data}, which is taken from the \code{data.list} instead,
#'and for the argument \code{file.name}, which is taken from \code{file.names} instead.
#'Note in particular that \code{historical.data} can still be supplied, which is then included in every
#'evaluated scenario, additionally to the scenario-specific data from the \code{data.list} entries.
#'@details Be aware that parallelization could also be performed at level of the chains
#'in Stan. However, we recommend to parallelize at level of the scenarios.
#'This way, also the generation of e.g. the plots and posterior toxicities are parallelized,
#'and, additionally, the initialization of the parallel backend needs to occur only once.
#'@returns Numbered list in which each entry contains the output of \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()}
#'for the corresponding trial in the \code{data.list}.
#'
#'
#' @seealso \code{\link[decider:scenario_covariate_jointBLRM]{scenario_covariate_jointBLRM}()},\code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()},
#'  \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()},
#' \code{\link[rstan:rstan]{rstan-package}},
#' \code{\link[ggplot2:ggplot]{ggplot2::ggplot}}, \code{\link[ggplot2:ggplot2-package]{ggplot2-package}}.
#' @export
scenario_list_covariate_jointBLRM <- function(
  data.list,
  n.cores=1,
  file.names=NULL,
  ...
){

  #avoid that checks complain about no visible binding for
  #loop variable of foreach loop.
  i <- NULL

  n.scen <- length(data.list)
  if("data"%in%names(list(...))){
    stop("Argument `data` cannot be supplied in \"...\", only in the `data.list`.")
  }
  if("file.name"%in%names(list(...))){
    stop("Argument file.name cannot be supplied in \"...\", use `file.names` instead.")
  }
  if(is.null(file.names)){
    file.names <- paste0("scenario_", 1:n.scen)
  }else{
    if(!is.character(file.names) | !length(file.names)==n.scen){
      stop("`file.names` must be a character type vector of the same length as `data.list`.")
    }
  }
  if(n.cores>1){
    cl_foreach <- makeCluster(n.cores)
    registerDoParallel(cl_foreach)
  }else{
    registerDoSEQ()
  }

  res.list <- foreach(i = 1:n.scen,
                      .packages = c("decider"),
                      .errorhandling = "pass",
                      .inorder = TRUE)%dopar%
    {

      return(scenario_covariate_jointBLRM(
        data = data.list[[i]],
        file.name = file.names[i],
        ...
      ))

    }

  #clusters need to be stopped
  if(n.cores>1){stopCluster(cl_foreach)}

  return(res.list)
}


#end of file
