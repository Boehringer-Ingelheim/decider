#'@title Evaluate data from oncology dose finding using a joint BLRM
#'@description Evaluates data scenarios consisting of observations of one or more
#'monotherapy or two-drug combination therapy dose-finding trials and computes
#'posterior toxicities for a trial of interest and a set of doses of interest.
#'
#'If multiple scenarios need to be evaluated, consider using the function
#'\code{\link[decider:scenario_list_jointBLRM]{scenario_list_jointBLRM}()}
#'instead, which is a parallelized wrapper that processes a list of data scenarios
#'within the same setting via \code{scenario_jointBLRM()}.
#'
#'A description of the underlying model and methods are given in the section Details.
#'@usage
#'scenario_jointBLRM(
#'    data=NULL,
#'    historical.data=NULL,
#'    doses.of.interest,
#'    dose.ref1,
#'    dose.ref2,
#'    trials.of.interest,
#'    types.of.interest=NULL,
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
#'    prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
#'                     tau_b1 =  c(log(0.125), log(2)/1.96),
#'                     tau_a2 =  c(log(0.25),  log(2)/1.96),
#'                     tau_b2 =  c(log(0.125), log(2)/1.96),
#'                     tau_eta = c(log(0.125), log(2)/1.96)
#'    ),
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
#'will be interpreted as the same trial.}}
#'@param historical.data Optional named list, must have the same structure as \code{data}. It is equivalent to include observations
#'as \code{data} or \code{historical.data}, the latter is included solely for convenience. Trial names across data and historical data must be consistent,
#'in the sense that observations with the same entries in \code{data$trial}, respectively \code{historical.data$trial}
#'are interpreted to belong to the same trial.
#'
#'Note that the argument is only included for convenience, as it allows to perform
#'multiple inferences using the same set of historical data with changing new data, without having to include the non-changing
#'(historical) observations in the data argument.
#'@param trials.of.interest Optional vector of numerical or character trial numbers/names, for which the posterior is to be computed.
#'If missing, all trials included in the \code{data} argument are used instead, or, if there are no trials in the \code{data}
#'and \code{historical.data}, the priors for the doses of interest are computed.
#'Trial names should be consistent with the numbers given in the \code{data} argument.
#'If trial names are included that do not have observations
#'in the \code{data} argument, the predictive DLT probabilities (based on a MAP approach)
#'will be computed instead for these trials
#'(i.e., the posterior for these trials will only depend on the data observed in
#'other trials).
#'@param doses.of.interest Numeric matrix with two rows and non-negative entries. Each column gives a dose combination
#'of interest, where the first row provides the dose level of compound 1, the second row of compound 2.
#'\code{NA} entries will be interpreted as "dose was not administered", i.e. are internally converted to \code{0}.
#'Setting both doses to zero in one column is not a valid dose combination and will throw an error.
#'@param types.of.interest Optional character vector with one entry for each entry of \code{trials.of.interest} which specifies the
#'trial type of the corresponding trial of interest. Supported trial types are \code{"mono1"}, \code{"mono2"},  \code{"combi"},
#'and \code{"all"}. If \code{types.of.interest} is not provided, the trial types are inferred from the data given for the trial.
#'If the type is neither specified by the user nor uniquely inferable from the data (e.g. trials which are not included in the data yet, or trials including observations from
#'multiple types), the trial type is implicitly set to \code{"all"}.
#'
#'The available trial types specify which of the dose levels from \code{doses.of.interest} are included in the summary for
#'the corresponding trial of interest. More specifically, if the type of a trial of interest is...
#'\itemize{
#'\item{...\code{"mono1"}: Only doses of interest for which the second compound is \code{0} are included in the summary and (if activated) plots.}
#'\item{...\code{"mono2"}: Only doses of interest for which the fist compound is \code{0} are included in the summary and (if activated) plots.}
#'\item{...\code{"combi"}: Only doses of interest for which none of the compound is \code{0} are included in the summary and (if activated) plots.}
#'\item{...\code{"all"}: All doses of interest are included in the summary. If plots are activated, a plot is generated for each trial type detected within the doses of interest.}
#'}
#'If the type is \code{"all"} and the escalation rule is \code{"dynamic.loss"}, the dynamic Bayes risk (expected loss) is computed based
#'on the reference dose of the trial type corresponding to the type of trial of interest. That is, if e.g. mono 1 and combination doses
#'are included in the doses of interest and the type is \code{"all"}, the dynamic Bayes risk will select for each dose whether it
#'needs to use mono reference doses or the reference dose combination.
#'@param dose.ref1 Numeric, must be positive. Reference dose for compound 1.
#'@param dose.ref2 Numeric, must be positive. Reference dose for compound 2.
#'@param esc.rule  Optional character. Can be either \code{"ewoc"}, \code{"loss"}, \code{"dynamic"} or \code{"dynamic.loss"}, where the
#'latter two are treated synonymously. Note: If \code{"loss"}
#'or \code{"dynamic.loss"} are chosen, the parameter \code{dosing.intervals} needs to contain three entries, otherwise it can either have
#'two or three entries. If \code{esc.rule} is \code{"ewoc"}, at most the first two entries of \code{dosing.intervals} are used.
#'See the section Details for a description of the available methods.
#'@param dosing.intervals Optional numeric with 1, 2, or 3 ascending positive entries. Must have three entries when the \code{esc.rule} is
#'set to \code{"loss"}, \code{"dynamic.loss"}, or \code{"dynamic"}, otherwise (i.e. when \code{esc.rule} is \code{"ewoc"}) lengths 1, 2, or
#'3 are permitted.\cr
#'If dosing interval has only one entry (requires that \code{esc.rule} is \code{"ewoc"}), the entry is interpreted as the boundary
#'between DLT rates considered as target dosing and overdosing, i.e., no underdosing interval is considered in this case.
#'Otherwise (including the default), entries 1 and 2 are interpreted
#'as lower and upper boundary of the target dosing interval. A third entry is ignored when \code{esc.rule} is \code{"ewoc"}, and
#'for other values of \code{esc.rule} required as boundary between the dosing intervals containing the excessive and unacceptable DLT rates.
#'The argument \code{dosing.intervals} defaults to \code{c(0.16, 0.33, 0.6)}, which corresponds to using \eqn{[0.16, 0.33)} as target interval, and
#'\eqn{[0.33,0.6)} as excessive toxicity interval for loss escalation. If \code{dosing.intervals} has only one entry, say \code{x}, the computation is equivalent
#'to setting the argument to \code{c(0, x)}, but the output tables and plots will be formatted slightly differently to reflect that the underdosing
#'interval is not considered in this case.
#'@param ewoc.threshold Optional numeric between 0 and 1. Overdosing thresholds for EWOC plots. Defaults to 0.25.
#'@param loss.weights Optional numerical vector with four entries (which can be arbitrary numbers), the default is \code{c(1,0,1,2)}. Specifies the interval weights/penalties
#'that are used for static loss escalation. This is only needed when \code{esc.rule} is \code{"loss"}, and the parameter is ignored otherwise.
#'More precisely, \code{loss.weights[1]} is the weight of the underdosing interval, \code{loss.weights[2]} the weight of the  target dosing interval,
#'\code{loss.weights[3]} the weight of the excessively toxic dosing interval, and \code{loss.weights[4]} is the weight of the unacceptably toxic dosing interval.
#'See the Details section below for explanations regarding the methodology.
#'@param dynamic.weights Optional numerical matrix with four rows and four columns, and arbitrary numbers as entries.  Specifies the interval weights/penalties
#'that are used for dynamic loss escalation. This is only needed when \code{esc.rule} is \code{"dynamic"} or \code{"dynamic.loss"}, and the parameter is ignored otherwise.
#'Each row of \code{dynamic.weights} specifies one of the static loss weight vectors that are interpolated during dynamic loss escalation. See the Details section below for explanations regarding the methodology.
#'
#'More precisely, \code{dynamic.weights[1, ]} is the static weight vector that is weighted by the reference probability of underdosing,
#'\code{dynamic.weights[2, ]} is the static weight vector that is weighted by the reference probability of target dosing,
#'\code{dynamic.weights[3, ]} is the static weight vector that is weighted by the reference probability of excessively toxic dosing, and
#'\code{dynamic.weights[4, ]} is the static weight vector that is weighted by the reference probability of unacceptably toxic dosing.
#'The default value for \code{dynamic.weights} is a matrix with the following four rows (in order):
#'\itemize{
#'\item{\code{c(0.32, 0, 0.32, 0.36)},}
#'\item{\code{c(0.29, 0, 0.31, 0.40)},}
#'\item{\code{c(0.27, 0, 0.33, 0.40)},}
#'\item{\code{c(0.20, 0, 0.30, 0.50)}.}
#'}
#'@param prior.mu Optional list that gives the prior distribution for the hyper means \eqn{\mu}. The list must
#'have named entries, which all need to be numeric vectors of length 2:
#'\itemize{
#'\item{\code{prior.mu$mu_a1}\cr Numeric with length two, defaults to \code{c(logit(0.33), 2)}. Specifies mean and SD of the hypermean \eqn{\mu_{\alpha_1}} of the parameter \eqn{log(\alpha_1)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.mu$mu_b1}\cr Numeric with length two, defaults to \code{c(0, 1)}. Specifies mean and SD of the hypermean \eqn{\mu_{\beta_1}} of the parameter \eqn{log(\beta_1)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.mu$mu_a2}\cr Numeric with length two, defaults to \code{c(logit(0.33), 2)}. Specifies mean and SD of the hypermean \eqn{\mu_{\alpha_2}} of the parameter \eqn{log(\alpha_2)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.mu$mu_b2}\cr Numeric with length two, defaults to \code{c(0, 1)}. Specifies mean and SD of the hypermean \eqn{\mu_{\beta_2}} of the parameter \eqn{log(\beta_2)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.mu$mu_eta}\cr Numeric with length two, defaults to \code{c(0, 1.121)}. Specifies mean and SD of the hypermean \eqn{\mu_{\eta}} of the parameter \eqn{\eta} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'}
#'@param prior.tau Optional list that gives the prior distribution for the between-trial heterogeneities (hyper SD) \eqn{\tau}. The list must
#'have named entries, which all need to be numeric vectors of length 2:
#'\itemize{
#'\item{\code{prior.tau$tau_a1}\cr Numeric with length two, defaults to \code{c(log(0.25), log(2)/1.96)}. Specifies mean and SD on log-scale of the heterogeneity \eqn{\tau_{\alpha_1}} of the parameter \eqn{log(\alpha_1)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.tau$tau_b1}\cr Numeric with length two, defaults to \code{c(log(0.125), log(2)/1.96)}. Specifies mean and SD on log-scale of the heterogeneity \eqn{\tau_{\beta_1}} of the parameter \eqn{log(\beta_1)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.tau$tau_a2}\cr Numeric with length two, defaults to \code{c(log(0.25), log(2)/1.96)}. Specifies mean and SD on log-scale of the heterogeneity \eqn{\tau_{\alpha_2}} of the parameter \eqn{log(\alpha_2)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.tau$tau_b2}\cr Numeric with length two, defaults to \code{c(log(0.125), log(2)/1.96)}. Specifies mean and SD on log-scale of the heterogeneity \eqn{\tau_{\beta_2}} of the parameter \eqn{log(\beta_2)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.tau$tau_eta}\cr Numeric with length two, defaults to \code{c(log(0.125), log(2)/1.96)}. Specifies mean and SD on log-scale of the heterogeneity \eqn{\tau_{\eta}} of the parameter \eqn{\eta} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'}
#'@param saturating Optional logical, defaults to \code{FALSE}. If \code{TRUE}, the BLRM will be using a saturating interaction term as described in
#'\code{\link[OncoBayes2:blrm_formula_saturating]{OncoBayes2::blrm_formula_saturating}()}. Also refer to the Details section below.
#'@param probs Optional numeric with arbitrary entries between 0 and 1. Provides the levels for the quantiles displayed in
#'the output. Defaults to \code{c(0.025, 0.5, 0.975)}.
#'@param path Optional character that specified the path to save the resulting output files. If \code{NULL} (the default), no
#'output is written (but still returned to R). Otherwise, it is checked whether \code{path} specifies a directory, and, if yes,
#'all output is saved there.
#'@param file.name Optional name for the output file. If \code{NULL} or missing, no output is saved (and also if no valid path is given). Results are only returned to R in this case. Note
#'that plots will not be returned to R unless the additional argument \code{plot.return} is specified.
#'@param plot.decisions Optional logical, defaults to \code{FALSE}. If \code{TRUE}, plots of escalation decisions according to the
#'specified escalation rule are created. If \code{plot.decisions} is \code{TRUE}, it is recommended to supply a writable \code{path} and valid \code{file.name},
#'as plots will by default only be saved and not returned to R. Plots are only returned to R if the additional argument \code{plot.return} is specified. In particular,
#'if \code{plot.decisions} is \code{TRUE} and \code{plot.return} is \code{FALSE}, the plots will not be accessible unless a valid \code{file.name}
#'is given to enable saving the plots.
#'
#'Depending on the value of \code{esc.rule}, the following plots are created:
#'\itemize{
#'\item{If \code{esc.rule=="ewoc"}, a diagram of the posterior probabilities that the DLT rate lies in the underdosing, target dosing and overdosing intervals
#'is created, which indicates the boundary given in \code{ewoc.threshold} with red color coding (and, therefore, implicitly the decision of EWOC-based
#'dose escalation rules). For trials of types \code{"mono1"} and \code{"mono2"} (cf. argument \code{types.of.interest}), a bar plot is created, and for trials of type \code{"combi"},
#'either a heatmap of the overdosing probabilities can be created (if the value of \code{plot.combi.heatmap} is \code{TRUE}) or alternatively
#'a bar plot in the same fashion as monotherapy trials.}
#'\item{If \code{esc.rule=="loss"}, a diagram of the expected loss (Bayes risk) based on the weight choice in \code{loss.weights}
#'is created. Further, the argument \code{plot.int.probs.loss} can be used to additionally create plots for the
#'underdosing, target dosing and overdosing intervals. For trials of types \code{"mono1"} and \code{"mono2"} (cf. argument \code{types.of.interest}), a bar plot is created, and for trials of type \code{"combi"},
#'either a heatmap of the expected loss can be created (if \code{plot.combi.heatmap==TRUE}) or alternatively
#'a bar plot in the same fashion as monotherapy trials.}
#'\item{If \code{esc.rule=="dynamic.loss"}, a diagram of the expected dynamic loss (dynamic Bayes risk) based on the weight choice in \code{dynamic.weights}
#'is created. Further, the argument \code{plot.int.probs.loss} can be used to additionally create plots for the
#'underdosing, target dosing and overdosing intervals. For trials of types \code{"mono1"} and \code{"mono2"} (cf. argument \code{types.of.interest}), a bar plot is created, and for trials of type \code{"combi"},
#'either a heatmap of the expected dynamic loss can be created (if \code{plot.combi.heatmap==TRUE}) or alternatively
#'a bar plot in the same fashion as monotherapy trials.}
#'}
#'@param plot.combi.heatmap Optional logical, defaults to \code{TRUE}. If the value is \code{TRUE}, combination therapy plots
#'are created as heatmaps instead of bar plots. This affects all escalation rules.
#'@param plot.int.probs.loss Optional logical, defaults to \code{FALSE}. Only has an effect if \code{esc.rule} is either
#'\code{"loss"} or \code{"dynamic.loss"}. In this case, if the value is \code{TRUE}, additional plots will
#'be created that display the interval probabilities to complement the (always created) plots that display the resulting expected loss.
#'@param plot.return Optional logical, defaults to \code{FALSE}. If set to \code{TRUE}, the functions return the created plots
#'to R in the result list. In this case, the result list will have an additional entry, \code{...$plots}, which contains the posterior plots
#'as a named list. The returned plots will be given as a list of \code{\link[ggplot2:ggplot]{ggplot2::ggplot}()} objects, which can be displayed in R
#'via the \code{\link[ggplot2:print.ggplot]{ggplot2::plot.ggplot}} method. Refer to the documentation of the \code{\link[ggplot2:ggplot2-package]{ggplot2-package}}
#'package for further detail.
#'@param plot.file.format Optional character, defaults to \code{"pdf"}. Can either be \code{"pdf"}, \code{"jpeg"} or \code{"png"}. Specifies the file format
#'of the created output plots (if plots are saved). Note that \code{"pdf"} will use \code{\link[grDevices:cairo]{grDevices::cairo_pdf}()}.
#'@param plot.unit Optional character string, can be "in", "cm", or "mm". Provides the unit of measurements for \code{plot.width} and \code{plot.height}. Has only an effect if all of
#' the parameters \code{plot.unit}, \code{plot.width}, and \code{plot.height} are defined. In this case, The output plots will be created with the
#' specified width and length (interpreted in the given unit of measurement). If one of these parameters is missing, the default device size is used.
#' @param plot.width Optional numerical value or vector, can have length 1 or 3 and must have positive entries.
#' Provides the width of the output plots measured in the unit given in \code{plot.unit}.
#' If \code{plot.width} is a single value, all plots will use the same width. Otherwise, if \code{plot.width} has length 3, the entry \code{plot.width[1]}
#' is used as width for plots for trials of type \code{"mono1"}, while \code{plot.width[2]} is used for trials of type \code{"mono2"} and \code{plot.width[3]}
#' for trials of type \code{"combi."}.\cr
#' \code{plot.width} has only an effect if all of
#' the parameters \code{plot.unit}, \code{plot.width}, and \code{plot.height} are defined. In this case, The output plots will be created with the
#' specified width and length (interpreted in the given unit of measurement). If one of these parameters is missing, the default device size is used.
#' @param plot.height Optional numerical value or vector, can have length 1 or 3 and must have positive entries.
#' Provides the height of the output plots measured in the unit given in \code{plot.unit}.
#' If \code{plot.height} is a single value, all plots will use the same height. Otherwise, if \code{plot.height} has length 3, the entry \code{plot.height[1]}
#' is used as height for plots for trials of type \code{"mono1"}, while \code{plot.height[2]} is used for trials of type \code{"mono2"} and \code{plot.height[3]}
#' for trials of type \code{"combi."}. \cr
#' \code{plot.height} has only an effect if all of
#' the parameters \code{plot.unit}, \code{plot.width}, and \code{plot.height} are defined. In this case, The output plots will be created with the
#' specified width and length (interpreted in the given unit of measurement). If one of these parameters is missing, the default device size is used.
##'@param iter Optional integer, number of total MCMC iterations per chain. Defaults to 26000. Note: Number of warmup iterations is counted towards \code{iter}, i.e.
#'of the \code{iter} many iterations, the first \code{warmup} many samples are not saved.
#'@param warmup Optional integer, number of warmup iterations discarded from total MCMC iterations per chain. Defaults to \code{1000}.
#'Note: Number of warmup iterations is counted towards \code{iter}, i.e. of the \code{iter} many iterations, the first \code{warmup} many samples
#'are not saved. Refer to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()} from the \code{\link[rstan:rstan]{rstan-package}} package for details.
#'@param chains Optional integer. Number of Markov chains constructed by Stan. Defaults to 4.
#'@param refresh Optional integer. Given to Stan's \code{refresh} argument for \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()}, defaults to \code{0}. A positive value indicates that updates on sampling progress are
#'printed every \code{refresh} many iterations. Set this to 0 to avoid output from Stan. Refer to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()} from the \code{\link[rstan:rstan]{rstan-package}} package for details.
#'@param adapt_delta Optional numeric between 0.6 and 1, default is 0.8. Given to Stan's \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()} method in the \code{control} argument of Stan. Translates to target Metropolis
#'acceptance probability during Hamiltonian Monte Carlo, respectively NUTS.
#'Can be used to influence the \code{stepsize} Stan uses for leapfrog steps during the NUTS algorithm. See \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()}, \code{\link[rstan:rstan]{rstan-package}} for details.
#'Note: One can increase the value of this parameter
#'in order to reduce the number of Stan warnings on "divergent transition". Larger values than 0.95 will slow down sampling. The function will not
#'permit to set \code{adapt_delta} to values below 0.6 (the \code{\link[rstan:rstan]{rstan-package}} default).
#'@param max_treedepth Optional integer, defaults to 15. Should not be lower than 10. This is a parameter of the NUTS algorithm. Roughly speaking, NUTS constructs
#'a search tree when generating a new proposal, until a stopping criterion (the NUTS criterion) is satisfied or, alternatively, until
#'the maximum depth of the search tree is reached (to avoid endless looping). The argument \code{max_treedepth} allows to control the latter.
#'The maximum treedepth will also automatically constrain the maximum number of leapfrog steps, and should therefore not be below 10.
#'See also \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()} for details on \code{max_treedepth}.
#'Note: The maximum treedepth should usually not be reached during sampling, otherwise Stan will print a warning.
#'In this case, it is recommended to increase its value.
#'@param seed Optional positive integer that specifies the seed to be used for the simulation. The default is \code{sample.int(.Machine$integer.max, 1)}.
#'Note that reproducibility can only be obtained when the function is executed on exactly the same computing architecture, using identical software versions
#'(e.g. of the compiler, Stan, and R), and the same input specifications. This is due to the internal use of \code{\link[rstan:rstan]{rstan-package}} for MCMC sampling, which is
#'only reproducible under these restrictions (refer to the Stan reference manual from \code{\link[rstan:rstan]{rstan-package}}'s \href{https://mc-stan.org}{homepage} for more detail).
#'@param output.scen.config Optional logical, defaults to \code{FALSE}. If \code{TRUE}, the results will contain the necessary specifications of
#'the evaluated scenario, i.e., the data scenario, priors, escalation rule specification, and seed. Otherwise, only the posterior summaries
#'are returned.
#'@details
#'The joint BLRM is defined according to (Neuenschwander et al., 2014 and 2016). It allows to perform Bayesian logistic regression
#'to estimate the dose-toxicity relationship of two different monotherapies and combination therapy with these compounds in
#'a joint model, which includes a hierarchical prior for robust borrowing across trials. The following gives a brief introduction to
#'the model, for more detail refer to the references provided at the end of this documentation page.
#'
#'## Model specification
#'The model assumes, first of all, that treatment with a dose or dose combination \eqn{d} is during trial \eqn{j} governed by
#'its toxicity parameter \eqn{\pi_j(d)}, so that
#'\deqn{r|\pi_j(d) \sim Binomial(n, \pi_j(d))}
#'for the number of DLTs, \eqn{r}, observed after treating a cohort of \eqn{n} patients at dose \eqn{d} during trial \eqn{j}.
#'The toxicity is then modeled based on the type of treatment. If \eqn{d_i} is a dose of compound \eqn{i}, it is assumed that
#'\deqn{logit(\pi_j(d)) = log(\alpha_{ij}) + \beta_{ij}\cdot log\left(\frac{d_i}{d^*_i}\right)}
#'for unknown parameters \eqn{\alpha_{ij}>0} and \eqn{\beta_{ij}>0} dependent on the trial \eqn{j} and compound \eqn{i}. Moreover,
#'\eqn{d^*_i} denotes a fixed reference dose for compound \eqn{i}.
#'
#'Further, if \eqn{d=(d_1, d_2)} is a dose combination of compounds 1 and 2, write
#'\eqn{\pi_j(d_i)} for their toxicity according to the monotherapy model, and define the
#'interaction-free toxicity in trial \eqn{j} as
#'\deqn{\pi_{0j}(d_1, d_2) = \pi_j(d_1) + \pi_j(d_2) - \pi_j(d_1)\cdot \pi_j(d_2).}
#'Based on this, an additional interaction parameter, \eqn{\eta_j}, is introduced,
#'and the combination therapy toxicity with interaction is modeled as
#'\deqn{logit(\pi_j(d_1, d_2)) = logit(\pi_{0j}(d_1, d_2)) + \eta_j \cdot \frac{d_1}{d^*_1} \cdot \frac{d_2}{d^*_2}.}
#'This concludes the definition of the likelihood of the model.
#'In the following, we denote the interaction term as the function
#'\deqn{
#'g(\eta, d_1, d_2) = \eta \cdot \frac{d_1}{d^*_1} \cdot \frac{d_2}{d^*_2}}.
#'This is the so-called linear (on logit scale) interaction term.
#'
#'### Saturating interaction term
#'Another possibility for modeling the interaction is the so-called
#'saturating interaction term, which performs additional scaling with the aim of
#'addressing potentially problematic behavior of the linear interaction model in case
#'of large dose levels. More specifically, if one has \eqn{\eta<0} and considers the
#'limit \eqn{d_i\rightarrow \infty} (for one or both \eqn{i}), it holds
#'\deqn{
#'g(\eta, d_1, d_2) \rightarrow -\infty
#'}
#'and, consequently, \eqn{\pi_j(d_1, d_2) \rightarrow 0}. Although this seems to
#'be mostly affecting doses larger than the reference dose, this property may not
#'be considered realistic, so that one may want to use modified models that
#'do not show this behavior. For this aim, the function allows to use a
#'saturating interaction term instead, which follows the ideas from \code{\link[OncoBayes2:blrm_formula_saturating]{OncoBayes2::blrm_formula_saturating}()},
#'to avoid the aforementioned phenomenon.
#'
#'The saturating interaction term is based on the usual interaction term, but
#'includes an additional scaling factor in dependence of the dose level, which ensures
#'a finite limit when letting \eqn{d\rightarrow \infty}. More precisely, the
#'linear interaction term,\eqn{g(\eta, d_1, d_2)},
#'is exchanged with the saturating term
#'\deqn{
#'\tilde g (\eta, d_1, d_2) = g(\eta, d_1, d_2) \cdot  \frac{2}{1+\frac{d_1}{d^*_1} \cdot \frac{d_2}{d^*_2}}
#'}
#'The scaling factor is defined so that whenever \eqn{d_i=d^*_i} for both \eqn{i},
#'the saturating interaction term is equal to the usual interaction term, while
#'it holds \eqn{\tilde g (\eta, d_1, d_2)\rightarrow 2\eta} for \eqn{d_i\rightarrow \infty}.
#'
#'## Prior structure and prior choice
#'Regarding the prior configuration of the joint BLRM, denote by
#'\deqn{\theta_j=(log(\alpha_{1j}), log(\beta_{1j}), log(\alpha_{2j}), log(\beta_{2j}), \eta_j)}
#'the trial-specific parameter vector for trial \eqn{j}.
#'For the general hierarchical prior assume that
#'\deqn{\theta_j|\mu,\Sigma \sim Normal_5(\mu, \Sigma)}
#'for a shared hyper mean vector \eqn{\mu} and a shared hyper covariance matrix \eqn{\Sigma}, with entries
#'\deqn{\mu=(\mu_1, \mu_2, \mu_3, \mu_4, \mu_5)}
#'and
#'\deqn{\Sigma = (\Sigma_{kl})}
#'for \eqn{k=1,...,5}, \eqn{l=1,...,5}. The entries \eqn{\Sigma_{kl}} are defined as
#'\deqn{\Sigma_{kk}=\tau_k\cdot\tau_k}
#'on the diagonal, and as
#'\deqn{\Sigma_{kl}=\rho_{kl}\cdot\tau_k\cdot\tau_l}
#'for \eqn{k} not equal to \eqn{l}. Here, \eqn{\rho_{kl}} are the correlation coefficients of parameters \eqn{k} and \eqn{l},
#'and \eqn{\tau_k} the standard deviations of their respective parameters. As the \eqn{\tau_k} control the standard deviation
#'across the parameters of different trials, these parameters are referred to as between-trial heterogeneities.
#'
#'Regarding the prior for the model, it suffices to specify the distributions of \eqn{\mu_k}, \eqn{\tau_k} and \eqn{\rho_{kl}}.
#'As a simplification, it is further assumed that all entries of \eqn{\theta_j} are uncorrelated, except for \eqn{log(\alpha_{ij})}
#'and \eqn{log(\beta_{ij})} (i.e., the parameters belonging to the same monotherapy compound).
#'That is, \eqn{\rho_{kl}=0} for all \eqn{kl}, except for \eqn{\rho_{12}} and \eqn{\rho_{34}}.
#'As hyperprior distribution of the correlations, the two remaining correlation coefficients \eqn{\rho_{12}} and \eqn{\rho_{34}}
#'are assumed to follow \eqn{Uniform(-1, 1)} distributions. Note that this cannot be changed in the functions \code{scenario_jointBLRM()},
#'\code{sim_jointBLRM()}, and \code{fit_jointBLRM()}.
#'
#'For the hyper means \eqn{\mu_k}, it is assumed that
#'\deqn{\mu_k \sim Normal(m_{\mu_k}, s_{\mu_k}^2)}
#'for fixed mean \eqn{m_{\mu_k}} and standard deviation \eqn{s_{\mu_k}}. These values are specified in the \code{prior.mu} argument of \code{scenario_jointBLRM()} and \code{sim_jointBLRM()}.
#'For simplicity, the entries of the list \code{prior.mu} are vectors of length two, which contain the fixed mean \eqn{m_{\mu_k}}
#'and fixed SD \eqn{s_{\mu_k}}. These entries are named for their corresponding parameters, as given in the following overview. For reasons of
#'readability, we will use the corresponding parameter names as subscript by writing e.g. \eqn{\mu_{\alpha_1}} for the hypermean of \eqn{log(\alpha_1)}.
#'This is summarized in the following table.
#'
#'| Entry            | \code{mu_a1}               | \code{mu_b1}              | \code{mu_a2}               | \code{mu_b2}               | \code{mu_eta}           |
#'| :--------------- | :------------------:       | :------------------:      | :------------------:       | :------------------:       | :------------------:    |
#'| Parameter        | \eqn{log(\alpha_{1j})}     | \eqn{log(\beta_{1j})}     | \eqn{log(\alpha_{2j})}     | \eqn{log(\beta_{2j})}      | \eqn{\eta_j}            |
#'| Hyper Mean       | \eqn{\mu_1=\mu_{\alpha_1}} | \eqn{\mu_2=\mu_{\beta_1}} | \eqn{\mu_3=\mu_{\alpha_2}} | \eqn{\mu_4==\mu_{\beta_2}} | \eqn{\mu_{\eta}}        |
#'
#'For the hyper SDs \eqn{\tau_k}, it is assumed that
#'\deqn{\tau_k \sim logNormal(m_{\tau_k}, s_{\tau_k}^2)}
#'for fixed mean \eqn{m_{\tau_k}} and standard deviation \eqn{s_{\tau_k}}, which are both to be given on log-scale (that is,
#'\eqn{\tau_k=exp(y)} for a random variable \eqn{y} which has a \eqn{Normal(m_{\tau_k}, s_{\tau_k}^2)} distribution).
#'These values are specified in the \code{prior.tau} argument of \code{scenario_jointBLRM()} and \code{sim_jointBLRM()}.
#'For simplicity, the entries of the list \code{prior.tau} are vectors of length two, which contain the fixed mean \eqn{m_{\tau_k}}
#'and fixed SD \eqn{s_{\tau_k}} on log-scale. These entries are named for their corresponding parameters, as given in the following overview.
#'
#'| Entry            | \code{tau_a1}                | \code{tau_b1}               | \code{tau_a2}                | \code{tau_b2}               | \code{tau_eta}           |
#'| :--------------- | :------------------:         | :------------------:        | :------------------:         | :------------------:        | :------------------:     |
#'| Parameter        | \eqn{log(\alpha_{1j})}       | \eqn{log(\beta_{1j})}       | \eqn{log(\alpha_{2j})}       | \eqn{log(\beta_{2j})}       | \eqn{\eta_j}             |
#'| Hyper SD         | \eqn{\tau_1=\tau_{\alpha_1}} | \eqn{\tau_2=\tau_{\beta_1}} | \eqn{\tau_3=\tau_{\alpha_2}} | \eqn{\tau_4=\tau_{\beta_2}} | \eqn{\tau_5=\tau_{\eta}} |
#'
#'See (Neuenschwander et al., 2014 and 2016) for details regarding the prior choice.
#'
#'## Escalation rules
#'The function supports in general three different escalation rules and corresponding displays of the results.
#'In this context, the term "escalation rule" refers to a method that is used to derive concrete dosing recommendations
#'from the posterior computed by the BLRM. The supported escalation rules are the escalation with overdose control (EWOC)
#'principle, static loss escalation, and dynamic loss escalation.
#'
#'### Escalation with overdose control (EWOC)
#'The EWOC principle goes back to Babb et al. (1998) and derives dosing recommendations based on the posterior probability for overdosing, where overdosing refers
#'to DLT rates above the specified target interval. By default, the overdosing interval starts therefore at 0.33.
#'Formally, EWOC considers the posterior distribution \eqn{\pi_j(d)} of the DLT rate for a dose \eqn{d} and demands that
#'the posterior probability of \eqn{\pi_j(d)} lying in the overdosing interval is lower than a prespecified feasibility bound
#'(which defaults to 0.25). As a formula, the EWOC condition is therefore for a dose \eqn{d}  during trial \eqn{j} (using the default values 0.33 and 0.25
#'for lower boundary of the overdosing interval and feasibility bound):
#'\deqn{
#'Prob(\pi_j(d) \ge 0.33) < 0.25.
#'}
#'Only doses that satisfy this criterion are then considered as recommendations for the next dose. To obtain a concrete dose among
#'the remaining levels, one typically demands an additional rule. Most commonly, the "optimal probability rule" (selecting the dose
#'that has the largest probability of having a DLT rate within the target interval among the levels that satisfy EWOC) or the "maximal dose rule" (selecting the
#'largest dose that satisfies EWOC) are chosen. For dose combinations, one of the compounds should additionally be selected to be
#'maximized first in case of draws.
#'
#'### Loss function-based escalation
#'Static loss escalation goes back to Neuenschwander et al. (2008) and applies Bayesian decision theory to obtain dosing recommendations.
#'That is, the selection of a dosing recommendation is viewed as a decision problem among different available dose levels and different decisions (dose levels)
#'are evaluated using a formal loss function. To define the loss function, the range of DLT rates is divided in four disjoint dosing intervals, namely,
#'underdosing (interval \eqn{I_1}), target dosing (\eqn{I_2}), excessively toxic dosing (\eqn{I_3}), and unacceptably toxic dosing (\eqn{I_4}). In contrast to EWOC,
#'this corresponds to dividing the overdosing interval into two distinct intervals to achieve a finer differentiation.
#'
#'To define a formal loss function based on these intervals, a pre-specified penalty value (or interval weight), \eqn{w_i}, is assigned to each of the dosing intervals.
#'The loss function is then defined to be \eqn{w_i} for the decision on a dose with DLT rate in interval \eqn{I_i}. As it is desired to select
#'doses with DLT rate in the target interval, one should typically assume that \eqn{w_2=0} and that the remaining penalties are positive numbers.
#'The dosing recommendation is now achieved by calculating for each fixed dose level the expected value of the loss function, the so-called Bayes risk
#'(or expected loss). By definition of the loss function, the Bayes risk of a dose \eqn{d} during trial \eqn{j} is
#'\deqn{
#'\sum_{i=1}^4 w_i \cdot Prob(\pi_j(d) \in I_i).
#'}
#'The recommended decision is then the dose with minimal Bayes risk among the possible doses for the next cohort.
#'
#'In terms of the pre-specified weights \eqn{w_i}, one should always set \eqn{w_2 = 0} for the target interval, but for the remaining weights
#'manual tuning may be required to cover different situations. Values that often work are
#'e.g. \eqn{w_1 = 1}, \eqn{w_3 = 1}, \eqn{w_4 = 2} for moderate to aggressive escalation decisions, or  \eqn{w_1 = 1}, \eqn{w_3 = 2}, \eqn{w_4 = 3}
#'for rather conservative escalation decisions. Note that \eqn{w_1 = 1}, \eqn{w_3 = 1}, \eqn{w_4 = 1} leads to a decision rule that simply maximizes the
#'probability to lie in the target interval, which was found to lead to overly aggressive decisions in simulations according to Zhou et al. (2018).
#'
#'### Dynamic adjustment of loss function-based escalation
#'Dynamic loss escalation refers to a heuristic that aims to adapt the pre-specified loss weights of the method described in the previous subsection.
#'The motivation for this are potential difficulties with respect to pre-specifying weights that ensure appropriate behavior
#'in different situations (larger or lower true DLT rates), as well as the attempt to obtain an escalation rule that behaves
#'more robustly across scenarios, in the sense that it obtains good performance in settings of very low DLT rates as well as in settings of very
#'large ones. Static loss escalation often needs to be tuned to one of these settings to obtain good on-trial behavior, but may perform poorly in other settings due to this.
#'Dynamic loss escalation tries to avoid this by adapting the loss weights to more aggressive or more conservative configurations
#'using a heuristic based on the observed data and posterior.
#'
#'To apply dynamic loss escalation, four fixed loss weight combinations
#'\deqn{W_k=(w_{k1}, w_{k2}, w_{k3}, w_{k4})}
#'for \eqn{k=1,...,4} with different degrees of aggressiveness need to be selected. The weights should be normalized
#'in the sense that the absolute values of each \eqn{W_k} must sum up to the same number (i.e., have identical \eqn{l_1} norms).
#'As decisions from static loss functions are invariant under linear rescaling of the weights, one can assume without loss of generality
#'that the weights sum up to 1.
#'
#'For the weight adjustment, new loss weights are computed prior to every decision by interpolating between the fixed weights \eqn{W_k}
#'based on the current posterior interval probabilities of the reference dose of the trial of interest.
#'Denote the reference dose as \eqn{d^*}, which can either be a monotherapy dose level (for monotherapy trials) or a dose combination.
#'Now, for each dosing interval \eqn{I_i}, the corresponding reference interval probabilities \eqn{p^*_i} are computed from the posterior of the BLRM, i.e.,
#'assuming that trial \eqn{j} is the trial of interest, compute for \eqn{i=1,...,4}
#'\deqn{
#'p^*_i = Prob(\pi_j(d^*) \in I_i).
#'}
#'Using these, define the current dynamic weight \eqn{w^*} as the interpolation
#'\deqn{
#'w^* = p^*_1 * W_1 + p^*_2 * W_2 + p^*_3 * W_3 + p^*_4 * W_4.
#'}
#'These dynamic weights can now used to define a formal loss function and arrive at a decision by computing the dose with
#'the minimal expected loss based on the dynamic penalty weights for the intervals.
#'
#'Intuitively, the heuristic assumes that the reference dose is a priori assumed to be a likely candidate for the MTD, and typically among the
#'larger planned dose levels. That is, when the reference dose has a large probability of being an underdose, the dose-toxicity scenario
#'was less toxic than expected a priori, and a more aggressive weight combination can be applied. Conversely, if the reference dose has a large
#'probability of being overly toxic, the planned dose levels (and therefore, the investigational treatment) may likely be more toxic
#'than expected a priori, so that a more conservative configuration of interval penalties (larger penalty for overdosing intervals) is appropriate.
#'Due to this, the weights \eqn{W_1} should represent the most aggressive configuration, \eqn{W_2} and \eqn{W_3} should be moderate to conservative
#'configurations, while \eqn{W_4} should be the most conservative among the prespecified weights. Refer to the parameter documentation of
#'\code{dynamic.weights} for a proposed default configuration that was found to work well in many situations.
#'
#'@returns
#'List. The output list will have an entry for each trial of interest that provides a summary of the posterior toxicities. If \code{output.scen.config}
#'is \code{TRUE}, additional entries that give the input data are included.
#'
#'If \code{plot.return} and \code{plot.decisions} are both \code{TRUE}, an additional entry is created that holds a list of all
#'output plots as \code{\link[ggplot2:ggplot]{ggplot2::ggplot}} objects.
#'
#'More precisely, the following list entries are always generated:
#'\itemize{
#'\item{\code{$trial-[...]}\cr
#'Here, \code{[...]} denotes the given trial name in \code{trials.of.interest}. The entry
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
#'    \item{\code{$plots$trial-[...]}\cr
#'    Here, \code{[...]} denotes the trial name of one of the trials of interest. The entry for a trial
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
#' result <- scenario_jointBLRM(
#'                    data=list(dose1 = c(1, 2, 4, 6, 8, 0,  0,  0,  1,  2),
#'                              dose2 = c(0, 0, 0, 0, 0, 10, 20, 30, 10, 10),
#'                              n.pat = c(3, 3, 3, 3, 3, 3,  6,  9,  3,  3),
#'                              n.dlt = c(0, 0, 0, 0, 1, 0,  0,  1,  0,  0),
#'                              trial = c(1, 1, 1, 1, 1, 2,  2,  3,  3,  3)
#'                              ),
#'                    trials.of.interest = c(1,       3),
#'                    types.of.interest =   c("mono1", "combi"),
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
#' Babb, J., Rogatko, A., & Zacks, S. (1998). Cancer phase I clinical trials: Efficient dose escalation with overdose control.
#' Statistics in medicine 17(10), 1103-1120.
#'
#' Zhou, H.,  Yuan, Y., & Nie, L. (2018). Accuracy, safety, and reliability of novel phase I designs.
#' Clinical Cancer Research, 24(21), 5483-5484 <doi: 10.1158/1078-0432.ccr-18-0168>.
#'
#' @seealso \code{\link[decider:scenario_list_jointBLRM]{scenario_list_jointBLRM}()},
#' \code{\link[decider:sim_jointBLRM]{sim_jointBLRM}()}, \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()},
#' \code{\link[rstan:rstan]{rstan-package}},
#' \code{\link[ggplot2:ggplot]{ggplot2::ggplot}}, \code{\link[ggplot2:ggplot2-package]{ggplot2-package}}.
#'
#'@md
#'@export scenario_jointBLRM

scenario_jointBLRM <- function(data=NULL,
                               historical.data=NULL,
                               doses.of.interest,
                               dose.ref1,
                               dose.ref2,
                               trials.of.interest,
                               types.of.interest=NULL,
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
                               prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
                                                tau_b1 =  c(log(0.125), log(2)/1.96),
                                                tau_a2 =  c(log(0.25),  log(2)/1.96),
                                                tau_b2 =  c(log(0.125), log(2)/1.96),
                                                tau_eta = c(log(0.125), log(2)/1.96)
                               ),
                               saturating = FALSE,
                               probs = c(0.025, 0.5, 0.975),
                               iter = 26000,
                               warmup = 1000,
                               refresh = 0, #floor(iter/10),
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
    if(!all(c("dose1", "dose2", "n.pat", "n.dlt", "trial")%in%names_dat)){

      stop("`data` must be a named list with entries `dose1`, `dose2`, `n.pat`, `n.dlt`, `trial`.")
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
    if(!all(c("dose1", "dose2", "n.pat", "n.dlt", "trial")%in%names_dat)){

      stop("`historical.data` must be a named list with entries `dose1`, `dose2`, `n.pat`, `n.dlt`, `trial`.")
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
      trial = as.array(c(1))
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
    data <- list(
      dose1=as.array(new_dose1),
      dose2=as.array(new_dose2),
      n.pat=as.array(new_n.pat),
      n.dlt=as.array(new_n.dlt),
      trial=as.array(new_trial)
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
        trial = as.array(c(1))
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

  if(!length(dose.ref1) ==1){
    stop("`dose.ref1` must be of length 1")
  }
  if(!length(dose.ref2) ==1){
    stop("`dose.ref2` must be of length 1")
  }
  if((!dose.ref1>0) | (!dose.ref2>0)){

    stop("`dose.ref1` and `dose.ref2` must be positive.")
  }

  if(!length(dose.ref1) ==1){
    stop("`dose.ref1` must be of length 1")
  }
  if(!length(dose.ref2) ==1){
    stop("`dose.ref2` must be of length 1")
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

  res_raw_all <- post_tox_jointBLRM(
    study.interest = std_interest_proc,
    type.interest = types.of.interest,
    dose1.interest=doses.of.interest[1,],
    dose2.interest=doses.of.interest[2,],
    dose1 = data$dose1,
    dose2 = data$dose2,
    dose.ref1 = dose.ref1,
    dose.ref2 = dose.ref2,
    n.pat = data$n.pat,
    n.dlt = data$n.dlt,
    n.study = study_proc,
    MAP.prior = MAP.prior,
    dosing.intervals = dosing.intervals.internal,
    probs = probs,
    prior.mu = prior.mu,
    prior.tau = prior.tau,
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

    if(isMAPcurr==1){
      #str_curr <- paste0("Trial ", std_curr, ": MAP toxicities")
      str_curr <- paste0("trial-", std_curr)
      plotstr_curr <- paste0("trial-", std_curr)
      std_curr_proc <- nstd+1
    }else{
      #str_curr <- paste0("Trial ", std_curr, ": posterior toxicities")
      str_curr <- paste0("trial-", std_curr)
      plotstr_curr <- paste0("trial-", std_curr)
      std_curr_proc <- new_studies[trials.of.interest[std]]
    }

    if(esc.rule=="ewoc"){
      if(!underint.deact){
        summ.current.std <- res_raw[[paste0(std)]]
      }else{
        dim_curr <- dim(res_raw[[paste0(std)]])
        summ.current.std <- matrix(NA, nrow=dim_curr[1], ncol = dim_curr[2]-1)
        rownames(summ.current.std) <- rownames(res_raw[[paste0(std)]])
        colnames(summ.current.std) <- c(colnames(res_raw[[paste0(std)]][, 1:(dim_curr[2]-3)]),
                                        colnames(res_raw[[paste0(std)]][, (dim_curr[2]-1):(dim_curr[2])]))

        summ.current.std[, 1:(dim_curr[2]-3)] <- res_raw[[paste0(std)]][, 1:(dim_curr[2]-3)]
        summ.current.std[, (dim_curr[2]-2):(dim_curr[2]-1)] <- res_raw[[paste0(std)]][, (dim_curr[2]-1):(dim_curr[2])]
      }
    }else{
      #perform loss or dynamic loss escalation
      dims <- dim(res_raw[[paste0(std)]])
      summ.current.std <- matrix(NA, nrow=dims[1], ncol=(dims[2]+1))
      rownames(summ.current.std)<-rownames(res_raw[[paste0(std)]])
      if(esc.rule=="loss"){
        losstype <- "Exp.Loss"
      }else{
        losstype <- "Exp.Loss"
      }
      colnames(summ.current.std)<-c(colnames(res_raw[[paste0(std)]]), losstype)
      summ.current.std[, 1:dims[2]] <- res_raw[[paste0(std)]]
      if(esc.rule=="loss"){
        for(d in 1:dims[1]){
          intprobs_curr <- res_raw[[paste0(std)]][d, (2+length(probs)+1):(2+length(probs)+length(dosing.intervals)+1)]
          summ.current.std[d, dims[2]+1] <- sum(intprobs_curr*loss.weights)
        }
      }else{
        #calculate dynamic loss weights
        ref.p.curr <- ref.probs.all[[paste0(std)]]
        dlwm1 <- ref.p.curr[1, 1]*dynamic.weights[1, ] + ref.p.curr[1, 2]*dynamic.weights[2, ]+
                 ref.p.curr[1, 3]*dynamic.weights[3, ] + ref.p.curr[1, 4]*dynamic.weights[4, ]
        dlwm2 <- ref.p.curr[2, 1]*dynamic.weights[1, ] + ref.p.curr[2, 2]*dynamic.weights[2, ]+
          ref.p.curr[2, 3]*dynamic.weights[3, ] + ref.p.curr[2, 4]*dynamic.weights[4, ]
        dlwc <- ref.p.curr[3, 1]*dynamic.weights[1, ] + ref.p.curr[3, 2]*dynamic.weights[2, ]+
          ref.p.curr[3, 3]*dynamic.weights[3, ] + ref.p.curr[3, 4]*dynamic.weights[4, ]

        doses_all <- rownames(res_raw[[paste0(std)]])

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
          intprobs_curr <- res_raw[[paste0(std)]][d, (2+length(probs)+1):(2+length(probs)+length(dosing.intervals)+1)]
          summ.current.std[d, dims[2]+1] <- sum(intprobs_curr*dlw_curr)
        }
      }

    }

    summary_list[[str_curr]] <- summ.current.std

    if(plot.decisions){
      if(is.null(file.name)){
        file.int.name <- NULL
      }else{
        file.int.name <- paste0(file.name, "_trial-", std_curr)
      }


      if(missing(plot.width)|missing(plot.height)|missing(plot.unit)){

        plot_list[[plotstr_curr]] <- plot_decisions_jointBLRM_int(
          summary=res_raw[[paste0(std)]],
          probs=probs,
          type=type,
          dosing.intervals=dosing.intervals.internal,
          esc.rule=esc.rule,
          ref.probs.all=ref.probs.all[[paste0(std)]],
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
                                       summary=res_raw[[paste0(std)]],
                                       probs=probs,
                                       type=type,
                                       dosing.intervals=dosing.intervals.internal,
                                       esc.rule=esc.rule,
                                       ref.probs.all=ref.probs.all[[paste0(std)]],
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

    }

  }

  if(output.scen.config){
    summary_list[["data"]] <- data_matrix_jointBLRM(data)
    summary_list[["prior"]] <- prior_mat_out(m=prior.mu, t=prior.tau)
    if(esc.rule%in%c("ewoc")){
      outconf <- matrix("-", nrow=3, ncol=1)
      rownames(outconf) <- c("seed","saturating", "ewoc.threshold")
      colnames(outconf) <- c("-")
      outconf[1, 1] <- seed
      outconf[2, 1] <- saturating
      outconf[3, 1] <- ewoc.threshold
      summary_list[["configuration"]] <- outconf
    }else if(esc.rule%in%c("loss")){
      outconf <- matrix("-", nrow=3, ncol=4)
      rownames(outconf) <- c("seed","saturating", "loss.weights")
      colnames(outconf) <- c("-", "-", "-", "-")
      outconf[1, 1] <- seed
      outconf[2, 1] <- saturating
      outconf[3, ] <- loss.weights
      summary_list[["configuration"]] <- outconf
    }else{
      outconf <- matrix("-", nrow=6, ncol=4)
      rownames(outconf) <- c("seed","saturating", "dynamic.weights[1,]", "dynamic.weights[2,]"
                             , "dynamic.weights[3,]", "dynamic.weights[4,]")
      colnames(outconf) <- c("-", "-", "-", "-")
      outconf[1, 1] <- seed
      outconf[2, 1] <- saturating
      outconf[3, ] <- dynamic.weights[1,]
      outconf[4, ] <- dynamic.weights[2,]
      outconf[5, ] <- dynamic.weights[3,]
      outconf[6, ] <- dynamic.weights[4,]
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


#'Process a list of data scenarios in parallel.
#'@rdname scenario_list_jointBLRM
#'@description Wrapper for calling \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()}
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
#' @seealso \code{\link[decider:scenario_jointBLRM]{scenario_jointBLRM}()}, \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}()},
#' \code{\link[rstan:rstan]{rstan-package}},
#' \code{\link[ggplot2:ggplot]{ggplot2::ggplot}}, \code{\link[ggplot2:ggplot2-package]{ggplot2-package}}.
#' @export
scenario_list_jointBLRM <- function(
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
                      .export = c("stanmodels"),
                      .errorhandling = "pass",
                      .inorder = TRUE)%dopar%
    {

      return(scenario_jointBLRM(
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
