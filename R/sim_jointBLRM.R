# --------------------------
# function: sim_jointBLRM
# --------------------------
#'Simulate dose-finding trials
#'@description Simulates dose-finding trials with up to six parallel monotherapy or two-drug combination therapy trials modeled together in a
#'joint BLRM. The function assumes that two different compounds are involved, compounds 1 and 2 (and their combination).
#'Up to two monotherapy trials for each compound can be actively simulated, and additionally two combination therapy trials.
#'Note that the term "trial" can in this context also refer to trial arms (e.g. monotherapy and combination therapy), and
#'that the underlying model will not differentiate between these notions.
#'The order in which cohorts for the simulated trials are enrolled can be specified freely. Further, besides the actively simulated
#'trials, the function supports the inclusion of historical data from arbitrarily many previous (already concluded)
#'dose-finding trials for the considered compounds.
#'
#'The function expects a specification of dose levels and their assumed DLT rates for each involved trial, based on which the course of
#'the dose-finding trial is simulated. Additional customization of the simulations is possible (separately for each simulated trial),
#'e.g. regarding the escalation rule (EWOC and loss-based methods), definition of MTD, the maximum increment of escalations,
#'the numbers of patients per cohort, and further specifications.
#'
#'The six possible trials are encoded with the following suffixes: \code{mono1.a}, \code{mono1.b}, \code{mono2.a},
#'\code{mono2.b}, \code{combi.a}, and \code{combi.b}. The ending (\code{.a} or \code{.b}) differentiates between two different
#'trials of the same kind, i.e., monotherapy for compound 1 (\code{mono1}), monotherapy for compound 2 (\code{mono2}), and combination therapy (\code{combi}).
#'For some function arguments, one version of the argument is available for each trial (differentiated via their suffixes),
#'so that different simulated trials may use different specifications in terms of e.g. dose levels, assumed DLT rate,
#'and starting dose.
#'
#'To differentiate between the up to six simulated trials and the actual simulation run that processes one realization
#'of the activated trials, we will in the following use the term "study" for a single simulation run that
#'determines a realization of each activated trial using the given specifications. For instance, if two
#'trials are active, say \code{mono1.a} and \code{combi.a}, the term "trial" refers to the simulated trials
#'with these names, while the notion "study" means a simulation run of both trials in the specified order
#'and under the given input assumptions. Simulating e.g. 1000 studies means in this context that the trials
#'\code{mono1.a} and \code{combi.a} are simulated in the specified order using the specified options for 1000 times.
#'
#'Most of the parameters have sensible defaults and will often not need to be specified explicitly.
#'The recommended workflow to set up the function is as follows:
#'\itemize{
#'\item{Specify reference doses.}
#'\item{Optional: Specify available historical data and, if needed, specific priors for the model (the default ones are weakly informative).}
#'\item{Specify dose levels, starting doses and assumed DLT rates for each trial that needs to be simulated, and
#'activate the corresponding trials via the \code{active.[...]} arguments.}
#'\item{Set the number of trials to be simulated.}
#'\item{Optional: Specify additional options for the simulations.}
#'\item{Optional but strongly recommended: Enable parallel simulation and specify the number of cores, and supply a path for a working directory
#'for saving and re-loading MCMC results. When a working directory is supplied, be aware that the function will delete or modify
#'existing .RData files containing "_tmp" and the supplied \code{file.name} in their name from the working directory to avoid
#'re-loading MCMC runs from previous simulations.
#'Refer to the documentation below for more detail.}
#'\item{Run the simulation. Note that simulations of e.g. 1000 trials may take multiple hours depending on the specifications.}
#'}
#'
#'@usage
#'sim_jointBLRM(
#'    active.mono1.a = FALSE,
#'    active.mono1.b = FALSE,
#'    active.mono2.a = FALSE,
#'    active.mono2.b = FALSE,
#'    active.combi.a = FALSE,
#'    active.combi.b = FALSE,
#'    doses.mono1.a,
#'    doses.mono2.a,
#'    doses.mono1.b,
#'    doses.mono2.b,
#'    doses.combi.a,
#'    doses.combi.b,
#'    dose.ref1,
#'    dose.ref2,
#'    tox.mono1.a,
#'    tox.mono2.a,
#'    tox.combi.a,
#'    tox.mono1.b,
#'    tox.mono2.b,
#'    tox.combi.b,
#'    start.dose.mono1.a,
#'    start.dose.mono2.a,
#'    start.dose.mono1.b,
#'    start.dose.mono2.b,
#'    start.dose.combi.a1,
#'    start.dose.combi.a2,
#'    start.dose.combi.b1,
#'    start.dose.combi.b2,
#'    cohort.queue = rep( c(1,2,3,4,5,6), times= 100),
#'    historical.data = NULL,
#'    esc.rule = "ewoc",
#'    esc.comp.max=1,
#'    dosing.intervals = c(0.16, 0.33, 0.6),
#'    ewoc.threshold = 0.25,
#'    loss.weights = c(1, 0, 1, 2),
#'    dynamic.weights =  rbind(
#'        c(0.32, 0, 0.32, 0.36),
#'        c(0.29, 0, 0.31, 0.40),
#'        c(0.27, 0, 0.33, 0.40),
#'        c(0.20, 0, 0.30, 0.50)),
#'    prior.mu = list(
#'        mu_a1 = c(logit(0.33), 2),
#'        mu_b1 = c(0, 1),
#'        mu_a2 = c(logit(0.33), 2),
#'        mu_b2 = c(0, 1),
#'        mu_eta = c(0, 1.121)),
#'    prior.tau = list(
#'        tau_a1 = c(log(0.25), log(4)/1.96),
#'        tau_b1 = c(log(0.125), log(4)/1.96),
#'        tau_a2 = c(log(0.25), log(4)/1.96),
#'        tau_b2 = c(log(0.125), log(4)/1.96),
#'        tau_eta = c(log(0.125),log(4)/1.96)),
#'    saturating = FALSE,
#'    esc.step = NULL,
#'    esc.step.mono1.a = esc.step,
#'    esc.step.mono2.a = esc.step,
#'    esc.step.mono1.b = esc.step,
#'    esc.step.mono2.b = esc.step,
#'    esc.step.combi.a1 = esc.step,
#'    esc.step.combi.b1 = esc.step,
#'    esc.step.combi.a2 = esc.step,
#'    esc.step.combi.b2 = esc.step,
#'    esc.constrain = FALSE,
#'    esc.constrain.mono1.a=esc.constrain,
#'    esc.constrain.mono2.a=esc.constrain,
#'    esc.constrain.mono1.b=esc.constrain,
#'    esc.constrain.mono2.b=esc.constrain,
#'    esc.constrain.combi.a1=esc.constrain,
#'    esc.constrain.combi.b1=esc.constrain,
#'    esc.constrain.combi.a2=esc.constrain,
#'    esc.constrain.combi.b2=esc.constrain,
#'    cohort.size = c(3),
#'    cohort.size.mono1.a = cohort.size,
#'    cohort.size.mono1.b = cohort.size,
#'    cohort.size.mono2.a = cohort.size,
#'    cohort.size.mono2.b = cohort.size,
#'    cohort.size.combi.a = cohort.size,
#'    cohort.size.combi.b = cohort.size,
#'    cohort.prob = NULL,
#'    cohort.prob.mono1.a = cohort.prob,
#'    cohort.prob.mono1.b = cohort.prob,
#'    cohort.prob.mono2.a = cohort.prob,
#'    cohort.prob.mono2.b = cohort.prob,
#'    cohort.prob.combi.a = cohort.prob,
#'    cohort.prob.combi.b = cohort.prob,
#'    max.n = 42,
#'    max.n.mono1.a = max.n,
#'    max.n.mono1.b = max.n,
#'    max.n.mono2.a = max.n,
#'    max.n.mono2.b = max.n,
#'    max.n.combi.a = max.n,
#'    max.n.combi.b = max.n,
#'    mtd.decision = list(
#'        target.prob = 0.5,
#'        pat.at.mtd = 6,
#'        min.pat = 12,
#'        min.dlt = 1,
#'        rule = 2),
#'    mtd.decision.combi.a = mtd.decision,
#'    mtd.decision.combi.b = mtd.decision,
#'    mtd.decision.mono1.a = mtd.decision,
#'    mtd.decision.mono1.b = mtd.decision,
#'    mtd.decision.mono2.a = mtd.decision,
#'    mtd.decision.mono2.b = mtd.decision,
#'    mtd.enforce = FALSE,
#'    mtd.enforce.mono1.a = mtd.enforce,
#'    mtd.enforce.mono2.a = mtd.enforce,
#'    mtd.enforce.mono1.b = mtd.enforce,
#'    mtd.enforce.mono2.b = mtd.enforce,
#'    mtd.enforce.combi.a = mtd.enforce,
#'    mtd.enforce.combi.b = mtd.enforce,
#'    backfill.mono1.a = FALSE,
#'    backfill.mono1.b = FALSE,
#'    backfill.mono2.a = FALSE,
#'    backfill.mono2.b = FALSE,
#'    backfill.combi.a = FALSE,
#'    backfill.combi.b = FALSE,
#'    backfill.size = c(3),
#'    backfill.prob = NULL,
#'    backfill.size.mono1.a = backfill.size,
#'    backfill.size.mono1.b = backfill.size,
#'    backfill.size.mono2.a = backfill.size,
#'    backfill.size.mono2.b = backfill.size,
#'    backfill.size.combi.a = backfill.size,
#'    backfill.size.combi.b = backfill.size,
#'    backfill.prob.mono1.a = backfill.prob,
#'    backfill.prob.mono1.b = backfill.prob,
#'    backfill.prob.mono2.a = backfill.prob,
#'    backfill.prob.mono2.b = backfill.prob,
#'    backfill.prob.combi.a = backfill.prob,
#'    backfill.prob.combi.b = backfill.prob,
#'    backfill.start.mono1.a = NULL,
#'    backfill.start.mono1.b = NULL,
#'    backfill.start.mono2.a = NULL,
#'    backfill.start.mono2.b = NULL,
#'    backfill.start.combi.a1 = NULL,
#'    backfill.start.combi.a2 = NULL,
#'    backfill.start.combi.b1 = NULL,
#'    backfill.start.combi.b2 = NULL,
#'    n.studies = 1,
#'    n.cores = 1,
#'    seed = sample.int(.Machine$integer.max, 1),
#'    chains = 4,
#'    iter = 6000,
#'    warmup = 1000,
#'    adapt_delta = 0.9,
#'    max_treedepth = 12,
#'    refresh=0,
#'    file.name = NULL,
#'    path = NULL,
#'    monitor.path = NULL,
#'    working.path = NULL,
#'    clean.working.path = FALSE,
#'    output.sim.config =TRUE
#')
#'@param active.mono1.a,active.mono1.b,active.mono2.a,active.mono2.b,active.combi.a,active.combi.b Logicals, default to \code{FALSE}. Each parameter activates the simulation of the
#'corresponding trial. If the parameter is \code{FALSE}, the corresponding trial is not simulated, all parameters with this suffix are ignored, and no historical information can be included for this trial.
#'@param doses.mono1.a,doses.mono1.b,doses.mono2.a,doses.mono2.b Numericals, one dimensional vectors with positive, strictly ascending entries. Each vector indicates the dose levels for the corresponding monotherapy trial.
#'Optional if \code{active.[...]} is \code{FALSE} for the corresponding suffix.
#'@param doses.combi.a,doses.combi.b Numericals, two dimensional arrays with two rows, arbitrarily many columns, and positive entries. Does not need to be ascending. Each column indicates a combination dose level for the corresponding combination therapy trial,
#'where the entry in the first row of a specific column gives the dose level of compound 1, and the entry in the second row the dose level of compound 2. Optional if \code{active.combi.[...]} is \code{FALSE} for the corresponding suffix.
#'@param dose.ref1 Positive numerical value. Reference dose of compound 1.
#'@param dose.ref2 Positive numerical value. Reference dose of compound 2.
#'@param tox.mono1.a,tox.mono1.b,tox.mono2.a,tox.mono2.b Numericals, one dimensional vectors with entries in \eqn{(0,1)}.  Optional if \code{active.mono.[...]} is \code{FALSE} for the corresponding suffix. Otherwise,
#'must have the same length as \code{doses.[...]} for the corresponding suffix. The entry \eqn{i} specifies the DLT rate of dose \eqn{i} of the corresponding trial.
#'@param tox.combi.a,tox.combi.b Numericals, one dimensional vectors with entries in \eqn{(0,1)}.  Optional if \code{active.combi.[...]} is \code{FALSE} for the corresponding suffix. Otherwise,
#'its length must equal the number of columns of \code{doses.combi.[...]}. The entry \eqn{i} specifies the DLT rate of the dose combination of row \eqn{i} of \code{doses.combi.[...]}.
#'@param start.dose.mono1.a,start.dose.mono2.a,start.dose.mono1.b,start.dose.mono2.b Positive numerical values that specify the starting dose for the simulated monotherapy trials.
#'The value of \code{start.dose.mono.[...]} must be one of the entries of \code{doses.mono.[...]} for the corresponding suffix. The starting doses are optional if the corresponding trial is not activated.
#'@param start.dose.combi.a1,start.dose.combi.a2,start.dose.combi.b1,start.dose.combi.b2 Positive numerical values that specify the starting dose for the simulated combination therapy trials.
#'The values with suffixes \code{.a1} and  \code{.a2} define the starting dose for first and second compound of the \code{combi.a} trial, and must appear as combination in one of the columns of \code{doses.combi.a}. The same applies for the trial with suffix \code{combi.b}. The starting doses are optional if the corresponding trial is not activated.
#'@param cohort.queue Optional numerical or character vector that specifies the order or pattern in which cohorts are enrolled in the simulated trials. By default, cohorts are enrolled sequentially in the order:
#'\code{mono1.a}, \code{mono2.a}, \code{combi.a}, \code{mono1.b}, \code{mono2.b}, \code{combi.b}, which is repeated until all trials concluded.
#' In general, entry \eqn{i} specifies to which simulated trial the \eqn{i}th simulated cohort is enrolled (provided that this trial has not concluded yet, in which case the entry would be ignored).
#'The \code{cohort.queue} can therefore be used to express e.g. that a certain number of monotherapy cohorts are completed prior to the start of a combination therapy trial.
#'The \code{cohort.queue} can either be a short pattern that is repeated automatically, or alternatively a long vector
#'that contains already the order of all cohorts that may occur in the trial. The latter can e.g. be used to ensure that a trial
#'has concluded prior to the start of some other trial by supplying enough cohorts to the former to reach its maximum sample size before the first
#'cohort for the latter trial is enrolled.
#'
#'The trials to which cohorts are assigned are encoded as follows:
#'\itemize{
#'\item{\code{"mono1.a"} or \code{1}\cr Corresponds to a cohort for the \code{mono1.a} trial (monotherapy trial for compound 1).}
#'\item{\code{"mono2.a"} or \code{2}\cr Corresponds to a cohort for the \code{mono2.a} trial (monotherapy trial for compound 2).}
#'\item{\code{"combi.a"} or \code{3}\cr Corresponds to a cohort for the \code{combi.a} trial (combination therapy trial).}
#'\item{\code{"mono1.b"} or \code{4}\cr Corresponds to a cohort for the \code{mono1.b} trial (monotherapy trial for compound 1).}
#'\item{\code{"mono2.b"} or \code{5}\cr Corresponds to a cohort for the \code{mono2.b} trial (monotherapy trial for compound 2).}
#'\item{\code{"combi.b"} or \code{6}\cr Corresponds to a cohort for the \code{combi.b} trial (combination therapy trial).}
#'}
#'Assuming that each cohort is assigned to the minimal possible cohort size (specified via the argument \code{cohort.size}),
#'the \code{cohort.queue} should contain sufficiently many entries for each activated trial so that it can in principle reach the
#'specified maximum patient number to be enrolled in the trial. For instance, if cohorts consist of 3 or more patients and the maximum sample size
#'(given in \code{max.n.[...]}) of this trial is 30, than the trial should appear at least 10 times in the cohort queue. If this
#'is not the case, the function will repeat the cohort queue until sufficiently many cohorts are contained.
#'@param historical.data Optional parameter that must be \code{NULL} (the default) or a named list that specifies the historical data to be included in the trial.
#'Historical data can be included from both purely historical trials (i.e., not actively simulated) and also directly to one of the
#'simulated trials (in which case the simulation assumes that this data was previously recorded in the corresponding trial).
#'
#'A value of \code{NULL} represents that no historical data is available. Otherwise, the historical data is specified via the named list entries \code{dose1}, \code{dose2}, \code{n.pat}, \code{n.dlt}, and \code{trial}
#'(where the names are not case-sensitive). The list entries should be vectors whose length is equal to the number of historical cohorts to be included.
#'Their format should be as follows:
#'\itemize{
#'\item{\code{historical.data$dose1}\cr Numerical vector with non-negative or \code{NA} entries. Entry \eqn{i} gives the dose level of compound 1 that was administered to the \eqn{i}th
#'historical cohort. If compound 1 was not administered to this cohort, the entry should be \code{0} or \code{NA}.}
#'\item{\code{historical.data$dose2}\cr Numerical vector with non-negative or \code{NA} entries. Entry \eqn{i} gives the dose level of compound 2 that was administered to the \eqn{i}th
#'historical cohort. If compound 2 was not administered to this cohort, the entry should be \code{0} or \code{NA}.}
#'\item{\code{historical.data$n.pat}\cr Numerical vector with non-negative entries. Entry \eqn{i} gives the number of patients of the \eqn{i}th
#'historical cohort. If the patient number is \code{0}, the historical cohort is removed automatically. \code{NA} entries are not permitted.}
#'\item{\code{historical.data$n.dlt}\cr Numerical vector with non-negative entries. Entry \eqn{i} gives the number of patients that experienced DLT among the \eqn{i}th
#'historical cohort. The number of DLTs is not allowed be larger than the number of patients in the cohort. \code{NA} entries are not permitted.}
#'\item{\code{historical.data$trial}\cr Numerical or character vector with arbitrary entries (except for \code{NA}
#'which is not permitted). Entry \eqn{i} gives the code for the trial from which the \eqn{i}th
#'historical cohort arises. That is, historical patients are grouped in different historical or concurrent trials according to these values,
#'where the entries are not case sensitive and do not differentiate between numbers and the converted strings (i.e., \code{42} and \code{"42"} represent the same trial).
#'The following entries are reserved to represent that the cohort belongs to one of the actively simulated trials:
#'\itemize{
#'\item{\code{"mono1.a"} or \code{1}\cr The cohort is assumed to belong to the actively simulated trial \code{mono1.a}.}
#'\item{\code{"mono2.a"} or \code{2}\cr The cohort is assumed to belong to the actively simulated trial \code{mono2.a}.}
#'\item{\code{"combi.a"} or \code{3}\cr The cohort is assumed to belong to the actively simulated trial \code{combi.a}.}
#'\item{\code{"mono1.b"} or \code{4}\cr The cohort is assumed to belong to the actively simulated trial \code{mono1.b}.}
#'\item{\code{"mono2.b"} or \code{5}\cr The cohort is assumed to belong to the actively simulated trial \code{mono2.b}.}
#'\item{\code{"combi.b"} or \code{6}\cr The cohort is assumed to belong to the actively simulated trial \code{combi.b}.}
#'}
#'
#'All other strings or numbers are interpreted to belong to purely historical (already completed) trials.
#'If data for the actively simulated trials are included, the corresponding trial must be activated (otherwise, the reserved
#'trial codes are not permitted), and the historical cohorts must match the type of the trial. The latter means that historical
#'data for e.g. the \code{mono1.a} trial must contain observations from monotherapy with compound 1.}
#'}
#'@param esc.rule Optional character string, must have one of the following values: \code{"ewoc"}, \code{"ewoc.opt"}, \code{"ewoc.max"},
#'\code{"loss"}, \code{"dynamic"}, or \code{"dynamic.loss"}, where the default is \code{"ewoc"}. The value of \code{esc.rule}
#'specifies which escalation rule is applied during the simulation, cf. Details. The character strings are
#'interpreted in a case-insensitive fashion. The possible values are interpreted as follows:
#'\itemize{
#'\item{\code{"ewoc"} or \code{"ewoc.opt"}\cr Applies the EWOC criterion with optimal probability rule and three dosing intervals (underdosing, target dosing, and overdosing). That is, the next dose is
#'the one with the largest probability that the true DLT rate lies in the target interval
#'among the doses whose probability of having a DLT rate in the overdosing interval is lower than the feasibility bound specified in the argument \code{ewoc.threshold}. If multiple doses or
#'dose combinations attain the maximum target probability and satisfy the EWOC principle, the escalation rule will additionally maximize the dose level of the compound that is
#'specified in the argument \code{esc.comp.max}.}
#'\item{\code{"ewoc.max"}\cr Applies the EWOC criterion with maximum dose rule and three dosing intervals (underdosing, target dosing, and overdosing). That is, the next dose is the largest one
#'among the doses whose probability of having true DLT rate in the overdosing interval is lower than the feasibility bound specified in the argument \code{ewoc.threshold}. For dose combinations,
#'the escalation rule will first maximize the compound specified in \code{esc.comp.max}, and afterwards the other compound.}
#'\item{\code{"loss"}\cr Applies the static loss escalation rule using the weights specified in the argument \code{loss.weights} and four dosing intervals
#'(under dosing, target dosing, excessive toxicity, and unacceptable toxicity). When multiple dose combinations
#'result in the same expected loss (Bayes risk), the function will additionally maximize the compound specified in the argument \code{esc.comp.max}
#'among the doses with minimum expected loss.}
#'\item{\code{"dynamic"} or \code{"dynamic.loss"}\cr Applies the dynamic loss escalation rule using the weights specified in the argument \code{dynamic.weights} and four dosing intervals
#'(under dosing, target dosing, excessive toxicity, and unacceptable toxicity). When multiple dose combinations
#'result in the same expected loss (Bayes risk), the function will additionally maximize the compound specified in the argument \code{esc.comp.max}
#'among the doses with minimum expected loss.}
#'}
#'@param esc.comp.max Optional integer, must be either \code{1} (the default) or \code{2}. The value specifies which compound (1 or 2) is maximized first when an escalation rule
#'leads to a draw among multiple dose combinations.
#'@param dosing.intervals Optional numeric vector with ascending entries between 0 and 1. The argument \code{dosing.intervals} must be of length 1, 2, or 3, and defines the boundaries
#'of the dosing intervals used for escalation decisions. The default is \code{c(0.16, 0.33, 0.6)} which can be used for any escalation rule. The interpretation and requirements vary slightly across different escalation rules (different values of \code{esc.rule}).
#'
#'More precisely, the following cases are differentiated. If \code{esc.rule} is...
#'\itemize{
#'\item{\code{"ewoc"},\code{"ewoc.max"}, or \code{"ewoc.opt"}\cr The length can be either 1, 2, or 3. In the latter two cases,
#'only the first two entries are used to specify the targeted dosing interval. That is, the target interval is
#'assumed to range from \code{dosing.intervals[1]} to \code{dosing.intervals[2]}, while a potential third entry is discarded.
#'If \code{dosing.intervals} has length 1 and value \code{x}, the argument is internally transformed to \code{c(0, x)}. That is,
#'the target interval is assumed to range from \code{0} to \code{x}, while the underdosing interval is not considered.\cr
#'The overdosing interval is implicitly defined to range from \code{dosing.intervals[2]} to \code{1}, while the underdosing
#'interval is defined to range from \code{0} to \code{dosing.intervals[1]}.}
#'\item{\code{"loss"},\code{"dynamic.loss"}, or \code{"dynamic"}\cr Must have three entries which are used to specify the targeted dosing interval and the interval for excessively toxic dosing.
#'That is, the target interval is assumed to range from \code{dosing.intervals[1]} to \code{dosing.intervals[2]}, while the interval
#'for excessively toxic dosing is assumed to range from \code{dosing.intervals[2]} to \code{dosing.intervals[3]}. The underdosing interval is
#'therefore implicitly defined to range from  \code{0} to \code{dosing.intervals[1]}, while the unacceptably toxic dosing interval is implicitly defined
#'to range from \code{dosing.intervals[3]} to \code{1}.}
#'}
#'@param ewoc.threshold Optional numerical value between 0 and 1 (excluding the boundaries), defaults to \code{0.25}. Defines the feasibility bound for the EWOC criterion,
#'i.e., when \code{esc.rule} is one of the EWOC-based rules, the next dose must have a probability of having DLT rate in the overdosing interval lower than \code{ewoc.threshold}.
#'The value is ignored when \code{esc.rule} specifies one of the loss-based escalation rules.
#'@param loss.weights Optional numerical vector with four entries (which can be arbitrary numbers), the default is \code{c(1,0,1,2)}. Specifies the interval weights/penalties
#'that are used for static loss escalation. This is only needed when \code{esc.rule} is \code{"loss"}, and the parameter is ignored otherwise.
#'More precisely, \code{loss.weights[1]} is the weight of the underdosing interval, \code{loss.weights[2]} the weight of the  target dosing interval,
#'\code{loss.weights[3]} the weight of the excessively toxic dosing interval, and \code{loss.weights[4]} is the weight of the unacceptably toxic dosing interval.
#'See the Details section below for explanations regarding the methodology.
#'@param dynamic.weights Optional numerical matrix with four rows and four columns, and arbitrary numbers as entries.  Specifies the interval weights/penalties
#'that are used for dynamic loss escalation. This is only needed when \code{esc.rule} is \code{"dynamic"} or \code{"dynamic.loss"}, and the parameter is ignored otherwise.
#'Each row of \code{dynamic.weights} specifies one of the static loss weight vectors that are interpolated during dynamic loss escalation. See the Details section below for explanations regarding the methodology.
#'
#'More precisely, \code{dynamic.weights[1, ]} is the static weight vector that is weighted by the posterior probability that the reference dose has true DLT rate in the underdose interval,
#'\code{dynamic.weights[2, ]} is the static weight vector that is weighted by the posterior probability that the reference dose has true DLT rate in the target interval,
#'\code{dynamic.weights[3, ]} is the static weight vector that is weighted by the posterior probability that the reference dose has true DLT rate in the excessively toxic dosing interval, and
#'\code{dynamic.weights[4, ]} is the static weight vector that is weighted by the posterior probability that the reference dose has true DLT rate in the unacceptably toxic dosing interval.
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
#'\item{\code{prior.mu$mu_a1}\cr Numeric with length two, defaults to \code{c(logit(0.33), 2)}. Specifies mean and SD of the hypermean \eqn{\mu_1} of the parameter \eqn{log(\alpha_1)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.mu$mu_b1}\cr Numeric with length two, defaults to \code{c(0, 1)}. Specifies mean and SD of the hypermean \eqn{\mu_2} of the parameter \eqn{log(\beta_1)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.mu$mu_a2}\cr Numeric with length two, defaults to \code{c(logit(0.33), 2)}. Specifies mean and SD of the hypermean \eqn{\mu_3} of the parameter \eqn{log(\alpha_2)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.mu$mu_b2}\cr Numeric with length two, defaults to \code{c(0, 1)}. Specifies mean and SD of the hypermean \eqn{\mu_4} of the parameter \eqn{log(\beta_2)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.mu$mu_eta}\cr Numeric with length two, defaults to \code{c(0, 1.121)}. Specifies mean and SD of the hypermean \eqn{\mu_5} of the parameter \eqn{\eta} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'}
#'@param prior.tau Optional list that gives the prior distribution for the between-trial heterogeneities (hyper SD) \eqn{\tau}. The list must
#'have named entries, which all need to be numeric vectors of length 2:
#'\itemize{
#'\item{\code{prior.tau$tau_a1}\cr Numeric with length two, defaults to \code{c(log(0.25), log(2)/1.96)}. Specifies mean and SD on log-scale of the heterogeneity \eqn{\tau_1} of the parameter \eqn{log(\alpha_1)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.tau$tau_b1}\cr Numeric with length two, defaults to \code{c(log(0.125), log(2)/1.96)}. Specifies mean and SD on log-scale of the heterogeneity \eqn{\tau_2} of the parameter \eqn{log(\beta_1)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.tau$tau_a2}\cr Numeric with length two, defaults to \code{c(log(0.25), log(2)/1.96)}. Specifies mean and SD on log-scale of the heterogeneity \eqn{\tau_3} of the parameter \eqn{log(\alpha_2)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.tau$tau_b2}\cr Numeric with length two, defaults to \code{c(log(0.125), log(2)/1.96)}. Specifies mean and SD on log-scale of the heterogeneity \eqn{\tau_4} of the parameter \eqn{log(\beta_2)} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'\item{\code{prior.tau$tau_eta}\cr Numeric with length two, defaults to \code{c(log(0.125), log(2)/1.96)}. Specifies mean and SD on log-scale of the heterogeneity \eqn{\tau_5} of the parameter \eqn{\eta} in the BLRM. The second entry must
#'therefore be positive. See the Details section below for more detail.}
#'}
#'@param saturating Optional logical, defaults to \code{FALSE}. If \code{TRUE}, the BLRM will be using a saturating interaction term as described in
#'\code{\link[OncoBayes2:blrm_formula_saturating]{OncoBayes2::blrm_formula_saturating}()}. Also refer to the Details section in the documentation
#'of \code{\link[OncoBLRM:scenario_jointBLRM]{scenario_jointBLRM}()}.
#'@param esc.step,esc.step.mono1.a,esc.step.mono2.a,esc.step.combi.a1,esc.step.combi.a2 Optional numerical values that specify
#'the maximum factor of dose escalations that is demanded additionally to the selected escalation rule. The default is \code{NULL}, in which case the escalation step
#'is determined automatically for the corresponding trial (if activated) as the maximum factor between the (sorted) available doses for the trial.
#'As an example, an value of \code{esc.step=2} indicates that the dose can at most be doubled, or, in other words, that the increment between subsequent dose
#'levels is at most 100%.
#'
#'The argument \code{esc.step} is a global variant, all escalation steps that were not specified by the user default to its value. In particular, if all trials shall use the same escalation factor,
#'only the value of \code{esc.step} needs to be defined, while the trial-specific variants can be used to express that a specific trial shall deviate from
#'the global parameter. Note that consistency checks are performed for these parameters, i.e., whether the specified escalation step allows to reach all available doses.
#'
#'The value of \code{esc.step.[...]} gives the escalation step that is used for the corresponding trial. Note that there are two escalation steps for combination therapy trials:
#'For instance, \code{esc.step.combi.a1} gives the maximum escalation factor in the first compound during combination therapy escalations, while
#'\code{esc.step.combi.a2} defines the escalation step for the second compound.
#'@param esc.step.mono1.b,esc.step.mono2.b,esc.step.combi.b1,esc.step.combi.b2
#'Same as \code{esc.step.[...].a} (where \code{[...]} is \code{mono1}, \code{mono2}, or \code{combi})
#'but for the second potentially simulated trial (suffix \code{.b}) of the respective trial type.
#'@param esc.constrain,esc.constrain.mono1.a,esc.constrain.mono2.a
#'Optional logicals, default to \code{FALSE}. When \code{TRUE}, escalation decisions in the corresponding trial
#'are constrained to the available dose levels, in the sense that only the next larger dose can be used for
#'the next cohort. This does not apply to de-escalation, where arbitrary lower dose levels are allowed.
#'Note that \code{esc.constrain[...]} overrules \code{esc.step[...]}, so that the given
#'escalation step is ignored when \code{esc.constrain[...]} is \code{TRUE}.
#'
#'As before, \code{esc.contrain} is the global variant and affects all trials. If a specific trial
#'shall deviate from the global setting, this needs to be activated specifically with
#'the trial specific variant \code{esc.constrain.[...]} using the corresponding suffix.
#'@param esc.constrain.combi.a1,esc.constrain.combi.a2
#'Optional logicals, default to \code{FALSE}. Work in the same way as \code{esc.constrain.mono[...]},
#'except that the compounds in a combination therapy trial are constrained separately from each other, where
#'e.g. \code{esc.constrain.combi.a1} affects compound 1, while  \code{esc.constrain.combi.a2} affects
#'compound 2. For instance, it is possible to activate the restriction to the next larger dose level only
#'for the first compound, and still use an escalation step for the second compound.
#'@param esc.constrain.mono1.b,esc.constrain.mono2.b,esc.constrain.combi.b1,esc.constrain.combi.b2
#'Same as \code{esc.constrain.[...].a} (where \code{[...]} is \code{mono1}, \code{mono2}, or \code{combi})
#'but for the second potentially simulated trial (suffix \code{.b}) of the respective trial type.
#'@param cohort.size,cohort.size.mono1.a,cohort.size.mono2.a,cohort.size.combi.a
#'Optional positive integer vectors that specify the available cohort sizes for simulated cohorts. The global variant, \code{cohort.size}, defaults to \code{3} (i.e., all cohorts will contain 3 patients),
#'the remaining trial-specific parameters default to the value of \code{cohort.size}. In particular, if all trials shall use the same cohort sizes,
#'only the value of \code{cohort.size} needs to be defined, while the trial-specific variants can be used to express that a specific trial shall deviate from
#'the global parameter.
#'
#'The \code{cohort.size} can either be a single value (all cohorts use the same size) or a vector of integers that give the available cohort sizes. In the latter case,
#'the simulated cohort sizes are drawn at random from the possible cohort sizes using the probability for each cohort size that is specified in \code{cohort.prob}.
#'As before, the trial-specific variant, \code{cohort.size.[...]}, only affect the corresponding trial.
#'@param cohort.size.mono1.b,cohort.size.mono2.b,cohort.size.combi.b
#'Same as \code{cohort.size.[...].a} (where \code{[...]} is \code{mono1}, \code{mono2}, or \code{combi})
#'but for the second potentially simulated trial (suffix \code{.b}) of the respective trial type.
#'@param cohort.prob,cohort.prob.mono1.a,cohort.prob.mono2.a,cohort.prob.combi.a
#'Optional positive numeric vectors with values between 0 and 1 that specify the probability for each of the available cohort sizes in the corresponding argument \code{cohort.size.[...]}. The global variant, \code{cohort.prob}, defaults to \code{NULL},
#'the remaining trial-specific parameters default to the value of \code{cohort.prob}. If the parameter is \code{NULL} for a trial,
#'the function will assign the same probability to each entry of \code{cohort.size}, i.e., the cohort size is chosen uniformly at random from
#'the provided possible cohort sizes.
#'If all trials shall use the same probabilities for selecting the cohort size,
#'only the value of \code{cohort.prob} needs to be defined, while the trial-specific variants can be used to express that a specific trial shall deviate from
#'the global parameter.
#'
#'Each of the arguments \code{cohort.prob.[...]} must be \code{NULL} or a vector of the same length as the argument \code{cohort.size.[...]}, where the lengths and contents may
#'differ across different trials/suffixes.
#'@param cohort.prob.mono1.b,cohort.prob.mono2.b,cohort.prob.combi.b
#'Same as \code{cohort.prob.[...].a} (where \code{[...]} is \code{mono1}, \code{mono2}, or \code{combi})
#'but for the second potentially simulated trial (suffix \code{.b}) of the respective trial type.
#'@param max.n,max.n.mono1.a,max.n.mono2.a,max.n.combi.a Optional positive integer values that specify
#'the maximum number of patients to be enrolled in simulated trials. If the simulated number of patients in a trial equals or surpasses the specified number in \code{max.n.[...]} and
#'has not stopped prematurely and not determined an MTD yet, the simulated trial is stopped without declaring an MTD due to reaching the maximum sample size.
#'
#'The global variant of the parameter, \code{max.n}, defaults to 42 (i.e., MTD must be reached after 14 cohorts when the default values for \code{cohort.size} and \code{cohort.prob} are used).
#'The trial-specific parameters, \code{max.n.[...]}, default to the global value \code{max.n} and only need to be adjusted if the corresponding trial shall deviate from the
#'global value.
#'@param max.n.mono1.b,max.n.mono2.b,max.n.combi.b
#'Same as \code{max.n.[...].a} (where \code{[...]} is \code{mono1}, \code{mono2}, or \code{combi})
#'but for the second potentially simulated trial (suffix \code{.b}) of the respective trial type.
#'@param mtd.decision,mtd.decision.mono1.a,mtd.decision.mono2.a,mtd.decision.combi.a
#'Optional named lists that specify the rules for MTD selection. The determination of the MTD may occur independently for each simulated trial as soon as the conditions that
#'are specified in the corresponding \code{mtd.decision[...]} are met. Each list must have five named entries, namely, \code{target.prob}, \code{pat.at.mtd}, \code{min.pat}, \code{min.dlt}, and \code{rule}.
#'The names of the list entries are not case sensitive. For the actual MTD selection, the function will always first demand that the current dose level of the current trial
#'is recommended again, i.e. is also the recommended dose for the next cohort in this trial (the so-called stabilization criterion). If this is the case,
#'the conditions specified in the corresponding \code{mtd.decision[...]} are checked. If these conditions are also satisfied, the function will declare the current dose
#'of this trial to be the MTD and will not simulate any further cohorts for this trial.
#'
#'The entries are single numerical values with the following specifications and interpretations.
#'\itemize{
#'\item{\code{mtd.decision[...]$rule}\cr Integer, can only be \code{1} or \code{2} (the latter is the default). Specifies which of the remaining entries of \code{mtd.decision[...]} are used as mandatory conditions for MTD selection. For both rules,
#'two conditions are demanded as necessary conditions, while only one of the remaining two conditions must hold. The precise specifications are:
#'\itemize{
#'\item{\code{mtd.decision[...]$rule=1}\cr The conditions posed by the entries \code{mtd.decision[...]$pat.at.mtd} and \code{mtd.decision[...]$min.dlt} are mandatory, while only one of the conditions that are
#'posed by \code{mtd.decision[...]$min.pat} and \code{mtd.decision[...]$target.prob} must hold.}
#'\item{\code{mtd.decision[...]$rule=2}\cr The conditions from \code{mtd.decision[...]$min.pat} and \code{mtd.decision[...]$min.dlt} are mandatory, while only one of the conditions that are
#'posed by \code{mtd.decision[...]$pat.at.mtd} and \code{mtd.decision[...]$target.prob} must hold.}
#'}
#'}
#'\item{\code{mtd.decision[...]$min.pat}\cr Non-negative integer, defaults to \code{12}. Specifies the minimum number of patients that need to have been enrolled in a trial before the MTD can be declared.}
#'\item{\code{mtd.decision[...]$pat.at.mtd}\cr Non-negative integer, defaults to \code{6}. Specifies the minimum number of patients that need to have been treated with the current dose
#'(i.e., the current candidate for MTD determination) before the MTD can be declared.}
#'\item{\code{mtd.decision[...]$min.dlt}\cr Non-negative integer, defaults to \code{1}. Specifies the minimum number of DLTs that need to have been observed in a trial before the MTD can be declared.}
#'\item{\code{mtd.decision[...]$target.prob}\cr Non-negative real number strictly smaller than 1, defaults to \code{0.5}. Specifies a minimum probability
#'of having true DLT rate in the target interval that needs to be reached by the dose that is
#'to be declared the MTD. Note that this condition is always one of the two optional conditions, i.e., in case of the first rule, the condition does not need to hold
#'when sufficiently many patients were enrolled in total, while in case of the second rule, the condition does not need to hold when sufficiently many patients
#'were enrolled at the MTD.}
#'}
#'As before, the global parameter \code{mtd.decision} can be used to select a rule that is used for all simulated trials.
#'The trial-specific parameters, \code{mtd.decision.[...]}, default to the global value and can be adjusted to use different, trial-specific decision rules
#'for some or all of the trials.
#'@param mtd.decision.mono1.b,mtd.decision.mono2.b,mtd.decision.combi.b
#'Same as \code{mtd.decision.[...].a} (where \code{[...]} is \code{mono1}, \code{mono2}, or \code{combi})
#'but for the second potentially simulated trial (suffix \code{.b}) of the respective trial type.
#'@param mtd.enforce,mtd.enforce.mono1.a,mtd.enforce.mono2.a,mtd.enforce.combi.a
#'Optional logicals, default to \code{FALSE}. These parameters indicate whether the declaration of the MTD should be enforced once the maximum number
#'of patients is reached and the trial was not stopped yet. That is, if the parameter \code{mtd.enforce.[...]} for a specific trial is \code{TRUE},
#'the function will declare the next recommended dose as MTD once the maximum number of patients (\code{max.n.[...]}) was reached (provided that
#'the trial was not stopped early or found the MTD before reaching the maximum sample size). In particular,
#'these trials are not counted to belong to the category of trials that reached the maximum sample size without declaring an MTD. Further,
#'note that the value of \code{mtd.decision.[...]$min} can be set to a value greater than the maximum sample size, which would result in
#'only declaring an MTD if the maximum sample size is reached, but not at earlier points (as the condition on the minimal patient number could
#'not be fulfilled then).
#'
#'As before, there is a global variant of this parameter, namely \code{mtd.enforce}, which can be used to set this option for all trials at once.
#'To set the escalation rule only for a specific trial, use the trial-specific variants \code{mtd.enforce.[...]} with the corresponding suffix.
#'@param mtd.enforce.mono1.b,mtd.enforce.mono2.b,mtd.enforce.combi.b
#'Same as \code{mtd.enforce.[...].a} (where \code{[...]} is \code{mono1}, \code{mono2}, or \code{combi})
#'but for the second potentially simulated trial (suffix \code{.b}) of the respective trial type.
#'@param backfill.mono1.a,backfill.mono2.a,backfill.combi.a,backfill.mono1.b,backfill.mono2.b,backfill.combi.b
#'Optional logicals, indicating whether back-fill cohorts are simulated for the trial.
#'@param backfill.size,backfill.size.mono1.a,backfill.size.mono2.a,backfill.size.combi.a,backfill.size.mono1.b,backfill.size.mono2.b,backfill.size.combi.b
#'Optional numericals, provide the size of simulated back-fill cohorts. Interpreted in the same fashion as \code{cohort.size}.
#'@param backfill.prob,backfill.prob.mono1.a,backfill.prob.mono2.a,backfill.prob.combi.a,backfill.prob.mono1.b,backfill.prob.mono2.b,backfill.prob.combi.b
#'Optional numericals, provide the probabilities if multiple possible back-fill cohort sizes are given.
#'Interpreted in the same fashion as \code{cohort.prob}.
#'@param backfill.start.mono1.a,backfill.start.mono1.b,backfill.start.mono2.a,backfill.start.mono2.b,backfill.start.combi.a1,backfill.start.combi.a2,backfill.start.combi.b1,backfill.start.combi.b2
#'Optional numericals. Specify the first dose on which back-fill cohorts are enrolled. If not provided, the
#'lowest available dose will be assumed to be the starting point of back-fill cohorts in the respective trial.
#'@param n.studies Positive integer that specifies the number of studies to be simulated, defaults to \code{1}. Due to the long simulation time, it is recommended
#'to first try lower numbers to obtain an estimation of the run-time for larger numbers of studies. Typically, about 1000 studies are recommended to obtain
#'acceptably accurate simulation results.
#'@param n.cores Optional positive integer, defaults to \code{1}. Enables parallelized simulation and specifies the number of cores to be used. If a value greater
#'than 1 is given, the function will process the simulations in parallel using the given number of cores.
#'@param seed Optional positive integer that specified the seed to be used for the simulation. The default is \code{sample.int(.Machine$integer.max, 1)}.
#'Note that reproducibility can only be obtained when the function is executed on exactly the same computing architecture, using identical software versions
#'(e.g. of the compiler, Stan, and R), and the same input specifications. This is due to the internal use of \code{\link[rstan:rstan-package]{rstan-package}} for MCMC sampling, which is
#'only reproducible under these restrictions (refer to the Stan reference manual from \code{\link[rstan:rstan-package]{rstan-package}}'s \href{https://mc-stan.org}{homepage} for more detail).
#'@param chains Optional positive integer that specifies the number of Markov chains to be used during MCMC sampling, defaults to \code{4}. The parameter is given to the method
#'\code{\link[rstan:sampling]{rstan::sampling}()} from the package \code{\link[rstan:rstan-package]{rstan-package}}.
#'@param iter Optional positive integer that specifies the total number of iterations per chain to be used during MCMC sampling, defaults to \code{6000}. The parameter is given to the method
#'\code{\link[rstan:sampling]{rstan::sampling}()} from the package \code{\link[rstan:rstan-package]{rstan-package}}. Be aware that the first \code{warmup} iterations among the \code{iter}-many samples are discarded.
#'The function enforces \code{iter} to be at least \code{2000}, and further that \code{iter} is larger or equal than \code{warmup + 1000} (to ensure
#'that at least 1000 samples from each chain are kept).
#'@param warmup Optional positive integer that specifies the number of iterations per chain that are discarded from the total number of iterations, \code{iter}.
#'Defaults to \code{1000}. The parameter is given to the method \code{\link[rstan:sampling]{rstan::sampling}()} from the package \code{\link[rstan:rstan-package]{rstan-package}}.
#'The function enforces \code{warmup} to be at least \code{1000}, and further that \code{iter} is larger or equal than \code{warmup + 1000} (to ensure
#'that at least 1000 samples from each chain are kept).
#'@param adapt_delta Optional numerical that must be at least \code{0.8} and smaller than \code{1}, defaults to \code{0.9}. The parameter is given to the argument \code{control} of the \code{\link[rstan:sampling]{rstan::sampling}()}
#'method from the package \code{\link[rstan:rstan-package]{rstan-package}}. Translates to target Metropolis acceptance probability during Hamiltonian Monte Carlo, respectively NUTS.
#'Can be used to influence the \code{stepsize} Stan uses for leapfrog steps during the NUTS algorithm. See \code{\link[rstan:sampling]{rstan::sampling}()}, \code{\link[rstan:rstan-package]{rstan-package}} for details and further reading.
#'
#'Note: Larger values than the default can be used to reduce the number of Stan warnings on "divergent transition", but will slow down sampling. The function will not
#'permit to set \code{adapt_delta} to values below 0.8 (the \code{\link[rstan:rstan-package]{rstan-package}} default).
#'@param max_treedepth Optional integer that must be at least \code{10}, defaults to \code{15}. The parameter is given to the argument \code{control} of the \code{\link[rstan:sampling]{rstan::sampling}()}
#'method from the package \code{\link[rstan:rstan-package]{rstan-package}}. This is a parameter of the NUTS algorithm. Roughly speaking, NUTS constructs
#'a search tree when generating a new proposal, until a stopping criterion (the NUTS criterion) is satisfied or, alternatively, until
#'the maximum depth of the search tree is reached (to avoid endless looping). The argument \code{max_treedepth} allows to control the latter.
#'The maximum treedepth will also automatically constrain the maximum number of leapfrog steps, and should therefore not be below 10.
#'See also \code{\link[rstan:sampling]{rstan::sampling}()} for details on \code{max_treedepth}.
#'
#'Note: The default maximum treedepth of 15 should usually not be reached during sampling (and never was in tests, although this cannot be excluded in general), otherwise Stan will print a warning.
#'In this case, it is recommended to increase its value.
#'@param refresh Optional integer, defaults to \code{0}. The parameter is given to the \code{\link[rstan:sampling]{rstan::sampling}()}
#'method from the package \code{\link[rstan:rstan-package]{rstan-package}} and controls the amount of intermediate output during calls to \code{\link[rstan:sampling]{rstan::sampling}()}.
#'Only affects non-parallelized simulations (otherwise it only affects the output of each worker, which is discarded anyway). If a positive value is given, the number is
#'interpreted as the number of samples after which intermediate progress is reported during MCMC sampling (e.g., when set to 1000, Stan prints a progress report every 1000 samples).
#'If refresh is \code{0} (the default) or negative, no output will be printed during sampling. Changing the value of \code{refresh} does not
#'affect the results of the simulation, although it is recommended to not change it to a positive number: the function can lead to many
#'calls to \code{\link[rstan:sampling]{rstan::sampling}()} (multiple thousand when 1000 trials are simulated), so that intermediate output would
#'not be helpful or interpretable in most settings anyway.
#'@param file.name Optional character string that provides a name for potential output files. No output will be written if \code{file.name} or \code{path}
#'are not specified, and in this case, all outputs are only returned to R. If both are given, an excel sheet with the simulation results is written to disk.
#'@param path Optional character string that specifies the path to the output directory. No output will be written if \code{file.name} or \code{path}
#'are not specified, and in this case, all outputs are only returned to R. If both are given and \code{path} points to a directory, an excel sheet
#'an excel sheet that contains the simulation results is written to the specified location.
#'@param monitor.path Optional character string that specifies a path to an additional output directory for monitoring simulation progress.
#'More precisely, if \code{monitor.path} is not \code{NULL} (the default) and a valid \code{file.name} was given,
#'the function will write an output file for each completed trial to \code{monitor.path} (provided that it specifies a writable directory).
#'This can be used to monitor the progress of (parallelized) simulations, which may take multiple hours to simulate.
#'@param working.path Optional character string that specifies a path to a directory for temporary results. This directory is then used
#'to save MCCM results for a given data scenario under unique (hashed) file names, which allows re-loading MCMC runs whenever the same
#'data scenario occurs in some later trial. This also works for parallelized simulations. File names of temporary results will use the given
#'\code{file.name} concatenated with "_tmp" as prefix. Do not modify such files manually to avoid errors. To avoid reloading results from previous simulations (with potentially different
#'MCMC settings or priors), the function will check for files whose name contains the \code{file.name} and "_tmp" from the directory
#'specified in \code{working.path} before the simulation.
#'If such files are detected, an error will be thrown. To allow the function
#'to automatically delete all such files from the \code{working.path} before and after the simulation, consider using the argument \code{clean.working.path}.
#'@param clean.working.path Optional logical, defaults to \code{FALSE}. Indicates whether the function is allowed to delete any existing files that contain
#'\code{file.name} and "_tmp" in their name from the directory
#'specified in \code{working.path}. If \code{TRUE}, such temporary files will be automatically removed from the working directory before
#'and after simulation. If \code{FALSE}, the function will check for such file names initially and throw an error when such files are found.
#'@param output.sim.config Optional logical that specifies whether the input parameters of the call to \code{sim_jointBLRM()} are added to
#'the output list (or, if activated, the output files written to disk). Defaults to \code{TRUE}. This can be used to save e.g. the seed, dose-toxicity scenarios, and
#'priors of the simulation to allow double-checking the configuration of the simulation at a later point in time.
#'
#'@returns List that contains a number of metrics that summarize the results of the simulation for each simulated trial,
#'and, depending on the specification, additional list entries that save the given input specification.
#'
#'By default, the following list entries are generated for each simulated trial (where \code{[...]}
#'is the corresponding suffix of the trial):
#'\itemize{
#'\item{\code{...$'results [...]'}\cr Overview of number of MTDs per dosing interval (under, target, over), and number of
#' trials that stopped without finding an MTD, either due to the stopping rule (when
#' EWOC-based escalation is used) or due to reaching the maximal patient number.}
#'
#'\item{\code{...$'summary [...]'}\cr Summary statistics (mean, median, min, max, 2.5% and 97.5% quantile) of
#'the number of patients, patients per dosing interval, DLTs, and
#'DLTs per dosing interval.}
#'
#'\item{\code{...$'MTDs [...]'}\cr Number of MTDs per dose level and their assumed DLT rates.}
#'
#'\item{\code{...$'#pat [...]'}\cr Mean, median, minimum and maximum number of patients enrolled at each dose level.}
#'}
#'Additionally, if \code{output.sim.config} is active, the list will include the following entries:
#'\itemize{
#'\item{\code{...$'historical.data'}\cr Only included when historical data was included in the simulations. Contains the specified cohorts
#' from historical data.}
#'
#'\item{\code{...$'prior'}\cr The prior distribution for the hyper-parameters that was used by the BLRM.}
#'
#'\item{\code{...$'specifications'}\cr All other simulation specifications. Includes e. g. reference doses, decision rules, escalation steps,
#'and the seed.}
#'
#'\item{\code{...$'Stan options'}\cr Options given to Stan, such as the number of MCMC iterations and chains.}
#'}
#'
#'
#'@details
#'The joint BLRM is defined according to (Neuenschwander et al., 2014 and 2016). It allows to perform Bayesian logistic regression
#'to estimate the dose-toxicity relationship of two different monotherapies and combination therapy with these compounds in
#'a joint model, which includes a hierarchical prior for robust borrowing across trials. Refer to the documentation of
#'\code{\link[OncoBLRM:scenario_jointBLRM]{scenario_jointBLRM}()} for a detailed model description.
#'
#'@references
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
#' Babb, J., Rogatko, A., & Zacks, S. (1998). Cancer phase I clinical trials: Efficient dose escalation with overdose control.
#' Statistics in medicine 17(10), 1103-1120.
#'
#' Zhou, H.,  Yuan, Y., & Nie, L. (2018). Accuracy, safety, and reliability of novel phase I designs.
#' Clinical Cancer Research, 24(21), 5483-5484 <doi: 10.1158/1078-0432.ccr-18-0168>.
#'
#'
#' @seealso \code{\link[OncoBLRM:scenario_jointBLRM]{scenario_jointBLRM}()}, \code{\link[rstan:sampling]{rstan::sampling}()},
#' \code{\link[rstan:rstan-package]{rstan-package}}.
#'
#'@export
sim_jointBLRM <- function(active.mono1.a = FALSE,
                          active.mono1.b = FALSE,
                          active.mono2.a = FALSE,
                          active.mono2.b = FALSE,
                          active.combi.a = FALSE,
                          active.combi.b = FALSE,

                          doses.mono1.a,
                          doses.mono2.a,
                          doses.mono1.b,
                          doses.mono2.b,
                          doses.combi.a,
                          doses.combi.b,
                          dose.ref1,
                          dose.ref2,
                          tox.mono1.a,
                          tox.mono2.a,
                          tox.combi.a,
                          tox.mono1.b,
                          tox.mono2.b,
                          tox.combi.b,

                          start.dose.mono1.a,
                          start.dose.mono2.a,
                          start.dose.mono1.b,
                          start.dose.mono2.b,
                          start.dose.combi.a1,
                          start.dose.combi.a2,
                          start.dose.combi.b1,
                          start.dose.combi.b2,

                          cohort.queue = rep( c(1,2,3,4,5,6), times= 100),
                          historical.data = NULL,

                          esc.rule = "ewoc",
                          esc.comp.max=1,
                          dosing.intervals = c(0.16, 0.33, 0.6),
                          ewoc.threshold = 0.25,
                          loss.weights = c(1, 0, 1, 2),
                          dynamic.weights =  rbind(c(0.32, 0, 0.32, 0.36),
                                                   c(0.29, 0, 0.31, 0.4),
                                                   c(0.27, 0, 0.33, 0.4),
                                                   c(0.2,  0, 0.3,  0.5)),

                          prior.mu = list(
                            mu_a1 = c(logit(0.33), 2),
                            mu_b1 = c(0, 1),
                            mu_a2 = c(logit(0.33), 2),
                            mu_b2 = c(0, 1),
                            mu_eta = c(0, 1.121)
                          ),

                          prior.tau = list(tau_a1 = c(log(0.25), log(4)/1.96),
                                           tau_b1 = c(log(0.125), log(4)/1.96),
                                           tau_a2 = c(log(0.25), log(4)/1.96),
                                           tau_b2 = c(log(0.125), log(4)/1.96),
                                           tau_eta = c(log(0.125),log(4)/1.96)),
                          saturating = FALSE,
                          esc.step = NULL,
                          esc.step.mono1.a = esc.step,
                          esc.step.mono2.a = esc.step,
                          esc.step.mono1.b = esc.step,
                          esc.step.mono2.b = esc.step,
                          esc.step.combi.a1 = esc.step,
                          esc.step.combi.b1 = esc.step,
                          esc.step.combi.a2 = esc.step,
                          esc.step.combi.b2 = esc.step,

                          esc.constrain = FALSE,
                          esc.constrain.mono1.a=esc.constrain,
                          esc.constrain.mono2.a=esc.constrain,
                          esc.constrain.mono1.b=esc.constrain,
                          esc.constrain.mono2.b=esc.constrain,
                          esc.constrain.combi.a1=esc.constrain,
                          esc.constrain.combi.b1=esc.constrain,
                          esc.constrain.combi.a2=esc.constrain,
                          esc.constrain.combi.b2=esc.constrain,

                          cohort.size = c(3),
                          cohort.size.mono1.a = cohort.size,
                          cohort.size.mono1.b = cohort.size,
                          cohort.size.mono2.a = cohort.size,
                          cohort.size.mono2.b = cohort.size,
                          cohort.size.combi.a = cohort.size,
                          cohort.size.combi.b = cohort.size,


                          cohort.prob = NULL,
                          cohort.prob.mono1.a = cohort.prob,
                          cohort.prob.mono1.b = cohort.prob,
                          cohort.prob.mono2.a = cohort.prob,
                          cohort.prob.mono2.b = cohort.prob,
                          cohort.prob.combi.a = cohort.prob,
                          cohort.prob.combi.b = cohort.prob,


                          max.n = 42,
                          max.n.mono1.a = max.n,
                          max.n.mono1.b = max.n,
                          max.n.mono2.a = max.n,
                          max.n.mono2.b = max.n,
                          max.n.combi.a = max.n,
                          max.n.combi.b = max.n,


                          mtd.decision = list(target.prob = 0.5,
                                              pat.at.mtd = 6,
                                              min.pat = 12,
                                              min.dlt = 1,
                                              rule = 2),
                          mtd.decision.combi.a = mtd.decision,
                          mtd.decision.combi.b = mtd.decision,
                          mtd.decision.mono1.a = mtd.decision,
                          mtd.decision.mono1.b = mtd.decision,
                          mtd.decision.mono2.a = mtd.decision,
                          mtd.decision.mono2.b = mtd.decision,

                          mtd.enforce = FALSE,
                          mtd.enforce.mono1.a = mtd.enforce,
                          mtd.enforce.mono2.a = mtd.enforce,
                          mtd.enforce.mono1.b = mtd.enforce,
                          mtd.enforce.mono2.b = mtd.enforce,
                          mtd.enforce.combi.a = mtd.enforce,
                          mtd.enforce.combi.b = mtd.enforce,

                          backfill.mono1.a = FALSE,
                          backfill.mono1.b = FALSE,
                          backfill.mono2.a = FALSE,
                          backfill.mono2.b = FALSE,
                          backfill.combi.a = FALSE,
                          backfill.combi.b = FALSE,

                          backfill.size = c(3),
                          backfill.prob = NULL,
                          backfill.size.mono1.a = backfill.size,
                          backfill.size.mono1.b = backfill.size,
                          backfill.size.mono2.a = backfill.size,
                          backfill.size.mono2.b = backfill.size,
                          backfill.size.combi.a = backfill.size,
                          backfill.size.combi.b = backfill.size,

                          backfill.prob.mono1.a = backfill.prob,
                          backfill.prob.mono1.b = backfill.prob,
                          backfill.prob.mono2.a = backfill.prob,
                          backfill.prob.mono2.b = backfill.prob,
                          backfill.prob.combi.a = backfill.prob,
                          backfill.prob.combi.b = backfill.prob,

                          backfill.start.mono1.a = NULL,
                          backfill.start.mono1.b = NULL,
                          backfill.start.mono2.a = NULL,
                          backfill.start.mono2.b = NULL,
                          backfill.start.combi.a1 = NULL,
                          backfill.start.combi.a2 = NULL,
                          backfill.start.combi.b1 = NULL,
                          backfill.start.combi.b2 = NULL,

                          n.studies = 1,
                          n.cores = 1,
                          seed = sample.int(.Machine$integer.max, 1),
                          chains = 4,
                          iter = 6000,
                          warmup = 1000,
                          adapt_delta = 0.9,
                          max_treedepth = 12,
                          refresh=0,
                          file.name = NULL,
                          #path to the output file
                          path = NULL,
                          monitor.path = NULL,
                          working.path = NULL,
                          clean.working.path = FALSE,
                          output.sim.config =TRUE
){
  #"**************************************************************************************"


  #---------------------------------------------------------------------------------------
  #Input checks before function start
  #---------------------------------------------------------------------------------------
  #Logical values are easily checked
  if(!is.logical(c(active.mono1.a, active.mono1.b, active.mono2.a, active.mono2.b,
                   active.combi.a,active.combi.b,  output.sim.config))){
    stop("The arguments `active.mono1.a`, `active.mono1.b`, `active.mono2.a`, `active.mono2.b`,\n",
         "`active.combi.a`, `active.combi.b`, and `output.sim.config` must \n",
         "be of type logical.")
  }
  if(!is.logical(c(esc.constrain, esc.constrain.mono1.a, esc.constrain.mono1.b,
                   esc.constrain.mono2.a, esc.constrain.mono2.b,
                   esc.constrain.combi.a1, esc.constrain.combi.a2,
                   esc.constrain.combi.b1, esc.constrain.combi.b2))){
    stop("The arguments `esc.constrain`, `esc.constrain.mono1.a`, `esc.constrain.mono1.b`,\n",
          "`esc.constrain.mono2.a`, `esc.constrain.mono2.b`, `esc.constrain.combi.a1`,\n",
          "`esc.constrain.combi.a2`, `esc.constrain.combi.b1`, `esc.constrain.combi.b2`\n",
          "must be of type logical.")
  }

  if(!is.logical(c(mtd.enforce, mtd.enforce.mono1.a, mtd.enforce.mono1.b,
                   mtd.enforce.mono2.a, mtd.enforce.mono2.b,
                   mtd.enforce.combi.a,
                   mtd.enforce.combi.b))){
    stop("The arguments `mtd.enforce`, `mtd.enforce.mono1.a`, `mtd.enforce.mono1.b`,\n",
         "`mtd.enforce.mono2.a`, `mtd.enforce.mono2.b`, `mtd.enforce.combi.a`,\n",
         "and `mtd.enforce.combi.b` must be of type logical.")
  }
  if(!any(c(active.mono1.a, active.mono1.b, active.mono2.a, active.mono2.b,
            active.combi.a,active.combi.b))){
    stop("None of the trials is activated (via arguments `active.[...]`).")
  }

  if(!is.logical(saturating)){
    stop("`saturating` must be of type logical.")
  }
  if(!is.logical(clean.working.path)){
    stop("`clean.working.path` must be of type logical.")
  }

  if(!is.logical(
    c(backfill.mono1.a,backfill.mono1.b,backfill.mono2.a,backfill.mono2.b,
      backfill.combi.a,backfill.combi.b)
  )){
    stop("The arguments `backfill.mono1.a`, `backfill.mono1.b`, `backfill.mono2.a`,\n",
         "`backfill.mono2.b`, `backfill.combi.a`, `backfill.combi.b` must be logical.")
  }
  #--------------------------------------------
  #Checks for integers and integer vectors
  #--------------------------------------------
  if(!is.wholenumber(refresh)){
    stop("`refresh` needs to be an integer.")
  }

  if(!is.wholenumber(max_treedepth)){
    stop("`max_treedepth` needs to be an integer.")
  }
  if(!max_treedepth>=10 | !length(max_treedepth)==1){
    stop("`max_treedepth`` must have length 1 and be larger than 10.")
  }

  if(!is.wholenumber(seed)){
    stop("`seed` needs to be a whole number.")
  }
  if(!length(seed)==1){
    stop("`seed` must have length 1.")
  }

  if(!is.wholenumber(chains)){
    stop("`chains` needs to be an integer.")
  }
  if(!chains>=1 | !length(chains)==1){
    stop("`chains` must have length 1 and be larger than 0.")
  }

  if(!is.wholenumber(iter)){
    stop("`iter` needs to be an integer.")
  }

  if(!iter>=2000 | !length(iter)==1){
    stop("`iter` must have length 1 and be larger or equal than 2000.")
  }

  if(!is.wholenumber(warmup)){
    stop("`warmup` needs to be an integer.")
  }
  if(!warmup>=1000 | !length(warmup)==1){
    stop("`warmup` must have length 1 and be larger or equal than 1000.")
  }
  if(!(iter-warmup >=1000)){
    stop("`iter` must be larger than `warmup` by at least 1000.")
  }

  if(!is.wholenumber(n.cores)){
    stop("`n.cores` needs to be an integer.")
  }
  if(!n.cores>=0 | !length(n.cores)==1){
    stop("`n.cores` must have length 1 and be larger or equal than 0.")
  }

  if(!is.wholenumber(max.n)){
    stop("`max.n` needs to be an integer.")
  }

  if(!max.n>=0 | !length(max.n)==1){
    stop("`max.n` must have length 1 and be larger than 0.")
  }

  if(!is.wholenumber(max.n.mono1.a)){
    stop("`max.n.mono1.a` needs to be an integer.")
  }
  if(!max.n.mono1.a>=0 | !length(max.n.mono1.a)==1){
    stop("`max.n.mono1.a` must have length 1 and be larger than 0.")
  }

  if(!is.wholenumber(max.n.mono1.b)){
    stop("`max.n.mono1.b` needs to be an integer.")
  }
  if(!max.n.mono1.b>=0 | !length(max.n.mono1.b)==1){
    stop("`max.n.mono1.b` must have length 1 and be larger than 0.")
  }

  if(!is.wholenumber(max.n.mono2.a)){
    stop("`max.n.mono2.a` needs to be an integer.")
  }

  if(!max.n.mono2.a>=0 | !length(max.n.mono2.a)==1){
    stop("`max.n.mono2.a` must have length 1 and be larger than 0.")
  }

  if(!is.wholenumber(max.n.mono2.b)){
    stop("`max.n.mono2.b` needs to be an integer.")
  }

  if(!max.n.mono2.b>=0 | !length(max.n.mono2.b)==1){
    stop("`max.n.mono2.b` must have length 1 and be larger than 0.")
  }


  if(!is.wholenumber(max.n.combi.a)){
    stop("`max.n.combi.a` needs to be an integer.")
  }

  if(!max.n.combi.a>=0 | !length(max.n.combi.a)==1){
    stop("`max.n.combi.a` must have length 1 and be larger than 0.")
  }

  if(!is.wholenumber(max.n.combi.b)){
    stop("`max.n.combi.b` needs to be an integer.")
  }

  if(!max.n.combi.b>=0 | !length(max.n.combi.b)==1){
    stop("`max.n.combi.b` must have length 1 and be larger than 0.")
  }


  if(!is.wholenumber(cohort.size)){
    stop("`cohort.size needs` to be an integer.")
  }
  if(!all(cohort.size>=1)){
    stop("Entries of `cohort.size` must be at least 1.")
  }

  if(!is.wholenumber(cohort.size.mono1.a)){
    stop("`cohort.size.mono1.a` needs to be an integer.")
  }
  if(!all(cohort.size.mono1.a>=1)){
    stop("Entries of `cohort.size.mono1.a` must be at least 1.")
  }

  if(!is.wholenumber(cohort.size.mono1.b)){
    stop("`cohort.size.mono1.b` needs to be an integer.")
  }
  if(!all(cohort.size.mono1.b>=1)){
    stop("Entries of `cohort.size.mono1.b` must be at least 1.")
  }

  if(!is.wholenumber(cohort.size.mono2.a)){
    stop("`cohort.size.mono2.a` needs to be an integer.")
  }
  if(!all(cohort.size.mono2.a>=1)){
    stop("Entries of `cohort.size.mono2.a` must be at least 1.")
  }

  if(!is.wholenumber(cohort.size.mono2.b)){
    stop("`cohort.size.mono2.b` needs to be an integer.")
  }
  if(!all(cohort.size.mono2.b>=1)){
    stop("Entries of `cohort.size.mono2.b` must be at least 1.")
  }

  if(!is.wholenumber(cohort.size.combi.a)){
    stop("`cohort.size.combi.a` needs to be an integer.")
  }
  if(!all(cohort.size.combi.a>=1)){
    stop("Entries of `cohort.size.combi.a` must be at least 1.")
  }

  if(!is.wholenumber(cohort.size.combi.b)){
    stop("`cohort.size.combi.b` needs to be an integer.")
  }
  if(!all(cohort.size.combi.b>=1)){
    stop("Entries of `cohort.size.combi.b` must be at least 1.")
  }

  #same for backfill.size
  if(!is.wholenumber(backfill.size)){
    stop("`backfill.size needs` to be an integer.")
  }
  if(!all(backfill.size>=0)){
    stop("Entries of `backfill.size` must be at least 0.")
  }

  if(!is.wholenumber(backfill.size.mono1.a)){
    stop("`backfill.size.mono1.a` needs to be an integer.")
  }
  if(!all(backfill.size.mono1.a>=0)){
    stop("Entries of `backfill.size.mono1.a` must be at least 0.")
  }

  if(!is.wholenumber(backfill.size.mono1.b)){
    stop("`backfill.size.mono1.b` needs to be an integer.")
  }
  if(!all(backfill.size.mono1.b>=0)){
    stop("Entries of `backfill.size.mono1.b` must be at least 0.")
  }

  if(!is.wholenumber(backfill.size.mono2.a)){
    stop("`backfill.size.mono2.a` needs to be an integer.")
  }
  if(!all(backfill.size.mono2.a>=0)){
    stop("Entries of `backfill.size.mono2.a` must be at least 0.")
  }

  if(!is.wholenumber(backfill.size.mono2.b)){
    stop("`backfill.size.mono2.b` needs to be an integer.")
  }
  if(!all(backfill.size.mono2.b>=0)){
    stop("Entries of `backfill.size.mono2.b` must be at least 0.")
  }

  if(!is.wholenumber(backfill.size.combi.a)){
    stop("`backfill.size.combi.a` needs to be an integer.")
  }
  if(!all(backfill.size.combi.a>=0)){
    stop("Entries of `backfill.size.combi.a` must be at least 0.")
  }

  if(!is.wholenumber(backfill.size.combi.b)){
    stop("`backfill.size.combi.b` needs to be an integer.")
  }
  if(!all(backfill.size.combi.b>=0)){
    stop("Entries of `backfill.size.combi.b` must be at least 0.")
  }

  if(!is.wholenumber(n.studies)){
    stop("`n.studies` needs to be an integer.")
  }
  if(!n.studies>=1){
    stop("`n.studies` must be at least 1.")
  }

  if(!is.wholenumber(esc.comp.max)){
    stop("`esc.comp.max` needs to be an integer.")
  }
  if(!esc.comp.max%in%c(1, 2)){
    stop("`esc.comp.max` must be either 1 or 2.")
  }
  #--------------------------------------------
  #Checks for strings
  #--------------------------------------------
  if(!is.null(path) & !is.null(file.name)){
    if(!dir.exists(file.path(path))){
      message("NOTE: specified `path` could not be opened, results will only be returned to R.")
    }
    if(!is.null(monitor.path)){
      if(!dir.exists(file.path(monitor.path))){
        message("NOTE: specified `monitor.path` for monitoring progress could not be opened.\n",
                "      No intermediate output will be written.")
      }
    }
  }

  if(!is.character(esc.rule)){
    stop("`esc.rule` must be a character type.")
  }
  if(!length(esc.rule)==1){
    stop("`esc.rule` must be of length 1.")
  }
  if(!esc.rule%in%c("ewoc", "ewoc.opt", "ewoc.max", "loss", "dynamic.loss", "dynamic")){
    stop("`esc.rule` can only be \"ewoc\", \"ewoc.max\", \"ewoc.opt\", \"loss\",\n",
         "\"dynamic\", or \"dynamic.loss\".")
  }
  if(esc.rule=="dynamic"){
    esc.rule <- "dynamic.loss"
  }
  if(!is.character(cohort.queue) & !is.numeric(cohort.queue)){
    stop("Only character or numeric are permitted as entries of `cohort.queue`.")
  }
  cohort.queue <- tolower(cohort.queue)
  if(!all(cohort.queue%in% c(1, 2, 3, 4, 5, 6,
                             "mono1.a", "mono1.b",
                             "mono2.a", "mono2.b",
                             "combi.a", "combi.b"))){
    stop("`cohort.queue` may only contain the following entries:\n",
         "1, 2, 3, 4, 5, 6, \"mono1.a\", \"mono1.b\", \"mono2.a\", \"mono2.b\"\n",
         "\"combi.a\", or \"combi.b\".")
  }
  #exchanges all character entries with the corresponding study numbers
  cohort.queue <- harmonize_vecnames_jointBLRM(cohort.queue)

  #--------------------------------------------
  #Checks for numeric and numeric vectors
  #--------------------------------------------
  #general numeric values
  if(!is.num(adapt_delta, len=1, low=0.8, up=1, uB=F, lB=T)){
    stop("`adapt_delta` must be a number that is at least 0.8 and smaller than 1")
  }

  if(is.null(cohort.prob)){
    cohort.prob <- rep(1/length(cohort.size), times = length(cohort.size))
  }

  if(is.null(cohort.prob.mono1.a)){
    cohort.prob.mono1.a <- rep(1/length(cohort.size.mono1.a),
                               times = length(cohort.size.mono1.a))
  }

  if(is.null(cohort.prob.mono1.b)){
    cohort.prob.mono1.b <- rep(1/length(cohort.size.mono1.b),
                               times = length(cohort.size.mono1.b))
  }

  if(is.null(cohort.prob.mono2.a)){
    cohort.prob.mono2.a <- rep(1/length(cohort.size.mono2.a),
                               times = length(cohort.size.mono2.a))
  }

  if(is.null(cohort.prob.mono2.b)){
    cohort.prob.mono2.b <- rep(1/length(cohort.size.mono2.b),
                               times = length(cohort.size.mono2.b))
  }

  if(is.null(cohort.prob.combi.a)){
    cohort.prob.combi.a <- rep(1/length(cohort.size.combi.a),
                               times = length(cohort.size.combi.a))
  }

  if(is.null(cohort.prob.combi.b)){
    cohort.prob.combi.b <- rep(1/length(cohort.size.combi.b),
                               times = length(cohort.size.combi.b))
  }

  if(!is.num(cohort.prob, len=length(cohort.size),
            low=0, up=1, uB=T, lB=F)){
    stop("`cohort.prob` must be of the same length as `cohort.size` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }

  if(!is.num(cohort.prob.mono1.a, len=length(cohort.size.mono1.a),
            low=0, up=1, uB=T, lB=F)){
    stop("`cohort.prob.mono1.a` must be of the same length as `cohort.size.mono1.a` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }


  if(!is.num(cohort.prob.mono1.b, len=length(cohort.size.mono1.b),
             low=0, up=1, uB=T, lB=F)){
    stop("`cohort.prob.mono1.b` must be of the same length as `cohort.size.mono1.b` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }

  if(!is.num(cohort.prob.mono2.a, len=length(cohort.size.mono2.a),
             low=0, up=1, uB=T, lB=F)){
    stop("`cohort.prob.mono2.a` must be of the same length as `cohort.size.mono2.a` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }

  if(!is.num(cohort.prob.mono2.b, len=length(cohort.size.mono2.b),
             low=0, up=1, uB=T, lB=F)){
    stop("`cohort.prob.mono2.b` must be of the same length as `cohort.size.mono2.b` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }

  if(!is.num(cohort.prob.combi.a, len=length(cohort.size.combi.a),
             low=0, up=1, uB=T, lB=F)){
    stop("`cohort.prob.combi.a` must be of the same length as `cohort.size.combi.a` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }

  if(!is.num(cohort.prob.combi.b, len=length(cohort.size.combi.b),
             low=0, up=1, uB=T, lB=F)){
    stop("`cohort.prob.combi.b` must be of the same length as `cohort.size.combi.b` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }


  #test for backfill probs
  if(is.null(backfill.prob)){
    backfill.prob <- rep(1/length(backfill.size), times = length(backfill.size))
  }

  if(is.null(backfill.prob.mono1.a)){
    backfill.prob.mono1.a <- rep(1/length(backfill.size.mono1.a),
                               times = length(backfill.size.mono1.a))
  }

  if(is.null(backfill.prob.mono1.b)){
    backfill.prob.mono1.b <- rep(1/length(backfill.size.mono1.b),
                               times = length(backfill.size.mono1.b))
  }

  if(is.null(backfill.prob.mono2.a)){
    backfill.prob.mono2.a <- rep(1/length(backfill.size.mono2.a),
                               times = length(backfill.size.mono2.a))
  }

  if(is.null(backfill.prob.mono2.b)){
    backfill.prob.mono2.b <- rep(1/length(backfill.size.mono2.b),
                               times = length(backfill.size.mono2.b))
  }

  if(is.null(backfill.prob.combi.a)){
    backfill.prob.combi.a <- rep(1/length(backfill.size.combi.a),
                               times = length(backfill.size.combi.a))
  }

  if(is.null(backfill.prob.combi.b)){
    backfill.prob.combi.b <- rep(1/length(backfill.size.combi.b),
                               times = length(backfill.size.combi.b))
  }

  if(!is.num(backfill.prob, len=length(backfill.size),
             low=0, up=1, uB=T, lB=F)){
    stop("`backfill.prob` must be of the same length as `backfill.size` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }

  if(!is.num(backfill.prob.mono1.a, len=length(backfill.size.mono1.a),
             low=0, up=1, uB=T, lB=F)){
    stop("`backfill.prob.mono1.a` must be of the same length as `backfill.size.mono1.a` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }


  if(!is.num(backfill.prob.mono1.b, len=length(backfill.size.mono1.b),
             low=0, up=1, uB=T, lB=F)){
    stop("`backfill.prob.mono1.b` must be of the same length as `backfill.size.mono1.b` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }

  if(!is.num(backfill.prob.mono2.a, len=length(backfill.size.mono2.a),
             low=0, up=1, uB=T, lB=F)){
    stop("`backfill.prob.mono2.a` must be of the same length as `backfill.size.mono2.a` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }

  if(!is.num(backfill.prob.mono2.b, len=length(backfill.size.mono2.b),
             low=0, up=1, uB=T, lB=F)){
    stop("`backfill.prob.mono2.b` must be of the same length as `backfill.size.mono2.b` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }

  if(!is.num(backfill.prob.combi.a, len=length(backfill.size.combi.a),
             low=0, up=1, uB=T, lB=F)){
    stop("`backfill.prob.combi.a` must be of the same length as `backfill.size.combi.a` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }

  if(!is.num(backfill.prob.combi.b, len=length(backfill.size.combi.b),
             low=0, up=1, uB=T, lB=F)){
    stop("`backfill.prob.combi.b` must be of the same length as `backfill.size.combi.b` and consist of \n",
         "numbers larger than 0 and smaller or equal to 1.")
  }

  #tests for dosing intervals
  if(!is.num(dosing.intervals, low=0, up=1,
             lB=T, uB=T)){
    stop("`dosing.intervals` must contain entries between 0 and 1.")
  }
  #if(all(dosing.intervals==dosing.intervals[1])){
  #  stop("dosing.intervals must contain at least two different numbers.")
  #}
  if(!length(dosing.intervals)%in%c(1, 2, 3)){
    stop("`dosing.intervals` must have length 1, 2, or 3.")
  }
  if(esc.rule%in%c("ewoc", "ewoc.opt", "ewoc.max")){
    if(length(dosing.intervals)==1){
      if(!dosing.intervals[1]>0){
        stop("When `dosing.intervals` has length 1, its entry must be strictly larger than 0.")
      }
      dosing.intervals <- c(0, dosing.intervals[1])
    }
    dosing.intervals <- dosing.intervals[1:2]
    if(dosing.intervals[2]==1){
      stop("1 is not permitted as upper boundary of the target interval specified in `dosing.intervals`.")
    }
    if(!dosing.intervals[1]<dosing.intervals[2]){
      stop("The entries of `dosing.interval` need to be strictly ascending when `esc.rule` is \"ewoc\".")
    }
  }else{
    if(!length(dosing.intervals)==3){
      stop("`dosing.intervals` must have length 3 when `esc.rule` is set to \"loss\" or \"dynamic.loss\".")
    }
    if(dosing.intervals[2]==1){
      stop("1 is not permitted as upper boundary of the target interval specified in `dosing.intervals`.")
    }
    if(dosing.intervals[1]>=dosing.intervals[2] |
       dosing.intervals[2]>=dosing.intervals[3] ){
      stop("The entries of `dosing.interval` need to be strictly ascending.")
    }
  }


  #----------------------------
  #consistency checks for doses
  #----------------------------

  #check whether sufficiently many cohorts are specified in cohort queue
  if(active.mono1.a){
    min.cs <- min(cohort.size.mono1.a)
    n.occ <- length(which(cohort.queue==1))
    if(!n.occ>0){
      stop("Trial `mono1.a` is active (`active.mono1.a`) but `cohort.queue` does not contain cohorts for it.")
    }
    if(!min.cs*n.occ >= max.n.mono1.a){
      message("NOTE:\n",
              "`cohort.queue` does not contain sufficiently many cohorts for trial `mono1.a`\n",
              "to guarantee that `max.n.mono1.a` can be reached (based on minimal cohort size).\n",
              "The pattern specified in the `cohort.queue `will be repeated until sufficiently\n",
              "many cohorts are included in the `cohort.queue`.\n",
              "To avoid this warning, supply a `cohort.queue` that contains as least ",
              paste0(ceiling(max.n.mono1.a/min.cs)), "\n cohorts for `mono1.a`.")
      cohort.queue <- rep(cohort.queue, times = ceiling(ceiling(max.n.mono1.a/min.cs)/n.occ)+1)
    }
  }
  if(active.mono1.b){
    min.cs <- min(cohort.size.mono1.b)
    n.occ <- length(which(cohort.queue==4))
    if(!n.occ>0){
      stop("Trial `mono1.b` is active (`active.mono1.b`) but `cohort.queue` does not contain cohorts for it.")
    }
    if(!min.cs*n.occ >= max.n.mono1.b){
      message("NOTE:\n",
              "`cohort.queue` does not contain sufficiently many cohorts for trial `mono1.b`\n",
              "to guarantee that `max.n.mono1.b` can be reached (based on minimal cohort size).\n",
              "The pattern specified in the `cohort.queue` will be repeated until sufficiently\n",
              "many cohorts are included in the `cohort.queue`.\n",
              "To avoid this warning, supply a `cohort.queue` that contains as least ",
              paste0(ceiling(max.n.mono1.b/min.cs)), "\n cohorts for `mono1.b`.")
      cohort.queue <- rep(cohort.queue, times = ceiling(ceiling(max.n.mono1.b/min.cs)/n.occ)+1)
    }
  }
  if(active.mono2.a){
    min.cs <- min(cohort.size.mono2.a)
    n.occ <- length(which(cohort.queue==2))
    if(!n.occ>0){
      stop("Trial `mono2.a` is active (`active.mono2.a`) but `cohort.queue` does not contain cohorts for it.")
    }
    if(!min.cs*n.occ >= max.n.mono2.a){
      message("NOTE:\n",
              "`cohort.queue` does not contain sufficiently many cohorts for trial `mono2.a`\n",
              "to guarantee that `max.n.mono2.a` can be reached (based on minimal cohort size).\n",
              "The pattern specified in the `cohort.queue` will be repeated until sufficiently\n",
              "many cohorts are included in the `cohort.queue`.\n",
              "To avoid this warning, supply a `cohort.queue` that contains as least ",
              paste0(ceiling(max.n.mono2.a/min.cs)), "\n cohorts for `mono2.a`.")
      cohort.queue <- rep(cohort.queue, times = ceiling(ceiling(max.n.mono2.a/min.cs)/n.occ)+1)
    }
  }

  if(active.mono2.b){
    min.cs <- min(cohort.size.mono2.b)
    n.occ <- length(which(cohort.queue==5))
    if(!n.occ>0){
      stop("Trial `mono2.b` is active (`active.mono2.b`) but the `cohort.queue` does not contain cohorts for it.")
    }
    if(!min.cs*n.occ >= max.n.mono2.b){
      message("NOTE:\n",
              "`cohort.queue` does not contain sufficiently many cohorts for trial `mono2.b`\n",
              "to guarantee that `max.n.mono2.b` can be reached (based on minimal cohort size).\n",
              "The pattern specified in the `cohort.queue` will be repeated until sufficiently\n",
              "many cohorts are included in the `cohort.queue`.\n",
              "To avoid this warning, supply a `cohort.queue` that contains as least ",
              paste0(ceiling(max.n.mono2.b/min.cs)), "\n cohorts for `mono2.b`.")
      cohort.queue <- rep(cohort.queue, times = ceiling(ceiling(max.n.mono2.b/min.cs)/n.occ)+1)
    }
  }

  if(active.combi.a){
    min.cs <- min(cohort.size.combi.a)
    n.occ <- length(which(cohort.queue==3))
    if(!n.occ>0){
      stop("Trial `combi.a` is active (`active.combi.a`) but the `cohort.queue` does not contain cohorts for it.")
    }
    if(!min.cs*n.occ >= max.n.combi.a){
      message("NOTE:\n",
              "`cohort.queue` does not contain sufficiently many cohorts for trial `combi.a`\n",
              "to guarantee that `max.n.combi.a` can be reached (based on minimal cohort size).\n",
              "The pattern specified in the `cohort.queue` will be repeated until sufficiently\n",
              "many cohorts are included in the `cohort.queue`.\n",
              "To avoid this warning, supply a `cohort.queue` that contains as least ",
              paste0(ceiling(max.n.combi.a/min.cs)), "\n cohorts for `combi.a`.")
      cohort.queue <- rep(cohort.queue, times = ceiling(ceiling(max.n.combi.a/min.cs)/n.occ)+1)
    }
  }
  if(active.combi.b){
    min.cs <- min(cohort.size.combi.b)
    n.occ <- length(which(cohort.queue==6))
    if(!n.occ>0){
      stop("Trial `combi.b` is active (`active.combi.b`) but the `cohort.queue` does not contain cohorts for it.")
    }
    if(!min.cs*n.occ >= max.n.combi.b){
      message("NOTE:\n",
              "`cohort.queue` does not contain sufficiently many cohorts for trial `combi.b`\n",
              "to guarantee that `max.n.combi.b` can be reached (based on minimal cohort size).\n",
              "The pattern specified in the `cohort.queue` will be repeated until sufficiently\n",
              "many cohorts are included in the `cohort.queue`.\n",
              "To avoid this warning, supply a `cohort.queue` that contains as least ",
              paste0(ceiling(max.n.combi.b/min.cs)), "\n cohorts for `combi.b`.")
      cohort.queue <- rep(cohort.queue, times = ceiling(ceiling(max.n.combi.b/min.cs)/n.occ)+1)
    }
  }

  #dummy doses to avoid error at the end of the function.
  if(!active.combi.a){
    doses.combi.a <- rbind(c(1), c(1))
    tox.combi.a <- c(0.1)
    start.dose.combi.a1 <- 1
    start.dose.combi.a2 <- 1
  }else{
    if(!is.num(doses.combi.a, low=0, lB=F, uB=F)){
      stop("`doses.combi.a` must contain positive real numbers.")
    }
    if(!length(dim(doses.combi.a))==2){
      stop("`doses.combi.a` must be a matrix with 2 rows.")
    }
    if(!dim(doses.combi.a)[1]==2){
      stop("`doses.combi.a` must be a matrix with 2 rows.")
    }
    if(!is.num(tox.combi.a, low=0, up=1, lB=F, uB=F,
               len=length(doses.combi.a[1,]))){
      stop("`tox.combi.a` must contain for each column of `doses.combi.a` one positive\n",
           "real number strictly below 1.")
    }
    if(!length(dim(tox.combi.a))==0){
      if(!dim(tox.combi.a)[1]==1){
        stop("`tox.combi.a` must be a vector or matrix with 1 row.")
      }
    }
    if(!is.num(start.dose.combi.a1, len=1, low=0, lB=F, uB=F)){
      stop("`start.dose.combi.a1` must be a positive real number.")
    }
    if(!is.num(start.dose.combi.a2, len=1, low=0, lB=F, uB=F)){
      stop("`start.dose.combi.a2` must be a positive real number.")
    }
    if(!paste(start.dose.combi.a1, start.dose.combi.a2, sep="+") %in%
       paste(doses.combi.a[1, ], doses.combi.a[2, ], sep="+")){
      stop("`start.dose.combi.a1` and `start.dose.combi.a2` must specify a dose level\n",
           "that was given in `doses.combi.a`.")
    }
  }
  if(!active.combi.b){
    doses.combi.b <- rbind(c(1), c(1))
    tox.combi.b <- c(0.1)
    start.dose.combi.b1 <- 1
    start.dose.combi.b2 <- 1
  }else{
    if(!is.num(doses.combi.b, low=0, lB=F, uB=F)){
      stop("`doses.combi.b` must contain positive real numbers.")
    }
    if(!length(dim(doses.combi.b))==2){
      stop("`doses.combi.b`` must be a matrix with 2 rows.")
    }
    if(!dim(doses.combi.b)[1]==2){
      stop("`doses.combi.b` must be a matrix with 2 rows.")
    }
    if(!is.num(tox.combi.b, low=0, up=1, lB=F, uB=F,
               len=length(doses.combi.b[1,]))){
      stop("`tox.combi.b` must contain for each column of `doses.combi.b` one positive\n",
           "real number strictly below 1.")
    }
    if(!length(dim(tox.combi.b))==0){
      if(!dim(tox.combi.b)[1]==1){
        stop("`tox.combi.b` must be a vector or matrix with 1 row.")
      }
    }
    if(!is.num(start.dose.combi.b1, len=1, low=0, lB=F, uB=F)){
      stop("`start.dose.combi.b1` must be a positive real number.")
    }
    if(!is.num(start.dose.combi.b2, len=1, low=0, lB=F, uB=F)){
      stop("`start.dose.combi.b2` must be a positive real number.")
    }
    if(!paste(start.dose.combi.b1, start.dose.combi.b2, sep="+") %in%
       paste(doses.combi.b[1, ], doses.combi.b[2, ], sep="+")){
      stop("`start.dose.combi.b1` and `start.dose.combi.b2` must specify a dose level\n",
           "that was given in `doses.combi.b`.")
    }
  }
  if(!active.mono1.a){
    doses.mono1.a <- c(1)
    tox.mono1.a <- c(0.1)
    start.dose.mono1.a <- 1
  }else{
    if(!is.num(doses.mono1.a, low=0, lB=F, uB=F)){
      stop("`doses.mono1.a` must contain positive real numbers.")
    }
    if(!length(dim(doses.mono1.a))==0){
      if(!dim(doses.mono1.a)[1]==1){
        stop("`doses.mono1.a` must be a vector or matrix with 1 row.")
      }
    }
    if(!is.num(tox.mono1.a, low=0, up=1, lB=F, uB=F,
               len=length(doses.mono1.a))){
      stop("`tox.mono1.a` must contain for each dose in `doses.mono1.a` one positive\n",
           "real number strictly below 1.")
    }
    if(!length(dim(tox.mono1.a))==0 ){
      if(!dim(tox.mono1.a)[1]==1){
        stop("`tox.mono1.a` must be a vector or matrix with 1 row.")
      }
    }
    if(!is.num(start.dose.mono1.a, len=1, low=0, lB=F, uB=F)){
      stop("`start.dose.mono1.a` must be a positive real number.")
    }
    if(!start.dose.mono1.a%in%doses.mono1.a){
      stop("`start.dose.mono1.a` must be a dose in `doses.mono1.a`.")
    }
  }
  if(!active.mono1.b){
    doses.mono1.b <- c(1)
    tox.mono1.b <- c(0.1)
    start.dose.mono1.b <- 1
  }else{
    if(!is.num(doses.mono1.b, low=0, lB=F, uB=F)){
      stop("`doses.mono1.b` must contain positive real numbers.")
    }
    if(!length(dim(doses.mono1.b))==0){
      if(!dim(doses.mono1.b)[1]==1){
        stop("`doses.mono1.b` must be a vector or matrix with 1 row.")
      }
    }
    if(!is.num(tox.mono1.b, low=0, up=1, lB=F, uB=F,
               len=length(doses.mono1.b))){
      stop("`tox.mono1.b` must contain for each dose in `doses.mono1.b` one positive\n",
           "real number strictly below 1.")
    }
    if(!length(dim(tox.mono1.b))==0){
      if(!dim(tox.mono1.b)[1]==1){
        stop("`tox.mono1.b` must be a vector or matrix with 1 row.")
      }
    }
    if(!is.num(start.dose.mono1.b, len=1, low=0, lB=F, uB=F)){
      stop("`start.dose.mono1.b` must be a positive real number.")
    }
    if(!start.dose.mono1.b%in%doses.mono1.b){
      stop("`start.dose.mono1.b` must be a dose in `doses.mono1.b`.")
    }
  }
  if(!active.mono2.a){
    doses.mono2.a <- c(1)
    tox.mono2.a <- c(0.1)
    start.dose.mono2.a <- 1
  }else{
    if(!is.num(doses.mono2.a, low=0, lB=F, uB=F)){
      stop("`doses.mono2.a` must contain positive real numbers.")
    }
    if(!length(dim(doses.mono2.a))==0){
      if(!dim(doses.mono2.a)[1]==1){
        stop("`doses.mono2.a` must be a vector or matrix with 1 row.")
      }
    }
    if(!is.num(tox.mono2.a, low=0, up=1, lB=F, uB=F,
               len=length(doses.mono2.a))){
      stop("`tox.mono2.a` must contain for each dose in `doses.mono2.a` one positive\n",
           "real number strictly below 1.")
    }
    if(!length(dim(tox.mono2.a))==0 ){
      if(!dim(tox.mono2.a)[1]==1){
        stop("`tox.mono2.a` must be a vector or matrix with 1 row.")
      }
    }
    if(!is.num(start.dose.mono2.a, len=1, low=0, lB=F, uB=F)){
      stop("`start.dose.mono2.a` must be a positive real number.")
    }
    if(!start.dose.mono2.a%in%doses.mono2.a){
      stop("`start.dose.mono2.a` must be a dose in `doses.mono2.a`.")
    }
  }
  if(!active.mono2.b){
    doses.mono2.b <- c(1)
    tox.mono2.b <- c(0.1)
    start.dose.mono2.b <- 1
  }else{
    if(!is.num(doses.mono2.b, low=0, lB=F, uB=F)){
      stop("`doses.mono2.b` must contain positive real numbers.")
    }
    if(!length(dim(doses.mono2.b))==0){
      if(!dim(doses.mono2.b)[1]==1){
        stop("`doses.mono2.b` must be a vector or matrix with 1 row.")
      }
    }
    if(!is.num(tox.mono2.b, low=0, up=1, lB=F, uB=F,
               len=length(doses.mono2.b))){
      stop("`tox.mono2.b` must contain for each dose in `doses.mono2.b` one positive\n",
           "real number strictly below 1.")
    }
    if(!length(dim(tox.mono2.b))==0){
      if(!dim(tox.mono2.b)[1]==1){
        stop("`tox.mono2.b` must be a vector or matrix with 1 row.")
      }
    }
    if(!is.num(start.dose.mono2.b, len=1, low=0, lB=F, uB=F)){
      stop("`start.dose.mono2.b` must be a positive real number.")
    }
    if(!start.dose.mono2.b%in%doses.mono2.b){
      stop("`start.dose.mono2.b` must be a dose in `doses.mono2.b`.")
    }
  }

  if(!is.num(dose.ref1, len=1, low=0, lB=F, uB=F)){
    stop("`dose.ref1` must be a positive real number.")
  }

  if(!is.num(dose.ref2, len=1, low=0, lB=F, uB=F)){
    stop("`dose.ref2` must be a positive real number.")
  }


  #checks and consistency of escalation steps
  if(is.null(esc.step.mono1.a) & active.mono1.a){
    esc.step.mono1.a <- max_step_BLRM(doses.mono1.a)
  }

  if(is.null(esc.step.mono1.b) & active.mono1.b){
    esc.step.mono1.b <- max_step_BLRM(doses.mono1.b)
  }

  if(is.null(esc.step.mono2.a) & active.mono2.a){
    esc.step.mono2.a <- max_step_BLRM(doses.mono2.a)
  }

  if(is.null(esc.step.mono2.b) & active.mono2.b){
    esc.step.mono2.b <- max_step_BLRM(doses.mono2.b)
  }

  if(is.null(esc.step.combi.a1) & active.combi.a){
    esc.step.combi.a1 <- max_step_BLRM(doses.combi.a[1, ])
  }

  if(is.null(esc.step.combi.b1) & active.combi.b){
    esc.step.combi.b1 <- max_step_BLRM(doses.combi.b[1, ])
  }

  if(is.null(esc.step.combi.a2) & active.combi.a){
    esc.step.combi.a2 <- max_step_BLRM(doses.combi.a[2, ])
  }

  if(is.null(esc.step.combi.b2) & active.combi.b){
    esc.step.combi.b2 <- max_step_BLRM(doses.combi.b[2, ])
  }

  if(!is.num(esc.step, len=1,
             low=0, lB=F) & !is.null(esc.step)){
    stop("`esc.step` must be NULL or a positive number.")
  }

  if(active.combi.a & !is.num(esc.step.combi.a1, len=1,
             low=0, lB=F)){
    stop("`esc.step.combi.a1` must be a positive number or NULL.")
  }

  if(active.combi.a & !is.num(esc.step.combi.a2, len=1,
             low=0, lB=F)){
    stop("`esc.step.combi.a2` must be a positive number or NULL.")
  }

  if(active.combi.b & !is.num(esc.step.combi.b1, len=1,
             low=0, lB=F)){
    stop("`esc.step.combi.b1` must be a positive number or NULL.")
  }

  if(active.combi.b & !is.num(esc.step.combi.b2, len=1,
             low=0, lB=F)){
    stop("`esc.step.combi.b2` must be a positive number or NULL.")
  }

  if(active.mono1.a & !is.num(esc.step.mono1.a, len=1,
             low=0, lB=F)){
    stop("`esc.step.mono1.a` must be a positive number or NULL.")
  }

  if(active.mono1.b &!is.num(esc.step.mono1.b, len=1,
             low=0, lB=F)){
    stop("`esc.step.mono1.b` must be a positive number or NULL.")
  }

  if(active.mono2.a &!is.num(esc.step.mono2.a, len=1,
             low=0, lB=F)){
    stop("`esc.step.mono2.a` must be a positive number or NULL.")
  }

  if(active.mono2.b &!is.num(esc.step.mono2.b, len=1,
             low=0, lB=F)){
    stop("`esc.step.mono2.b` must be a positive number or NULL.")
  }

  if(active.mono1.a){
    if(!esc.step.mono1.a>=max_step_BLRM(doses.mono1.a)){
    stop("`esc.step.mono1.a` does not allow to reach all doses given in `doses.mono1.a`.\n",
         "Please specify a larger value for this argument.")
    }
  }

  if(active.mono1.b){
    if(!esc.step.mono1.b>=max_step_BLRM(doses.mono1.b)){
    stop("`esc.step.mono1.b` does not allow to reach all doses given in `doses.mono1.b`.\n",
         "Please specify a larger value for this argument.")

    }
  }

  if(active.mono2.a){
    if(!esc.step.mono2.a>=max_step_BLRM(doses.mono2.a)){
    stop("`esc.step.mono2.a` does not allow to reach all doses given in `doses.mono2.a`.\n",
         "Please specify a larger value for this argument.")
    }
  }

  if(active.mono2.b){
    if(!esc.step.mono2.b>=max_step_BLRM(doses.mono2.b)){
    stop("`esc.step.mono2.b` does not allow to reach all doses given in `doses.mono2.b`.\n",
         "Please specify a larger value for this argument.")
    }
  }

  if(active.combi.a){
    if(!esc.step.combi.a1>=max_step_BLRM(doses.combi.a[1,])){
    stop("`esc.step.combi.a1` does not allow to reach all doses given in `doses.combi.a[1,]`.\n",
         "Please specify a larger value for this argument.")
    }
  }

  if(active.combi.b){
    if(!esc.step.combi.b1>=max_step_BLRM(doses.combi.b[1,])){
    stop("`esc.step.combi.b1` does not allow to reach all doses given in `doses.combi.b[1,]`.\n",
         "Please specify a larger value for this argument.")
    }
  }

  if(active.combi.a){
    if(!esc.step.combi.a2>=max_step_BLRM(doses.combi.a[2,])){
    stop("`esc.step.combi.a2` does not allow to reach all doses given in `doses.combi.a[2,]`.\n",
         "Please specify a larger value for this argument.")
    }
  }

  if(active.combi.b){
    if(!esc.step.combi.b2>=max_step_BLRM(doses.combi.b[2,])){
    stop("`esc.step.combi.b2` does not allow to reach all doses given in `doses.combi.b[2,]`.\n",
         "Please specify a larger value for this argument.")
    }
  }

  #-------------------------------
  #Loss weights
  #-------------------------------
  if(esc.rule%in%c("loss")){
    if(!is.num(loss.weights, len=4, lB=F, uB=T)){
      stop("`loss.weights` must be a vector of real numbers and have length 4.")
    }
  }
  if(esc.rule%in%c("dynamic.loss", "dynamic")){
    if(!is.num(dynamic.weights, lB=F, uB=F)){
      stop("`dynamic.weights` must contain real numbers.")
    }
    if(!(all(dim(dynamic.weights)==c(4, 4)))){
      stop("`dynamic.weights` must be a 4x4 matrix.")
    }
    norm_l1 <- sum(abs(dynamic.weights[1,]))
    if(!sum(abs(dynamic.weights[2,]))==norm_l1){
      stop("`dynamic.weights` must be normalized: the sum of the abolute values of\n",
           "each row must be the same for each row (failed in row 2).")
    }

    if(!sum(abs(dynamic.weights[3,]))==norm_l1){
      stop("`dynamic.weights` must be normalized: the sum of the abolute values of\n",
           "each row must be the same for each row (failed in row 3).")
    }

    if(!sum(abs(dynamic.weights[4,]))==norm_l1){
      stop("`dynamic.weights` must be normalized: the sum of the abolute values of\n",
           "each row must be the same for each row (failed in row 4).")
    }
  }


  if(esc.rule%in%c("ewoc", "ewoc.opt", "ewoc.max")){
    if(!is.num(ewoc.threshold, len=1, low=0, up=1, lB=F, uB=F)){
      stop("`ewoc.threshold` must be a number between 0 and 1.")
    }
  }

  #---------------------
  #check decision rules
  #---------------------
  if(!is.list(mtd.decision)){
    stop("`mtd.decision` must be a named list.")
  }
  if(is.null(names(mtd.decision))){
    stop("`mtd.decision` must be a named list.")
  }
  names(mtd.decision) <- toupper(names(mtd.decision))
  if(!is.dec.rule(mtd.decision)){
    stop("`mtd.decision` does not follow the format specified in the documentation.")
  }


  if(!is.list(mtd.decision.combi.a)){
    stop("`mtd.decision.combi.a` must be a named list.")
  }
  if(is.null(names(mtd.decision.combi.a))){
    stop("`mtd.decision.combi.a` must be a named list.")
  }
  names(mtd.decision.combi.a) <- toupper(names(mtd.decision.combi.a))
  if(!is.dec.rule(mtd.decision.combi.a)){
    stop("`mtd.decision.combi.a` does not follow the format specified in the documentation.")
  }

  if(!is.list(mtd.decision.combi.b)){
    stop("`mtd.decision.combi.b` must be a named list.")
  }
  if(is.null(names(mtd.decision.combi.b))){
    stop("`mtd.decision.combi.b` must be a named list.")
  }
  names(mtd.decision.combi.b) <- toupper(names(mtd.decision.combi.b))
  if(!is.dec.rule(mtd.decision.combi.b)){
    stop("`mtd.decision.combi.b` does not follow the format specified in the documentation.")
  }

  if(!is.list(mtd.decision.mono1.a)){
    stop("`mtd.decision.mono1.a` must be a named list.")
  }
  if(is.null(names(mtd.decision.mono1.a))){
    stop("`mtd.decision.mono1.a` must be a named list.")
  }
  names(mtd.decision.mono1.a) <- toupper(names(mtd.decision.mono1.a))
  if(!is.dec.rule(mtd.decision.mono1.a)){
    stop("`mtd.decision.mono1.a` does not follow the format specified in the documentation.")
  }

  if(!is.list(mtd.decision.mono1.b)){
    stop("`mtd.decision.mono1.b` must be a named list.")
  }
  if(is.null(names(mtd.decision.mono1.b))){
    stop("`mtd.decision.mono1.b` must be a named list.")
  }
  names(mtd.decision.mono1.b) <- toupper(names(mtd.decision.mono1.b))
  if(!is.dec.rule(mtd.decision.mono1.b)){
    stop("`mtd.decision.mono1.b` does not follow the format specified in the documentation.")
  }

  if(!is.list(mtd.decision.mono2.a)){
    stop("`mtd.decision.mono2.a` must be a named list.")
  }
  if(is.null(names(mtd.decision.mono2.a))){
    stop("`mtd.decision.mono2.a` must be a named list.")
  }
  names(mtd.decision.mono2.a) <- toupper(names(mtd.decision.mono2.a))
  if(!is.dec.rule(mtd.decision.mono2.a)){
    stop("`mtd.decision.mono2.a` does not follow the format specified in the documentation.")
  }

  if(!is.list(mtd.decision.mono2.b)){
    stop("`mtd.decision.mono2.b` must be a named list.")
  }
  if(is.null(names(mtd.decision.mono2.b))){
    stop("`mtd.decision.mono2.b` must be a named list.")
  }
  names(mtd.decision.mono2.b) <- toupper(names(mtd.decision.mono2.b))
  if(!is.dec.rule(mtd.decision.mono2.b)){
    stop("`mtd.decision.mono2.b` does not follow the format specified in the documentation.")
  }

  #----------------------
  #check priors
  #----------------------
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

  #--------------------------------
  #check backfill-starting doses
  #--------------------------------
  #set backfill start to first available dose if not given
  if(is.null(backfill.start.mono1.a)){
    backfill.start.mono1.a <- doses.mono1.a[1]
  }

  if(is.null(backfill.start.mono1.b)){
    backfill.start.mono1.b <- doses.mono1.b[1]
  }

  if(is.null(backfill.start.mono2.a)){
    backfill.start.mono2.a <- doses.mono2.a[1]
  }

  if(is.null(backfill.start.mono2.b)){
    backfill.start.mono2.b <- doses.mono2.b[1]
  }

  if(is.null(backfill.start.combi.a1)){
    backfill.start.combi.a1 <- doses.combi.a[1, 1]
  }

  if(is.null(backfill.start.combi.b1)){
    backfill.start.combi.b1 <- doses.combi.b[1, 1]
  }
  if(is.null(backfill.start.combi.a2)){
    backfill.start.combi.a2 <- doses.combi.a[2, 1]
  }

  if(is.null(backfill.start.combi.b2)){
    backfill.start.combi.b2 <- doses.combi.b[2, 1]
  }


  #checks of backfill starting doses
  if(backfill.mono1.a & active.mono1.a){
    if(!is.numeric(backfill.start.mono1.a)){
      stop("`backfill.start.mono1.a` must be numeric")
    }
    if(!length(backfill.start.mono1.a)==1){
      stop("`backfill.start.mono1.a` must have length 1")
    }
    if(!backfill.start.mono1.a %in% doses.mono1.a){
      stop("`backfill.start.mono1.a` must be an element of `doses.mono1.a`")
    }
  }

  if(backfill.mono1.b & active.mono1.b){
    if(!is.numeric(backfill.start.mono1.b)){
      stop("`backfill.start.mono1.b` must be numeric")
    }
    if(!length(backfill.start.mono1.b)==1){
      stop("`backfill.start.mono1.b` must have length 1")
    }
    if(!backfill.start.mono1.b %in% doses.mono1.b){
      stop("`backfill.start.mono1.b` must be an element of `doses.mono1.b`")
    }
  }

  if(backfill.mono2.a & active.mono2.a){
    if(!is.numeric(backfill.start.mono2.a)){
      stop("`backfill.start.mono2.a` must be numeric")
    }
    if(!length(backfill.start.mono2.a)==1){
      stop("`backfill.start.mono2.a` must have length 1")
    }
    if(!backfill.start.mono2.a %in% doses.mono2.a){
      stop("`backfill.start.mono2.a` must be an element of `doses.mono2.a`")
    }
  }

  if(backfill.mono2.b & active.mono2.b){
    if(!is.numeric(backfill.start.mono2.b)){
      stop("`backfill.start.mono2.b` must be numeric")
    }
    if(!length(backfill.start.mono2.b)==1){
      stop("`backfill.start.mono2.b` must have length 1")
    }
    if(!backfill.start.mono2.b %in% doses.mono2.b){
      stop("`backfill.start.mono2.b` must be an element of `doses.mono2.b`")
    }
  }

  if(backfill.combi.a & active.combi.a){
    if(!is.numeric(backfill.start.combi.a1)){
      stop("`backfill.start.combi.a1` must be numeric")
    }
    if(!length(backfill.start.combi.a1)==1){
      stop("`backfill.start.combi.a1` must have length 1")
    }
    if(!is.numeric(backfill.start.combi.a2)){
      stop("`backfill.start.combi.a2` must be numeric")
    }
    if(!length(backfill.start.combi.a2)==1){
      stop("`backfill.start.combi.a2` must have length 1")
    }
    if(!length(which(
      doses.combi.a[1, ]== backfill.start.combi.a1 &
      doses.combi.a[2, ]== backfill.start.combi.a2
    ))>0){
      stop("`backfill.start.combi.a1` must be a dose given in `doses.combi.a`")
    }
  }

  if(backfill.combi.b & active.combi.b){
    if(!is.numeric(backfill.start.combi.b1)){
      stop("`backfill.start.combi.b1` must be numeric")
    }
    if(!length(backfill.start.combi.b1)==1){
      stop("`backfill.start.combi.b1` must have length 1")
    }
    if(!is.numeric(backfill.start.combi.b2)){
      stop("`backfill.start.combi.b2` must be numeric")
    }
    if(!length(backfill.start.combi.b2)==1){
      stop("`backfill.start.combi.b2` must have length 1")
    }
    if(!length(which(
      doses.combi.b[1, ]== backfill.start.combi.b1 &
      doses.combi.b[2, ]== backfill.start.combi.b2
    ))>0){
      stop("`backfill.start.combi.b1` must be a dose given in `doses.combi.b`")
    }
  }

  #------------------------------------------------------
  #Check historical data dyntax and consistency
  #------------------------------------------------------

  if(!is.list(historical.data) & !is.null(historical.data)){
    stop("`historical.data` must be NULL or a named list.")
  }
  if(!is.null(historical.data)){
    if(is.null(names(historical.data))){
      stop("`historical.data` must be NULL or a named list.")
    }
    names(historical.data) <- tolower(names(historical.data))
    if(!is.historical.data(historical.data)){
      stop("`historical.data` does not follow the format specified in the documentation.")
    }
    if(!all(historical.data$n.dlt <= historical.data$n.pat)){
      stop("`historical.data` may not contain observations with more DLT than patients.")
    }
    historical.data$trial <- tolower(historical.data$trial)

    #strip all NA doses, and detect which cohorts can be removed
    historical.data <- clean.na.hist(historical.data)
    idx_noninf_obs <- which((historical.data$dose1==0 &
                               historical.data$dose2==0) |
                              (historical.data$n.pat==0)
    )
    remove.hist <- FALSE
    if(!length(idx_noninf_obs)==0){
      message("NOTE: `historical.data` cohorts were detected for which both doses are 0 and/or \n",
              " no patients were treated. These observations were removed.")
      if(length(idx_noninf_obs)==length(historical.data$dose1)){
        remove.hist <- TRUE
      }else{
        historical.data <- remove.noninf.obs(historical.data, idx_noninf_obs)
      }
    }
    #continue consistency checks
    if(!remove.hist){
      stds <- historical.data$trial
      idx_m1 <- which(stds %in% c(1, "mono1.a"))
      if(!length(idx_m1)==0){
        if(!active.mono1.a){
          stop("`mono1.a` is not active, but `historical.data` is given for this trial.")
        }
        if(any(historical.data$dose2[idx_m1]>0)){
          stop("`historical.data` for `mono1.a`/trial 1 must be mono 1 observations.")
        }
      }
      idx_m4 <- which(stds %in% c(4, "mono1.b"))
      if(!length(idx_m4)==0){
        if(!active.mono1.b){
          stop("`mono1.b` is not active, but `historical.data` is given for this trial.")
        }
        if(any(historical.data$dose2[idx_m4]>0)){
          stop("`historical.data` for `mono1.b`/trial 4 must be mono 1 observations.")
        }
      }
      idx_m2 <- which(stds %in% c(2, "mono2.a"))
      if(!length(idx_m2)==0){
        if(!active.mono2.a){
          stop("`mono2.a` is not active, but `historical.data` is given for this trial.")
        }
        if(any(historical.data$dose1[idx_m2]>0)){
          stop("`historical.data` for `mono2.a`/trial 2 must be mono 2 observations.")
        }
      }
      idx_m5 <- which(stds %in% c(5, "mono2.b"))
      if(!length(idx_m5)==0){
        if(!active.mono2.b){
          stop("`mono2.b` is not active, but `historical.data` is given for this trial.")
        }
        if(any(historical.data$dose1[idx_m2]>0)){
          stop("`historical.data` for `mono2.b`/trial 5 must be mono 2 observations.")
        }
      }
      idx_c1 <- which(stds %in% c(3, "combi.a"))
      if(!length(idx_c1)==0){
        if(!active.combi.a){
          stop("`combi.a` is not active, but `historical.data` is given for this trial.")
        }
        if(any(historical.data$dose1[idx_c1]==0 | historical.data$dose2[idx_c1]==0)){
          stop("`historical.data` for `combi.a`/trial 3 must be combination observations.")
        }
      }
      idx_c2 <- which(stds %in% c(6, "combi.b"))
      if(!length(idx_c2)==0){
        if(!active.combi.b){
          stop("`combi.b` is not active, but `historical.data` is given for this trial.")
        }
        if(any(historical.data$dose1[idx_c2]==0 | historical.data$dose2[idx_c2]==0)){
          stop("`historical.data` for `combi.b`/trial 6 must be combination observations.")
        }
      }
      historical.data$trial <- harmonize_vecnames_jointBLRM(historical.data$trial)
      fx <- levels(factor(historical.data$trial))
      idx_hist <- which(!fx%in%c(1, 2, 3, 4, 5, 6))
      if(!length(idx_hist)==0){
        for(n in 1:length(idx_hist)){
          idx_curr_hstd <- which(historical.data$trial%in%c(fx[idx_hist[n]]))
          if(!(all(historical.data$dose1[idx_curr_hstd]==0 & historical.data$dose2[idx_curr_hstd]>0) |
               all(historical.data$dose2[idx_curr_hstd]==0 & historical.data$dose1[idx_curr_hstd]>0) |
               all(historical.data$dose1[idx_curr_hstd]>0 & historical.data$dose2[idx_curr_hstd]>0))){
            stop("A trial in `historical.data` contains observations of more than one trial type (mono 1,\n",
                 "mono 2, combination). This is not permitted.")
          }
        }
      }

    }else{
      historical.data <- NULL
    }
  }

  #check whether dose levels are sorted
  if(is.unsorted(doses.mono1.a, strictly = TRUE)){
    stop("`doses.mono1.a` must be sorted (have strictly ascending entries).")
  }
  if(is.unsorted(doses.mono1.b, strictly = TRUE)){
    stop("`doses.mono1.b` must be sorted (have strictly ascending entries).")
  }
  if(is.unsorted(doses.mono2.a, strictly = TRUE)){
    stop("`doses.mono2.a` must be sorted (have strictly ascending entries).")
  }
  if(is.unsorted(doses.mono2.b, strictly = TRUE)){
    stop("`doses.mono2.b` must be sorted (have strictly ascending entries).")
  }


  #-----------
  #Working Dir
  if(!is.null(working.path) & !is.null(file.name)){
    if(!is.character(working.path)){
      stop("`working.path` must be of character type")
    }else{
      if(!dir.exists(working.path)){
        stop("`working.path` was supplied but does not point to a directory")
      }
      if(clean.working.path){
        #check for temp files using naming conventions and delete automatically
         if(length(list.files(working.path, full.names = T, pattern = paste0("^", file.name, "_tmp.")))>0){
          message("Cleaning existing temporary data from working directory.")
          #message(paste0("Removing files starting with '", paste0(file.name, "_tmp"), "' from working directory."))
          #initial cleanup of working directory
          file.remove(list.files(working.path, full.names = T, pattern = paste0("^", file.name, "_tmp.")))
        }

      }else{
        #check for temp files and ask user to clean directory manually if detected
        if(length(list.files(working.path, full.names = T, pattern = paste0("^", file.name, "_tmp.")))>0){
          stop(paste0("`working.path` specifies a working directory that contains files that include `",
                      paste0(file.name, "_tmp"),"` in the file name. Such files may be accessed, modified, or deleted by `sim_jointBLRM()`. ",
                      "Please supply a working directory that does not contain such files or provide the argument ",
                      "'clean.working.path' to allow the function to remove all such files prior to simulating."))
        }
      }
    }
  }

  #-----------------------------------------------------------------------------------------------------------------
  #Checks done, start the simulation
  #-----------------------------------------------------------------------------------------------------------------

  set.seed(seed=seed)

  dose.name.vector.combi1 <- paste(doses.combi.a[1,],doses.combi.a[2,], sep="+")
  dose.name.vector.combi2 <- paste(doses.combi.b[1,],doses.combi.b[2,], sep="+")

  #names for the toxicities,in order to be able to read out the values for a given dose
  names(tox.combi.a) <- dose.name.vector.combi1
  names(tox.mono1.a)<- doses.mono1.a
  names(tox.mono2.a)<- doses.mono2.a
  names(tox.combi.b) <- dose.name.vector.combi2
  names(tox.mono1.b)<- doses.mono1.b
  names(tox.mono2.b)<- doses.mono2.b


  #sample seeds for each trial
  trial_seeds = sample.int(.Machine$integer.max, n.studies)
  #Note: this ensures that the trials will be based on the same seeds, regardless
  #of whether the simulation is parallelized or not.

  #parallelization
  if(n.cores>1){
    cl_foreach <- makeCluster(n.cores)
    registerDoParallel(cl_foreach)
    message("Parallel simulation with ", n.cores, " cores started.")
  }else{
    registerDoSEQ()
    message("Sequential simulation started.")
  }


  # ----------------------------------------------------------------------------
  # ACTUAL SIMULATION
  # ----------------------------------------------------------------------------
  #foreach loop over the desired number of trials.
     res.list <- foreach( i = 1:n.studies,
                          .packages = c("OncoBLRM"),
                          #.packages = c("openxlsx", "flock", "rstan", "openssl"),
                          # .export = c("get_int_probs_BLRM", "inv_logit", "logit",
                          #             "main_jointBLRM", "post_tox_jointBLRM_sim",
                          #             "trial_jointBLRM", "harmonize_vecnames_jointBLRM",
                          #             "get_dloss_dec_jointBLRM", "get_loss_dec_jointBLRM",
                          #             "get_ewoc_dec_jointBLRM"),
                          .errorhandling = "pass",
                          .inorder = FALSE)%dopar%
  # res.list <- list()
  # for(i in 1:n.studies)
      {

        res<- OncoBLRM:::trial_jointBLRM(  doses.mono1.a = doses.mono1.a,
              #OncoBLRM:::trial_jointBLRM(  doses.mono1.a = doses.mono1.a,
                                doses.mono2.a = doses.mono2.a,
                                doses.combi.a = doses.combi.a,
                                doses.mono1.b = doses.mono1.b,
                                doses.mono2.b = doses.mono2.b,
                                doses.combi.b = doses.combi.b,
                                start.dose.mono1.a= start.dose.mono1.a,
                                start.dose.mono2.a= start.dose.mono2.a,
                                start.dose.mono1.b= start.dose.mono1.b,
                                start.dose.mono2.b= start.dose.mono2.b,
                                start.dose.combi.a1= start.dose.combi.a1,
                                start.dose.combi.a2= start.dose.combi.a2,
                                start.dose.combi.b1 = start.dose.combi.b1,
                                start.dose.combi.b2 = start.dose.combi.b2,
                                seed = trial_seeds[i],
              #BLRM = BLRM,

                                historical.data= historical.data,

                                active.mono1.a = active.mono1.a,
                                active.mono1.b = active.mono1.b,
                                active.mono2.a = active.mono2.a,
                                active.mono2.b = active.mono2.b,
                                active.combi.a = active.combi.a,
                                active.combi.b = active.combi.b,

                                cohort.queue = cohort.queue,

                                tox.combi.a = tox.combi.a,
                                tox.mono1.a = tox.mono1.a,
                                tox.mono2.a = tox.mono2.a,
                                tox.combi.b = tox.combi.b,
                                tox.mono1.b = tox.mono1.b,
                                tox.mono2.b = tox.mono2.b,

                                cohort.size.mono1.a= cohort.size.mono1.a,
                                cohort.prob.mono1.a= cohort.prob.mono1.a,
                                cohort.size.mono2.a= cohort.size.mono2.a,
                                cohort.prob.mono2.a= cohort.prob.mono2.a,
                                cohort.size.mono1.b= cohort.size.mono1.b,
                                cohort.prob.mono1.b= cohort.prob.mono1.b,
                                cohort.size.mono2.b= cohort.size.mono2.b,
                                cohort.prob.mono2.b= cohort.prob.mono2.b,
                                cohort.size.combi.a= cohort.size.combi.a,
                                cohort.prob.combi.a= cohort.prob.combi.a,
                                cohort.size.combi.b= cohort.size.combi.b,
                                cohort.prob.combi.b= cohort.prob.combi.b,

                                esc.rule = esc.rule,
                                esc.comp.max = esc.comp.max,

                                esc.step.mono1.a=esc.step.mono1.a,
                                esc.step.mono2.a=esc.step.mono2.a,
                                esc.step.mono1.b=esc.step.mono1.b,
                                esc.step.mono2.b=esc.step.mono2.b,
                                esc.step.combi.a1=esc.step.combi.a1,
                                esc.step.combi.b1=esc.step.combi.b1,
                                esc.step.combi.a2=esc.step.combi.a2,
                                esc.step.combi.b2=esc.step.combi.b2,

                                esc.constrain.mono1.a=esc.constrain.mono1.a,
                                esc.constrain.mono2.a=esc.constrain.mono2.a,
                                esc.constrain.mono1.b=esc.constrain.mono1.b,
                                esc.constrain.mono2.b=esc.constrain.mono2.b,
                                esc.constrain.combi.a1=esc.constrain.combi.a1,
                                esc.constrain.combi.b1=esc.constrain.combi.b1,
                                esc.constrain.combi.a2=esc.constrain.combi.a2,
                                esc.constrain.combi.b2=esc.constrain.combi.b2,

                                prior.tau = prior.tau,
                                prior.mu = prior.mu,
                                saturating = saturating,

                                dose.ref1 = dose.ref1,
                                dose.ref2 = dose.ref2,

                                dosing.intervals = dosing.intervals,
                                loss.weights = loss.weights,
                                dynamic.weights = dynamic.weights,
                                ewoc = ewoc.threshold,

                                max.n.mono1.a=max.n.mono1.a,
                                max.n.mono1.b=max.n.mono1.b,
                                max.n.mono2.a=max.n.mono2.a,
                                max.n.mono2.b=max.n.mono2.b,
                                max.n.combi.a=max.n.combi.a,
                                max.n.combi.b=max.n.combi.b,
                                decision.combi.a=mtd.decision.combi.a,
                                decision.combi.b=mtd.decision.combi.b,
                                decision.mono1.a=mtd.decision.mono1.a,
                                decision.mono1.b=mtd.decision.mono1.b,
                                decision.mono2.a=mtd.decision.mono2.a,
                                decision.mono2.b=mtd.decision.mono2.b,

                                mtd.enforce.mono1.a = mtd.enforce.mono1.a,
                                mtd.enforce.mono1.b = mtd.enforce.mono1.b,
                                mtd.enforce.mono2.a = mtd.enforce.mono2.a,
                                mtd.enforce.mono2.b = mtd.enforce.mono2.b,
                                mtd.enforce.combi.a = mtd.enforce.combi.a,
                                mtd.enforce.combi.b = mtd.enforce.combi.b,

                                backfill.mono1.a = backfill.mono1.a,
                                backfill.mono1.b = backfill.mono1.b,
                                backfill.size.mono1.a = backfill.size.mono1.a,
                                backfill.size.mono1.b = backfill.size.mono1.b,
                                backfill.prob.mono1.a = backfill.prob.mono1.a,
                                backfill.prob.mono1.b = backfill.prob.mono1.b,
                                backfill.mono2.a = backfill.mono2.a,
                                backfill.mono2.b = backfill.mono2.b,
                                backfill.size.mono2.a = backfill.size.mono2.a,
                                backfill.size.mono2.b = backfill.size.mono2.b,
                                backfill.prob.mono2.a = backfill.prob.mono2.a,
                                backfill.prob.mono2.b = backfill.prob.mono2.b,
                                backfill.combi.a = backfill.combi.a,
                                backfill.combi.b = backfill.combi.b,
                                backfill.size.combi.a = backfill.size.combi.a,
                                backfill.size.combi.b = backfill.size.combi.b,
                                backfill.prob.combi.a = backfill.prob.combi.a,
                                backfill.prob.combi.b = backfill.prob.combi.b,
                                backfill.start.mono1.a = backfill.start.mono1.a,
                                backfill.start.mono1.b = backfill.start.mono1.b,
                                backfill.start.mono2.a = backfill.start.mono2.a,
                                backfill.start.mono2.b = backfill.start.mono2.b,
                                backfill.start.combi.a1 = backfill.start.combi.a1,
                                backfill.start.combi.a2 = backfill.start.combi.a2,
                                backfill.start.combi.b1 = backfill.start.combi.b1,
                                backfill.start.combi.b2 = backfill.start.combi.b2,

                                working.path = working.path,
                                file.name = file.name,
                                chains = chains,
                                iter = iter,
                                refresh = refresh,
                                adapt_delta = adapt_delta,
                                warmup = warmup,
                                max_treedepth = max_treedepth
      )

      if(!is.null(file.name)&!is.null(monitor.path)){
        if(dir.exists(file.path(monitor.path))){
          write.table(matrix("Done", nrow=1, ncol=1),
                     file=paste(file.path(monitor.path),
                                "/trial-",i, "_", file.name,".txt", sep=""),
                     row.names = F, col.names = F, append = FALSE)
        }
      }

      return(res)
    }


  #clusters need to be stopped
  if(n.cores>1){stopCluster(cl_foreach)}
  #---------------------------------------------------------------------------
  #Post-Processing the results
  #---------------------------------------------------------------------------

  #--------------------------
  #setup output structures.
  #--------------------------
  #for number of under/target/over MTDs and stopped trials
  results.combi1 <- matrix(0, nrow=3, ncol=5)
  colnames(results.combi1) <- c('underdose', 'target dose', 'overdose', 'max n reached before MTD', 'all doses too toxic')
  rownames(results.combi1) <- c('Number of trials', 'Percentage', 'Percentage not all too toxic')
  results.combi1['Percentage not all too toxic', 'all doses too toxic']<- NA

  results.combi2 <- matrix(0, nrow=3, ncol=5)
  colnames(results.combi2) <- c('underdose', 'target dose', 'overdose', 'max n reached before MTD', 'all doses too toxic')
  rownames(results.combi2) <- c('Number of trials', 'Percentage', 'Percentage not all too toxic')
  results.combi2['Percentage not all too toxic', 'all doses too toxic']<- NA

  results.mono.1 <- matrix(0, nrow=3, ncol=5)
  colnames(results.mono.1) <- c('underdose', 'target dose', 'overdose', 'max n reached before MTD', 'all doses too toxic')
  rownames(results.mono.1) <- c('Number of trials', 'Percentage', 'Percentage not all too toxic')
  results.mono.1['Percentage not all too toxic', 'all doses too toxic']<- NA

  results.mono.2 <- matrix(0, nrow=3, ncol=5)
  colnames(results.mono.2) <- c('underdose', 'target dose', 'overdose', 'max n reached before MTD', 'all doses too toxic')
  rownames(results.mono.2) <- c('Number of trials', 'Percentage', 'Percentage not all too toxic')
  results.mono.2['Percentage not all too toxic', 'all doses too toxic']<- NA

  results.mono.4 <- matrix(0, nrow=3, ncol=5)
  colnames(results.mono.4) <- c('underdose', 'target dose', 'overdose', 'max n reached before MTD', 'all doses too toxic')
  rownames(results.mono.4) <- c('Number of trials', 'Percentage', 'Percentage not all too toxic')
  results.mono.4['Percentage not all too toxic', 'all doses too toxic']<- NA

  results.mono.5 <- matrix(0, nrow=3, ncol=5)
  colnames(results.mono.5) <- c('underdose', 'target dose', 'overdose', 'max n reached before MTD', 'all doses too toxic')
  rownames(results.mono.5) <- c('Number of trials', 'Percentage', 'Percentage not all too toxic')
  results.mono.5['Percentage not all too toxic', 'all doses too toxic']<- NA


  #Overview of treatments and number of patients.
  summary.combi1 <- matrix(NA, nrow=10, ncol=6)
  colnames(summary.combi1) <- c('Median', 'Mean', 'Min.', 'Max.', '2.5%', '97.5%')
  rownames(summary.combi1) <- c('#Pat underdose', '#Pat target dose', '#Pat overdose', '#Pat (all)', '#DLT underdose', '#DLT target dose', '#DLT overdose', '#DLT (all)', '% overdose', '% DLT')

  summary.combi2 <- matrix(NA, nrow=10, ncol=6)
  colnames(summary.combi2) <- c('Median', 'Mean', 'Min.', 'Max.', '2.5%', '97.5%')
  rownames(summary.combi2) <- c('#Pat underdose', '#Pat target dose', '#Pat overdose', '#Pat (all)', '#DLT underdose', '#DLT target dose', '#DLT overdose', '#DLT (all)', '% overdose', '% DLT')


  summary.mono.1 <- matrix(NA, nrow=10, ncol=6)
  colnames(summary.mono.1) <- c('Median', 'Mean', 'Min.', 'Max.', '2.5%', '97.5%')
  rownames(summary.mono.1) <- c('#Pat underdose', '#Pat target dose', '#Pat overdose', '#Pat (all)', '#DLT underdose', '#DLT target dose', '#DLT overdose', '#DLT (all)', '% overdose', '% DLT')

  summary.mono.2 <- matrix(NA, nrow=10, ncol=6)
  colnames(summary.mono.2) <- c('Median', 'Mean', 'Min.', 'Max.', '2.5%', '97.5%')
  rownames(summary.mono.2) <- c('#Pat underdose', '#Pat target dose', '#Pat overdose', '#Pat (all)', '#DLT underdose', '#DLT target dose', '#DLT overdose', '#DLT (all)', '% overdose', '% DLT')

  summary.mono.4 <- matrix(NA, nrow=10, ncol=6)
  colnames(summary.mono.4) <- c('Median', 'Mean', 'Min.', 'Max.', '2.5%', '97.5%')
  rownames(summary.mono.4) <- c('#Pat underdose', '#Pat target dose', '#Pat overdose', '#Pat (all)', '#DLT underdose', '#DLT target dose', '#DLT overdose', '#DLT (all)', '% overdose', '% DLT')

  summary.mono.5 <- matrix(NA, nrow=10, ncol=6)
  colnames(summary.mono.5) <- c('Median', 'Mean', 'Min.', 'Max.', '2.5%', '97.5%')
  rownames(summary.mono.5) <- c('#Pat underdose', '#Pat target dose', '#Pat overdose', '#Pat (all)', '#DLT underdose', '#DLT target dose', '#DLT overdose', '#DLT (all)', '% overdose', '% DLT')


  #Matrices that show MTD per dose and their input tox
  input.combi1 <- matrix(NA, ncol = length(tox.combi.a), nrow = 4)
  rownames(input.combi1) <- c('Dose', 'True P(DLT)', 'True category', 'MTD declared (n)')
  dose.name.vector.combi1 <- paste(doses.combi.a[1,],doses.combi.a[2,], sep="+")
  colnames(input.combi1) <- dose.name.vector.combi1

  input.combi2 <- matrix(NA, ncol = length(tox.combi.b), nrow = 4)
  rownames(input.combi2) <- c('Dose', 'True P(DLT)', 'True category', 'MTD declared (n)')
  dose.name.vector.combi2 <- paste(doses.combi.b[1,],doses.combi.b[2,], sep="+")
  colnames(input.combi2) <- dose.name.vector.combi2

  input.mono.1 <- matrix(NA, ncol = length(doses.mono1.a), nrow = 4)
  rownames(input.mono.1) <- c('Dose', 'True P(DLT)', 'True category', 'MTD declared (n)')
  colnames(input.mono.1) <- doses.mono1.a

  input.mono.2 <- matrix(NA, ncol = length(doses.mono2.a), nrow = 4)
  rownames(input.mono.2) <- c('Dose', 'True P(DLT)', 'True category', 'MTD declared (n)')
  colnames(input.mono.2) <- doses.mono2.a

  input.mono.4 <- matrix(NA, ncol = length(doses.mono1.b), nrow = 4)
  rownames(input.mono.4) <- c('Dose', 'True P(DLT)', 'True category', 'MTD declared (n)')
  colnames(input.mono.4) <- doses.mono1.b

  input.mono.5 <- matrix(NA, ncol = length(doses.mono2.b), nrow = 4)
  rownames(input.mono.5) <- c('Dose', 'True P(DLT)', 'True category', 'MTD declared (n)')
  colnames(input.mono.5) <- doses.mono2.b

  #----------------------------------------------------------------------------------------------------
  #already fill up parts of the input-data (toxicities, dose levels,....)

  input.combi1[1,] <- dose.name.vector.combi1
  input.combi1[2,] <- tox.combi.a
  input.combi1[4,] <- rep(0, times = length(tox.combi.a))

  input.combi2[1,] <- dose.name.vector.combi2
  input.combi2[2,] <- tox.combi.b
  input.combi2[4,] <- rep(0, times = length(tox.combi.b))

  input.mono.1[1,] <- doses.mono1.a
  input.mono.1[2,] <- tox.mono1.a
  input.mono.1[4,] <- rep(0, times = length(doses.mono1.a))

  input.mono.2[1,] <- doses.mono2.a
  input.mono.2[2,] <- tox.mono2.a
  input.mono.2[4,] <- rep(0, times = length(doses.mono2.a))

  input.mono.4[1,] <- doses.mono1.b
  input.mono.4[2,] <- tox.mono1.b
  input.mono.4[4,] <- rep(0, times = length(doses.mono1.b))

  input.mono.5[1,] <- doses.mono2.b
  input.mono.5[2,] <- tox.mono2.b
  input.mono.5[4,] <- rep(0, times = length(doses.mono2.b))

  #labels "under", "target", "over" are gvien to the doses. (For a "nicer" output format :) )
  for (l in 1:(length(tox.combi.a))){
    if (as.numeric(input.combi1[2,l])<dosing.intervals[1]){
      input.combi1[3,l] <- "underdose"
    } else if (dosing.intervals[1]<=as.numeric(input.combi1[2,l]) & as.numeric(input.combi1[2,l])<dosing.intervals[2]){
      input.combi1[3,l] <- "target dose"
    } else if (dosing.intervals[2]<=as.numeric(input.combi1[2,l])){
      input.combi1[3,l] <- "overdose"
    }
  }
  #...more labels
  for (l in 1:(length(tox.combi.b))){
    if (as.numeric(input.combi2[2,l])<dosing.intervals[1]){
      input.combi2[3,l] <- "underdose"
    } else if (dosing.intervals[1]<=as.numeric(input.combi2[2,l]) & as.numeric(input.combi2[2,l])<dosing.intervals[2]){
      input.combi2[3,l] <- "target dose"
    } else if (dosing.intervals[2]<=as.numeric(input.combi2[2,l])){
      input.combi2[3,l] <- "overdose"
    }
  }
  #...even more labels
  for (l in 1:(length(tox.mono1.a))){
    if (as.numeric(input.mono.1[2,l])<dosing.intervals[1]){
      input.mono.1[3,l] <- "underdose"
    } else if (dosing.intervals[1]<=as.numeric(input.mono.1[2,l]) & as.numeric(input.mono.1[2,l])<dosing.intervals[2]){
      input.mono.1[3,l] <- "target dose"
    } else if (dosing.intervals[2]<=as.numeric(input.mono.1[2,l])){
      input.mono.1[3,l] <- "overdose"
    }
  }
  #...labels, again
  for (l in 1:(length(tox.mono1.b))){
    if (as.numeric(input.mono.4[2,l])<dosing.intervals[1]){
      input.mono.4[3,l] <- "underdose"
    } else if (dosing.intervals[1]<=as.numeric(input.mono.4[2,l]) & as.numeric(input.mono.4[2,l])<dosing.intervals[2]){
      input.mono.4[3,l] <- "target dose"
    } else if (dosing.intervals[2]<=as.numeric(input.mono.4[2,l])){
      input.mono.4[3,l] <- "overdose"
    }
  }

  #...and labels
  for (l in 1:(length(tox.mono2.a))){
    if (as.numeric(input.mono.2[2,l])<dosing.intervals[1]){
      input.mono.2[3,l] <- "underdose"
    } else if (dosing.intervals[1]<=as.numeric(input.mono.2[2,l]) & as.numeric(input.mono.2[2,l])<dosing.intervals[2]){
      input.mono.2[3,l] <- "target dose"
    } else if (dosing.intervals[2]<=as.numeric(input.mono.2[2,l])){
      input.mono.2[3,l] <- "overdose"
    }
  }

  #...and labels
  for (l in 1:(length(tox.mono2.b))){
    if (as.numeric(input.mono.5[2,l])<dosing.intervals[1]){
      input.mono.5[3,l] <- "underdose"
    } else if (dosing.intervals[1]<=as.numeric(input.mono.5[2,l]) & as.numeric(input.mono.5[2,l])<dosing.intervals[2]){
      input.mono.5[3,l] <- "target dose"
    } else if (dosing.intervals[2]<=as.numeric(input.mono.5[2,l])){
      input.mono.5[3,l] <- "overdose"
    }
  }




  #----------------------------------------------------------
  #set up the other frequently used data structures
  #----------------------------------------------------------

  #vectors for saving the number of patients and dlts in each of the following simulated trials (these numbers will be needed to calculate the output later)
  #i.e. these vectors are just counters
  trials.n.all <- matrix(NA, nrow = 6, ncol = n.studies)
  trials.n.target <- matrix(NA, nrow = 6, ncol = n.studies)
  trials.n.under <- matrix(NA, nrow = 6, ncol = n.studies)
  trials.n.over <- matrix(NA, nrow = 6, ncol = n.studies)
  trials.dlt.all <- matrix(NA, nrow = 6, ncol = n.studies)
  trials.dlt.under <- matrix(NA, nrow = 6, ncol = n.studies)
  trials.dlt.target <- matrix(NA, nrow = 6, ncol = n.studies)
  trials.dlt.over <- matrix(NA, nrow = 6, ncol = n.studies)
  trials.perc.over <- matrix(NA, nrow = 6, ncol = n.studies)
  trials.perc.dlt <- matrix(NA, nrow = 6, ncol = n.studies)

  #also save number of patients at each dose for the trials.
  if(active.mono1.a){
    trials.p.dose.m1 <- matrix(NA, nrow = length(doses.mono1.a), ncol = n.studies)
    summ.ppd.m1a <- matrix(NA, nrow = 5, ncol= length(doses.mono1.a))
    summ.ppd.m1a[1,] <- paste0(doses.mono1.a)
    colnames(summ.ppd.m1a) <- paste0(doses.mono1.a)
    rownames(summ.ppd.m1a) <- c("Dose", "mean #pat", "median #pat", "min #pat", "max #pat")
  }
  if(active.mono2.a){
    trials.p.dose.m2 <- matrix(NA, nrow = length(doses.mono2.a), ncol = n.studies)
    summ.ppd.m2a <- matrix(NA, nrow = 5, ncol= length(doses.mono2.a))
    summ.ppd.m2a[1,] <- paste0(doses.mono2.a)
    colnames(summ.ppd.m2a) <- paste0(doses.mono2.a)
    rownames(summ.ppd.m2a) <- c("Dose", "mean #pat", "median #pat", "min #pat", "max #pat")
  }
  if(active.mono1.b){
    trials.p.dose.m4 <- matrix(NA, nrow = length(doses.mono1.b), ncol = n.studies)
    summ.ppd.m1b <- matrix(NA, nrow = 5, ncol= length(doses.mono1.b))
    summ.ppd.m1b[1,] <- paste0(doses.mono1.b)
    colnames(summ.ppd.m1b) <- paste0(doses.mono1.b)
    rownames(summ.ppd.m1b) <- c("Dose",  "mean #pat", "median #pat", "min #pat", "max #pat")
  }
  if(active.mono2.b){
    trials.p.dose.m5 <- matrix(NA, nrow = length(doses.mono2.b), ncol = n.studies)
    summ.ppd.m2b <- matrix(NA, nrow = 5, ncol= length(doses.mono2.b))
    summ.ppd.m2b[1,] <- paste0(doses.mono2.b)
    colnames(summ.ppd.m2b) <- paste0(doses.mono2.b)
    rownames(summ.ppd.m2b) <- c("Dose", "mean #pat", "median #pat", "min #pat", "max #pat")
  }
  if(active.combi.a){
    trials.p.dose.c1 <- matrix(NA, nrow = length(doses.combi.a[1,]), ncol = n.studies)
    summ.ppd.c1 <- matrix(NA, nrow = 5, ncol= length(doses.combi.a[1,]))
    summ.ppd.c1[1,] <- dose.name.vector.combi1
    colnames(summ.ppd.c1) <- dose.name.vector.combi1
    rownames(summ.ppd.c1) <- c("Dose", "mean #pat", "median #pat", "min #pat", "max #pat")
  }
  if(active.combi.b){
    trials.p.dose.c2 <- matrix(NA, nrow = length(doses.combi.b[1,]), ncol = n.studies)
    summ.ppd.c2 <- matrix(NA, nrow = 5, ncol= length(doses.combi.b[1,]))
    summ.ppd.c2[1,] <- dose.name.vector.combi2
    colnames(summ.ppd.c2) <- dose.name.vector.combi2
    rownames(summ.ppd.c2) <- c("Dose", "mean #pat", "median #pat", "min #pat", "max #pat")
  }


  #--------------------------------------------------------------------------------------------------------
  #Fill in output data in these stuctures
  #--------------------------------------------------------------------------------------------------------

  #extract MTDs, number of patients, etc. of completed simulations
  for(i in 1:length(res.list)){

    #take one of the results
    res <- res.list[[i]]

    #read out the numbers of patients and dlts it needed
    trials.n.all[,i] <- res$'n.pat'
    trials.dlt.all[,i] <- res$'n.dlt'

    #how many patients were treated at underdoses? how many dlts were caused by these treatments?
    trials.n.under[,i] <- res$'trials.n.under'
    trials.dlt.under[,i] <- res$'trials.dlt.under'

    #how many patients were treated at target doses? how many dlts were caused by these treatments?
    trials.n.target[,i] <- res$'trials.n.target'
    trials.dlt.target[,i] <- res$'trials.dlt.target'

    #how many patients were treated at overdoses? how many dlts were caused by these treatments?
    trials.n.over[,i] <- res$'trials.n.over'
    trials.dlt.over[,i] <- res$'trials.dlt.over'

    #ratio of patients treated at overdoses
    trials.perc.over[,i] <- res$'trials.n.over'/res$'n.pat'
    trials.perc.dlt[,i] <- res$'n.dlt'/res$'n.pat'

    #read out patient numbers per dose level
    if(active.mono1.a){
      trials.p.dose.m1[,i] <- res$'pat.d.mono.1'
    }
    if(active.mono2.a){
      trials.p.dose.m2[,i] <- res$'pat.d.mono.2'
    }
    if(active.mono1.b){
      trials.p.dose.m4[,i] <- res$'pat.d.mono.4'
    }
    if(active.mono2.b){
      trials.p.dose.m5[,i] <- res$'pat.d.mono.5'
    }
    if(active.combi.a){
      trials.p.dose.c1[,i] <- res$'pat.d.combi.1'
    }
    if(active.combi.b){
      trials.p.dose.c2[,i] <- res$'pat.d.combi.2'
    }

    #determine the reason for the end of these trials
    #(stopped early or found MTD?)
    #-------------------------------------------------------------------------
    #------
    #Mono 1
    #------
    if(active.mono1.a==TRUE){

      if(res$'mono.1.stopped' == TRUE){

        results.mono.1["Number of trials", "all doses too toxic"] <- results.mono.1["Number of trials", "all doses too toxic"] +1

      } else if(res$'mono.1.is.MTD' == TRUE){

        #get toxicity and dose name
        name.of.dose <- toString(res$'current.dose.mono.1')
        curr.tox.mono.1 <- res$'curr.tox.mono.1'
        if(curr.tox.mono.1 < dosing.intervals[1]){

          #save which dose was declared MTD
          input.mono.1["MTD declared (n)", name.of.dose ] <- as.numeric(input.mono.1["MTD declared (n)", name.of.dose])+1
          #increase the underdose-counter
          results.mono.1["Number of trials", "underdose"] <- results.mono.1["Number of trials", "underdose"] + 1  #MTD is underdose

        }else if (dosing.intervals[1] <= curr.tox.mono.1 & curr.tox.mono.1 < dosing.intervals[2]){

          #save which dose was declared MTD
          input.mono.1["MTD declared (n)", name.of.dose] <- as.numeric(input.mono.1["MTD declared (n)", name.of.dose])+1
          #increase the targetdose-counter
          results.mono.1["Number of trials", "target dose"] <- results.mono.1["Number of trials", "target dose"] +1    #MTD is in target toxicity interval

        }else if(dosing.intervals[2]<= curr.tox.mono.1){

          #save which dose was declared MTD
          input.mono.1["MTD declared (n)", name.of.dose] <- as.numeric(input.mono.1["MTD declared (n)", name.of.dose])+1
          #increase the overdose counter
          results.mono.1["Number of trials", "overdose"] <- results.mono.1["Number of trials", "overdose"] + 1    #MTD is over dose

        }
      } else if (trials.n.all[1,i] >= max.n.mono1.a){
        #increase the "max number of patients reached before MTD"-counter
        results.mono.1["Number of trials", "max n reached before MTD"] <- results.mono.1["Number of trials", "max n reached before MTD"] +1 #no MTD was found, the trial stopped after the max amount of patients had been treated

      }

    }

    #----------------------------------------------------------------------------------------------------------
    #Mono 1
    if(active.mono1.b==TRUE){
      if(res$'mono.4.stopped' == TRUE){

        results.mono.4["Number of trials", "all doses too toxic"] <- results.mono.4["Number of trials", "all doses too toxic"] +1

      } else if(res$'mono.4.is.MTD' == TRUE){

        name.of.dose <- toString(res$'current.dose.mono.4')
        curr.tox.mono.4 <- res$'curr.tox.mono.4'
        # the MTD was found
        if(curr.tox.mono.4 < dosing.intervals[1]){

          #save which dose was declared MTD
          input.mono.4["MTD declared (n)", name.of.dose ] <- as.numeric(input.mono.4["MTD declared (n)", name.of.dose])+1
          results.mono.4["Number of trials", "underdose"] <- results.mono.4["Number of trials", "underdose"] + 1  #MTD is underdose

        }else if (dosing.intervals[1] <= curr.tox.mono.4 & curr.tox.mono.4 < dosing.intervals[2]){

          #save, which dose was declared MTD
          input.mono.4["MTD declared (n)", name.of.dose] <- as.numeric(input.mono.4["MTD declared (n)", name.of.dose])+1
          results.mono.4["Number of trials", "target dose"] <- results.mono.4["Number of trials", "target dose"] +1    #MTD is in target toxicity interval

        }else if(dosing.intervals[2]<= curr.tox.mono.4){

          #save, which dose was declared MTD
          input.mono.4["MTD declared (n)", name.of.dose] <- as.numeric(input.mono.4["MTD declared (n)", name.of.dose])+1

          results.mono.4["Number of trials", "overdose"] <- results.mono.4["Number of trials", "overdose"] + 1    #MTD is over dose

        }
      } else if (trials.n.all[4,i] >= max.n.mono1.b){

        results.mono.4["Number of trials", "max n reached before MTD"] <- results.mono.4["Number of trials", "max n reached before MTD"] +1 #no MTD was found, the trial stopped after the max amount of patients had been treated

      }

    }


    #----------------------------------------------------------------------------------------------------------
    #------
    #Mono 2
    #------
    if(active.mono2.a==TRUE){
      if(res$'mono.2.stopped' == TRUE){

        results.mono.2["Number of trials", "all doses too toxic"] <- results.mono.2["Number of trials", "all doses too toxic"] +1

      } else if(res$'mono.2.is.MTD' == TRUE){
        name.of.dose <- toString(res$'current.dose.mono.2')
        curr.tox.mono.2 <- res$'curr.tox.mono.2'
        # the MTD was found
        if(curr.tox.mono.2 < dosing.intervals[1]){

          #save which dose was declared MTD
          input.mono.2["MTD declared (n)", name.of.dose] <- as.numeric(input.mono.2["MTD declared (n)", name.of.dose])+1

          results.mono.2["Number of trials", "underdose"] <- results.mono.2["Number of trials", "underdose"] + 1  #MTD is underdose

        }else if (dosing.intervals[1] <= curr.tox.mono.2 & curr.tox.mono.2 < dosing.intervals[2]){

          #save which dose was declared MTD
          input.mono.2["MTD declared (n)", name.of.dose] <- as.numeric(input.mono.2["MTD declared (n)", name.of.dose])+1

          results.mono.2["Number of trials", "target dose"] <- results.mono.2["Number of trials", "target dose"] +1    #MTD is in target toxicity interval

        }else if(dosing.intervals[2]<= curr.tox.mono.2){

          #save which dose was declared MTD
          input.mono.2["MTD declared (n)", name.of.dose] <- as.numeric(input.mono.2["MTD declared (n)", name.of.dose])+1


          results.mono.2["Number of trials", "overdose"] <- results.mono.2["Number of trials", "overdose"] + 1    #MTD is over dose

        }
      } else if (trials.n.all[2,i] >= max.n.mono2.a){

        results.mono.2["Number of trials", "max n reached before MTD"] <- results.mono.2["Number of trials", "max n reached before MTD"] +1 #no MTD was found, the trial stopped after the max amount of patients had been treated

      }

    }

    #----------------------------------------------------------------------------------------------------------
    #Mono 2
    if(active.mono2.b==TRUE){
      if(res$'mono.5.stopped' == TRUE){

        results.mono.5["Number of trials", "all doses too toxic"] <- results.mono.5["Number of trials", "all doses too toxic"] +1

      } else if(res$'mono.5.is.MTD' == TRUE){
        name.of.dose <- toString(res$'current.dose.mono.5')
        curr.tox.mono.5 <- res$'curr.tox.mono.5'
        # the MTD was found
        if(curr.tox.mono.5 < dosing.intervals[1]){

          #save which dose was declared MTD
          input.mono.5["MTD declared (n)", name.of.dose] <- as.numeric(input.mono.5["MTD declared (n)", name.of.dose])+1

          results.mono.5["Number of trials", "underdose"] <- results.mono.5["Number of trials", "underdose"] + 1  #MTD is underdose

        }else if (dosing.intervals[1] <= curr.tox.mono.5 & curr.tox.mono.5 < dosing.intervals[2]){

          #save which dose was declared MTD
          input.mono.5["MTD declared (n)", name.of.dose] <- as.numeric(input.mono.5["MTD declared (n)", name.of.dose])+1

          results.mono.5["Number of trials", "target dose"] <- results.mono.5["Number of trials", "target dose"] +1    #MTD is in target toxicity interval

        }else if(dosing.intervals[2]<= curr.tox.mono.5){

          #save which dose was declared MTD
          input.mono.5["MTD declared (n)", name.of.dose] <- as.numeric(input.mono.5["MTD declared (n)", name.of.dose])+1


          results.mono.5["Number of trials", "overdose"] <- results.mono.5["Number of trials", "overdose"] + 1    #MTD is over dose

        }
      } else if (trials.n.all[5,i] >= max.n.mono2.b){

        results.mono.5["Number of trials", "max n reached before MTD"] <- results.mono.5["Number of trials", "max n reached before MTD"] +1 #no MTD was found, the trial stopped after the max amount of patients had been treated

      }

    }



    #----------------------------------------------------------------------------------------------------------
    #------
    #Combi
    #------
    if(active.combi.a == TRUE){
      if(res$'combi.stopped.1' == TRUE){

        results.combi1["Number of trials", "all doses too toxic"] <- results.combi1["Number of trials", "all doses too toxic"] +1

      } else if(res$'combi.is.MTD.1' == TRUE){

        name.of.dose <- paste(res$'current.dose.1.combi.1', res$'current.dose.2.combi.1', sep = "+")
        curr.tox.combi <- res$'curr.tox.combi.1'
        # the MTD was found
        if(curr.tox.combi < dosing.intervals[1]){

          #save which dose was declared MTD
          input.combi1["MTD declared (n)", name.of.dose] <- as.numeric(input.combi1["MTD declared (n)", name.of.dose])+1

          results.combi1["Number of trials", "underdose"] <- results.combi1["Number of trials", "underdose"] + 1  #MTD is underdose

        }else if (dosing.intervals[1] <= curr.tox.combi & curr.tox.combi < dosing.intervals[2]){

          #save which dose was declared MTD
          input.combi1["MTD declared (n)", name.of.dose] <- as.numeric(input.combi1["MTD declared (n)", name.of.dose]) +1
          results.combi1["Number of trials", "target dose"] <- results.combi1["Number of trials", "target dose"] +1    #MTD is in target toxicity interval

        }else if(dosing.intervals[2]<= curr.tox.combi){

          #save which dose was declared MTD
          input.combi1["MTD declared (n)", name.of.dose] <- as.numeric(input.combi1["MTD declared (n)", name.of.dose])+1
          results.combi1["Number of trials", "overdose"] <- results.combi1["Number of trials", "overdose"] + 1    #MTD is over dose

        }
      } else if (trials.n.all[3,i] >= max.n.combi.a){

        results.combi1["Number of trials", "max n reached before MTD"] <- results.combi1["Number of trials", "max n reached before MTD"] +1 #no MTD was found, the trial stopped after the max amount of patients had been treated

      }

    }

    #----------------------------------------------------------------------------------------------------------
    #Combi
    if(active.combi.b == TRUE){
      if(res$'combi.stopped.2' == TRUE){

        results.combi2["Number of trials", "all doses too toxic"] <- results.combi2["Number of trials", "all doses too toxic"] +1

      } else if(res$'combi.is.MTD.2' == TRUE){

        name.of.dose <- paste(res$'current.dose.1.combi.2', res$'current.dose.2.combi.2', sep = "+")
        curr.tox.combi.2 <- res$'curr.tox.combi.2'
        # the MTD was found
        if(curr.tox.combi.2 < dosing.intervals[1]){

          #save which dose was declared MTD
          input.combi2["MTD declared (n)", name.of.dose] <- as.numeric(input.combi2["MTD declared (n)", name.of.dose])+1

          results.combi2["Number of trials", "underdose"] <- results.combi2["Number of trials", "underdose"] + 1  #MTD is underdose

        }else if (dosing.intervals[1] <= curr.tox.combi.2 & curr.tox.combi.2 < dosing.intervals[2]){

          #save which dose was declared MTD
          input.combi2["MTD declared (n)", name.of.dose] <- as.numeric(input.combi2["MTD declared (n)", name.of.dose]) +1
          results.combi2["Number of trials", "target dose"] <- results.combi2["Number of trials", "target dose"] +1    #MTD is in target toxicity interval

        }else if(dosing.intervals[2]<= curr.tox.combi.2){

          #save which dose was declared MTD
          input.combi2["MTD declared (n)", name.of.dose] <- as.numeric(input.combi2["MTD declared (n)", name.of.dose])+1
          results.combi2["Number of trials", "overdose"] <- results.combi2["Number of trials", "overdose"] + 1    #MTD is over dose

        }
      } else if (trials.n.all[6,i] >= max.n.combi.b){

        results.combi2["Number of trials", "max n reached before MTD"] <- results.combi2["Number of trials", "max n reached before MTD"] +1 #no MTD was found, the trial stopped after the max amount of patients had been treated

      }

    }

  }
  #end of for-loop

  #---------------------------------------------------------------------------
  #SAVE the data for the output
  #---------------------------------------------------------------------------

  #list for saving each of the simulated trials
  simulation.results <- list()

  #the following is pretty straightforward bookkeeping.
  n.studies <- length(res.list)
  #-----
  #COMBO
  #-----
  #prepare the output data for the combo trial
  if(active.combi.a){
    #fill the results matrix (straightforward)
    results.combi1["Percentage", "underdose"] <- 100*results.combi1["Number of trials", "underdose"]/n.studies
    results.combi1["Percentage", "target dose"] <- 100*results.combi1["Number of trials", "target dose"]/n.studies
    results.combi1["Percentage", "overdose"] <- 100*results.combi1["Number of trials", "overdose"]/n.studies
    results.combi1["Percentage", "max n reached before MTD"] <- 100*results.combi1["Number of trials", "max n reached before MTD"]/n.studies
    results.combi1["Percentage", "all doses too toxic"] <- 100*results.combi1["Number of trials", "all doses too toxic"]/n.studies

    #how many trials were not stopped for being to toxic? (used to compute percentages below)
    n.studies.not.tox <- n.studies - results.combi1["Number of trials", "all doses too toxic"]

    if(n.studies.not.tox == 0){
      results.combi1["Percentage not all too toxic", "underdose"] <- ""
      results.combi1["Percentage not all too toxic", "target dose"] <- ""
      results.combi1["Percentage not all too toxic", "overdose"] <- ""
      results.combi1["Percentage not all too toxic", "max n reached before MTD"] <- ""
      results.combi1["Percentage not all too toxic", "all doses too toxic"] <- ""
    } else {
      results.combi1["Percentage not all too toxic", "underdose"] <- 100*results.combi1["Number of trials", "underdose"] / n.studies.not.tox
      results.combi1["Percentage not all too toxic", "target dose"] <- 100*results.combi1["Number of trials", "target dose"] / n.studies.not.tox
      results.combi1["Percentage not all too toxic", "overdose"] <- 100*results.combi1["Number of trials", "overdose"] / n.studies.not.tox
      results.combi1["Percentage not all too toxic", "max n reached before MTD"] <- 100*results.combi1["Number of trials", "max n reached before MTD"] / n.studies.not.tox
      results.combi1["Percentage not all too toxic", "all doses too toxic"] <- NA
    }

    #fill the summary matrix
    summary.combi1["#Pat underdose", "Median"] <- median(trials.n.under[3,], na.rm = TRUE)
    summary.combi1["#Pat target dose", "Median"] <- median(trials.n.target[3,], na.rm = TRUE)
    summary.combi1["#Pat overdose", "Median"] <- median(trials.n.over[3,], na.rm = TRUE)
    summary.combi1["#Pat (all)", "Median"] <- median(trials.n.all[3,], na.rm = TRUE)
    summary.combi1["#DLT underdose", "Median"] <- median(trials.dlt.under[3,], na.rm = TRUE)
    summary.combi1["#DLT target dose", "Median"] <- median(trials.dlt.target[3,], na.rm = TRUE)
    summary.combi1["#DLT overdose", "Median"] <- median(trials.dlt.over[3,], na.rm = TRUE)
    summary.combi1["#DLT (all)", "Median"] <- median(trials.dlt.all[3,], na.rm = TRUE)
    summary.combi1["% overdose", "Median"] <- median(trials.perc.over[3,], na.rm = TRUE)
    summary.combi1["% DLT", "Median"] <- median(trials.perc.dlt[3,], na.rm = TRUE)


    summary.combi1["#Pat underdose", "Mean"] <- mean(trials.n.under[3,], na.rm = TRUE)
    summary.combi1["#Pat target dose", "Mean"] <- mean(trials.n.target[3,], na.rm = TRUE)
    summary.combi1["#Pat overdose", "Mean"] <- mean(trials.n.over[3,], na.rm = TRUE)
    summary.combi1["#Pat (all)", "Mean"] <- mean(trials.n.all[3,], na.rm = TRUE)
    summary.combi1["#DLT underdose", "Mean"] <- mean(trials.dlt.under[3,], na.rm = TRUE)
    summary.combi1["#DLT target dose", "Mean"] <- mean(trials.dlt.target[3,], na.rm = TRUE)
    summary.combi1["#DLT overdose", "Mean"] <- mean(trials.dlt.over[3,], na.rm = TRUE)
    summary.combi1["#DLT (all)", "Mean"] <- mean(trials.dlt.all[3,], na.rm = TRUE)
    summary.combi1["% overdose", "Mean"] <- mean(trials.perc.over[3,], na.rm = TRUE)
    summary.combi1["% DLT", "Mean"] <- mean(trials.perc.dlt[3,], na.rm = TRUE)

    summary.combi1["#Pat underdose", "Min."] <- min(trials.n.under[3,], na.rm = TRUE)
    summary.combi1["#Pat target dose", "Min."] <- min(trials.n.target[3,], na.rm = TRUE)
    summary.combi1["#Pat overdose", "Min."] <- min(trials.n.over[3,], na.rm = TRUE)
    summary.combi1["#Pat (all)", "Min."] <- min(trials.n.all[3,], na.rm = TRUE)
    summary.combi1["#DLT underdose", "Min."] <- min(trials.dlt.under[3,], na.rm = TRUE)
    summary.combi1["#DLT target dose", "Min."] <- min(trials.dlt.target[3,], na.rm = TRUE)
    summary.combi1["#DLT overdose", "Min."] <- min(trials.dlt.over[3,], na.rm = TRUE)
    summary.combi1["#DLT (all)", "Min."] <- min(trials.dlt.all[3,], na.rm = TRUE)
    summary.combi1["% overdose", "Min."] <- min(trials.perc.over[3,], na.rm = TRUE)
    summary.combi1["% DLT", "Min."] <- min(trials.perc.dlt[3,], na.rm = TRUE)

    summary.combi1["#Pat underdose", "Max."] <- max(trials.n.under[3,], na.rm = TRUE)
    summary.combi1["#Pat target dose", "Max."] <- max(trials.n.target[3,], na.rm = TRUE)
    summary.combi1["#Pat overdose", "Max."] <- max(trials.n.over[3,], na.rm = TRUE)
    summary.combi1["#Pat (all)", "Max."] <- max(trials.n.all[3,], na.rm = TRUE)
    summary.combi1["#DLT underdose", "Max."] <- max(trials.dlt.under[3,], na.rm = TRUE)
    summary.combi1["#DLT target dose", "Max."] <- max(trials.dlt.target[3,], na.rm = TRUE)
    summary.combi1["#DLT overdose", "Max."] <- max(trials.dlt.over[3,], na.rm = TRUE)
    summary.combi1["#DLT (all)", "Max."] <- max(trials.dlt.all[3,], na.rm = TRUE)
    summary.combi1["% overdose", "Max."] <- max(trials.perc.over[3,], na.rm = TRUE)
    summary.combi1["% DLT", "Max."] <- max(trials.perc.dlt[3,], na.rm = TRUE)


    summary.combi1["#Pat underdose", "2.5%"] <- quantile(trials.n.under[3,], na.rm = TRUE, probs = 0.025)
    summary.combi1["#Pat target dose", "2.5%"] <- quantile(trials.n.target[3,], na.rm = TRUE, probs = 0.025)
    summary.combi1["#Pat overdose", "2.5%"] <- quantile(trials.n.over[3,], na.rm = TRUE, probs = 0.025)
    summary.combi1["#Pat (all)", "2.5%"] <- quantile(trials.n.all[3,], na.rm = TRUE, probs = 0.025)
    summary.combi1["#DLT underdose", "2.5%"] <- quantile(trials.dlt.under[3,], na.rm = TRUE, probs = 0.025)
    summary.combi1["#DLT target dose", "2.5%"] <- quantile(trials.dlt.target[3,], na.rm = TRUE, probs = 0.025)
    summary.combi1["#DLT overdose", "2.5%"] <- quantile(trials.dlt.over[3,], na.rm = TRUE, probs = 0.025)
    summary.combi1["#DLT (all)", "2.5%"] <- quantile(trials.dlt.all[3,], na.rm = TRUE, probs = 0.025)
    summary.combi1["% overdose", "2.5%"] <- quantile(trials.perc.over[3,], na.rm = TRUE, probs = 0.025)
    summary.combi1["% DLT", "2.5%"] <- quantile(trials.perc.dlt[3,], na.rm = TRUE, probs = 0.025)


    summary.combi1["#Pat underdose", "97.5%"] <- quantile(trials.n.under[3,], na.rm = TRUE, probs = 0.975)
    summary.combi1["#Pat target dose", "97.5%"] <- quantile(trials.n.target[3,], na.rm = TRUE, probs = 0.975)
    summary.combi1["#Pat overdose", "97.5%"] <- quantile(trials.n.over[3,], na.rm = TRUE, probs = 0.975)
    summary.combi1["#Pat (all)", "97.5%"] <- quantile(trials.n.all[3,], na.rm = TRUE, probs = 0.975)
    summary.combi1["#DLT underdose", "97.5%"] <- quantile(trials.dlt.under[3,], na.rm = TRUE, probs = 0.975)
    summary.combi1["#DLT target dose", "97.5%"] <- quantile(trials.dlt.target[3,], na.rm = TRUE, probs = 0.975)
    summary.combi1["#DLT overdose", "97.5%"] <- quantile(trials.dlt.over[3,], na.rm = TRUE, probs = 0.975)
    summary.combi1["#DLT (all)", "97.5%"] <- quantile(trials.dlt.all[3,], na.rm = TRUE, probs = 0.975)
    summary.combi1["% overdose", "97.5%"] <- quantile(trials.perc.over[3,], na.rm = TRUE, probs = 0.975)
    summary.combi1["% DLT", "97.5%"] <- quantile(trials.perc.dlt[3,], na.rm = TRUE, probs = 0.975)

    for(d in 1:length(tox.combi.a)){
      summ.ppd.c1[2, d] <- mean(trials.p.dose.c1[d, ], na.rm=T)
      summ.ppd.c1[3, d] <- median(trials.p.dose.c1[d, ], na.rm=T)
      summ.ppd.c1[4, d] <- min(trials.p.dose.c1[d, ], na.rm=T)
      summ.ppd.c1[5, d] <- max(trials.p.dose.c1[d, ], na.rm=T)
    }



    #append the matrices to the list of output data.
    simulation.results[['results combi.a']] <- results.combi1
    simulation.results[['summary combi.a']] <- summary.combi1
    simulation.results[['MTDs combi.a']] <- input.combi1
    simulation.results[['#pat. combi.a']] <- summ.ppd.c1
  }

  #------------------------------------------------------------------------------------------------------------
  #same as for the other trials
  if(active.combi.b){
    #fill the results matrix (straightforward)

    results.combi2["Percentage", "underdose"] <- 100*results.combi2["Number of trials", "underdose"]/n.studies
    results.combi2["Percentage", "target dose"] <- 100*results.combi2["Number of trials", "target dose"]/n.studies
    results.combi2["Percentage", "overdose"] <- 100*results.combi2["Number of trials", "overdose"]/n.studies
    results.combi2["Percentage", "max n reached before MTD"] <- 100*results.combi2["Number of trials", "max n reached before MTD"]/n.studies
    results.combi2["Percentage", "all doses too toxic"] <- 100*results.combi2["Number of trials", "all doses too toxic"]/n.studies

    n.studies.not.tox <- n.studies - results.combi2["Number of trials", "all doses too toxic"]

    if(n.studies.not.tox == 0){
      results.combi2["Percentage not all too toxic", "underdose"] <- ""
      results.combi2["Percentage not all too toxic", "target dose"] <- ""
      results.combi2["Percentage not all too toxic", "overdose"] <- ""
      results.combi2["Percentage not all too toxic", "max n reached before MTD"] <- ""
      results.combi2["Percentage not all too toxic", "all doses too toxic"] <- ""
    } else {
      results.combi2["Percentage not all too toxic", "underdose"] <- 100*results.combi2["Number of trials", "underdose"] / n.studies.not.tox
      results.combi2["Percentage not all too toxic", "target dose"] <- 100*results.combi2["Number of trials", "target dose"] / n.studies.not.tox
      results.combi2["Percentage not all too toxic", "overdose"] <- 100*results.combi2["Number of trials", "overdose"] / n.studies.not.tox
      results.combi2["Percentage not all too toxic", "max n reached before MTD"] <- 100*results.combi2["Number of trials", "max n reached before MTD"] / n.studies.not.tox
      results.combi2["Percentage not all too toxic", "all doses too toxic"] <- NA
    }

    #fill the summary matrix

    summary.combi2["#Pat underdose", "Median"] <- median(trials.n.under[6,], na.rm = TRUE)
    summary.combi2["#Pat target dose", "Median"] <- median(trials.n.target[6,], na.rm = TRUE)
    summary.combi2["#Pat overdose", "Median"] <- median(trials.n.over[6,], na.rm = TRUE)
    summary.combi2["#Pat (all)", "Median"] <- median(trials.n.all[6,], na.rm = TRUE)
    summary.combi2["#DLT underdose", "Median"] <- median(trials.dlt.under[6,], na.rm = TRUE)
    summary.combi2["#DLT target dose", "Median"] <- median(trials.dlt.target[6,], na.rm = TRUE)
    summary.combi2["#DLT overdose", "Median"] <- median(trials.dlt.over[6,], na.rm = TRUE)
    summary.combi2["#DLT (all)", "Median"] <- median(trials.dlt.all[6,], na.rm = TRUE)
    summary.combi2["% overdose", "Median"] <- median(trials.perc.over[6,], na.rm = TRUE)
    summary.combi2["% DLT", "Median"] <- median(trials.perc.dlt[6,], na.rm = TRUE)

    summary.combi2["#Pat underdose", "Mean"] <- mean(trials.n.under[6,], na.rm = TRUE)
    summary.combi2["#Pat target dose", "Mean"] <- mean(trials.n.target[6,], na.rm = TRUE)
    summary.combi2["#Pat overdose", "Mean"] <- mean(trials.n.over[6,], na.rm = TRUE)
    summary.combi2["#Pat (all)", "Mean"] <- mean(trials.n.all[6,], na.rm = TRUE)
    summary.combi2["#DLT underdose", "Mean"] <- mean(trials.dlt.under[6,], na.rm = TRUE)
    summary.combi2["#DLT target dose", "Mean"] <- mean(trials.dlt.target[6,], na.rm = TRUE)
    summary.combi2["#DLT overdose", "Mean"] <- mean(trials.dlt.over[6,], na.rm = TRUE)
    summary.combi2["#DLT (all)", "Mean"] <- mean(trials.dlt.all[6,], na.rm = TRUE)
    summary.combi2["% overdose", "Mean"] <- mean(trials.perc.over[6,], na.rm = TRUE)
    summary.combi2["% DLT", "Mean"] <- mean(trials.perc.dlt[6,], na.rm = TRUE)

    summary.combi2["#Pat underdose", "Min."] <- min(trials.n.under[6,], na.rm = TRUE)
    summary.combi2["#Pat target dose", "Min."] <- min(trials.n.target[6,], na.rm = TRUE)
    summary.combi2["#Pat overdose", "Min."] <- min(trials.n.over[6,], na.rm = TRUE)
    summary.combi2["#Pat (all)", "Min."] <- min(trials.n.all[6,], na.rm = TRUE)
    summary.combi2["#DLT underdose", "Min."] <- min(trials.dlt.under[6,], na.rm = TRUE)
    summary.combi2["#DLT target dose", "Min."] <- min(trials.dlt.target[6,], na.rm = TRUE)
    summary.combi2["#DLT overdose", "Min."] <- min(trials.dlt.over[6,], na.rm = TRUE)
    summary.combi2["#DLT (all)", "Min."] <- min(trials.dlt.all[6,], na.rm = TRUE)
    summary.combi2["% overdose", "Min."] <- min(trials.perc.over[6,], na.rm = TRUE)
    summary.combi2["% DLT", "Min."] <- min(trials.perc.dlt[6,], na.rm = TRUE)

    summary.combi2["#Pat underdose", "Max."] <- max(trials.n.under[6,], na.rm = TRUE)
    summary.combi2["#Pat target dose", "Max."] <- max(trials.n.target[6,], na.rm = TRUE)
    summary.combi2["#Pat overdose", "Max."] <- max(trials.n.over[6,], na.rm = TRUE)
    summary.combi2["#Pat (all)", "Max."] <- max(trials.n.all[6,], na.rm = TRUE)
    summary.combi2["#DLT underdose", "Max."] <- max(trials.dlt.under[6,], na.rm = TRUE)
    summary.combi2["#DLT target dose", "Max."] <- max(trials.dlt.target[6,], na.rm = TRUE)
    summary.combi2["#DLT overdose", "Max."] <- max(trials.dlt.over[6,], na.rm = TRUE)
    summary.combi2["#DLT (all)", "Max."] <- max(trials.dlt.all[6,], na.rm = TRUE)
    summary.combi2["% overdose", "Max."] <- max(trials.perc.over[6,], na.rm = TRUE)
    summary.combi2["% DLT", "Max."] <- max(trials.perc.dlt[6,], na.rm = TRUE)


    summary.combi2["#Pat underdose", "2.5%"] <- quantile(trials.n.under[6,], na.rm = TRUE, probs = 0.025)
    summary.combi2["#Pat target dose", "2.5%"] <- quantile(trials.n.target[6,], na.rm = TRUE, probs = 0.025)
    summary.combi2["#Pat overdose", "2.5%"] <- quantile(trials.n.over[6,], na.rm = TRUE, probs = 0.025)
    summary.combi2["#Pat (all)", "2.5%"] <- quantile(trials.n.all[6,], na.rm = TRUE, probs = 0.025)
    summary.combi2["#DLT underdose", "2.5%"] <- quantile(trials.dlt.under[6,], na.rm = TRUE, probs = 0.025)
    summary.combi2["#DLT target dose", "2.5%"] <- quantile(trials.dlt.target[6,], na.rm = TRUE, probs = 0.025)
    summary.combi2["#DLT overdose", "2.5%"] <- quantile(trials.dlt.over[6,], na.rm = TRUE, probs = 0.025)
    summary.combi2["#DLT (all)", "2.5%"] <- quantile(trials.dlt.all[6,], na.rm = TRUE, probs = 0.025)
    summary.combi2["% overdose", "2.5%"] <- quantile(trials.perc.over[6,], na.rm = TRUE, probs = 0.025)
    summary.combi2["% DLT", "2.5%"] <- quantile(trials.perc.dlt[6,], na.rm = TRUE, probs = 0.025)


    summary.combi2["#Pat underdose", "97.5%"] <- quantile(trials.n.under[6,], na.rm = TRUE, probs = 0.975)
    summary.combi2["#Pat target dose", "97.5%"] <- quantile(trials.n.target[6,], na.rm = TRUE, probs = 0.975)
    summary.combi2["#Pat overdose", "97.5%"] <- quantile(trials.n.over[6,], na.rm = TRUE, probs = 0.975)
    summary.combi2["#Pat (all)", "97.5%"] <- quantile(trials.n.all[6,], na.rm = TRUE, probs = 0.975)
    summary.combi2["#DLT underdose", "97.5%"] <- quantile(trials.dlt.under[6,], na.rm = TRUE, probs = 0.975)
    summary.combi2["#DLT target dose", "97.5%"] <- quantile(trials.dlt.target[6,], na.rm = TRUE, probs = 0.975)
    summary.combi2["#DLT overdose", "97.5%"] <- quantile(trials.dlt.over[6,], na.rm = TRUE, probs = 0.975)
    summary.combi2["#DLT (all)", "97.5%"] <- quantile(trials.dlt.all[6,], na.rm = TRUE, probs = 0.975)
    summary.combi2["% overdose", "97.5%"] <- quantile(trials.perc.over[6,], na.rm = TRUE, probs = 0.975)
    summary.combi2["% DLT", "97.5%"] <- quantile(trials.perc.dlt[6,], na.rm = TRUE, probs = 0.975)



    for(d in 1:length(tox.combi.b)){
      summ.ppd.c2[2, d] <- mean(trials.p.dose.c2[d, ], na.rm=T)
      summ.ppd.c2[3, d] <- median(trials.p.dose.c2[d, ], na.rm=T)
      summ.ppd.c2[4, d] <- min(trials.p.dose.c2[d, ], na.rm=T)
      summ.ppd.c2[5, d] <- max(trials.p.dose.c2[d, ], na.rm=T)
    }

    simulation.results[['results combi.b']] <- results.combi2
    simulation.results[['summary combi.b']] <- summary.combi2
    simulation.results[['MTDs combi.b']] <- input.combi2
    simulation.results[['#pat. combi.b']] <- summ.ppd.c2
  }


  #------------------------------------------------------------------------------------------------------------
  #same as for the other trials
  #-------------------
  #MONO compound 1
  #-------------------
  #prepare the results of the mono trial for compound 1
  if(active.mono1.a){
    #fill the results matrix (straightforward)
    results.mono.1["Percentage", "underdose"] <- 100*results.mono.1["Number of trials", "underdose"]/n.studies
    results.mono.1["Percentage", "target dose"] <- 100*results.mono.1["Number of trials", "target dose"]/n.studies
    results.mono.1["Percentage", "overdose"] <- 100*results.mono.1["Number of trials", "overdose"]/n.studies
    results.mono.1["Percentage", "max n reached before MTD"] <- 100*results.mono.1["Number of trials", "max n reached before MTD"]/n.studies
    results.mono.1["Percentage", "all doses too toxic"] <- 100*results.mono.1["Number of trials", "all doses too toxic"]/n.studies

    n.studies.not.tox <- n.studies - results.mono.1["Number of trials", "all doses too toxic"]

    if(n.studies.not.tox == 0){
      results.mono.1["Percentage not all too toxic", "underdose"] <- ""
      results.mono.1["Percentage not all too toxic", "target dose"] <- ""
      results.mono.1["Percentage not all too toxic", "overdose"] <- ""
      results.mono.1["Percentage not all too toxic", "max n reached before MTD"] <- ""
      results.mono.1["Percentage not all too toxic", "all doses too toxic"] <- ""
    } else {
      results.mono.1["Percentage not all too toxic", "underdose"] <- 100*results.mono.1["Number of trials", "underdose"] / n.studies.not.tox
      results.mono.1["Percentage not all too toxic", "target dose"] <- 100*results.mono.1["Number of trials", "target dose"] / n.studies.not.tox
      results.mono.1["Percentage not all too toxic", "overdose"] <- 100*results.mono.1["Number of trials", "overdose"] / n.studies.not.tox
      results.mono.1["Percentage not all too toxic", "max n reached before MTD"] <- 100*results.mono.1["Number of trials", "max n reached before MTD"] / n.studies.not.tox
      results.mono.1["Percentage not all too toxic", "all doses too toxic"] <- NA
    }


    #fill the summary matrix

    summary.mono.1["#Pat underdose", "Median"] <- median(trials.n.under[1,], na.rm = TRUE)
    summary.mono.1["#Pat target dose", "Median"] <- median(trials.n.target[1,], na.rm = TRUE)
    summary.mono.1["#Pat overdose", "Median"] <- median(trials.n.over[1,], na.rm = TRUE)
    summary.mono.1["#Pat (all)", "Median"] <- median(trials.n.all[1,], na.rm = TRUE)
    summary.mono.1["#DLT underdose", "Median"] <- median(trials.dlt.under[1,], na.rm = TRUE)
    summary.mono.1["#DLT target dose", "Median"] <- median(trials.dlt.target[1,], na.rm = TRUE)
    summary.mono.1["#DLT overdose", "Median"] <- median(trials.dlt.over[1,], na.rm = TRUE)
    summary.mono.1["#DLT (all)", "Median"] <- median(trials.dlt.all[1,], na.rm = TRUE)
    summary.mono.1["% overdose", "Median"] <- median(trials.perc.over[1,], na.rm = TRUE)
    summary.mono.1["% DLT", "Median"] <- median(trials.perc.dlt[1,], na.rm = TRUE)

    summary.mono.1["#Pat underdose", "Mean"] <- mean(trials.n.under[1,], na.rm = TRUE)
    summary.mono.1["#Pat target dose", "Mean"] <- mean(trials.n.target[1,], na.rm = TRUE)
    summary.mono.1["#Pat overdose", "Mean"] <- mean(trials.n.over[1,], na.rm = TRUE)
    summary.mono.1["#Pat (all)", "Mean"] <- mean(trials.n.all[1,], na.rm = TRUE)
    summary.mono.1["#DLT underdose", "Mean"] <- mean(trials.dlt.under[1,], na.rm = TRUE)
    summary.mono.1["#DLT target dose", "Mean"] <- mean(trials.dlt.target[1,], na.rm = TRUE)
    summary.mono.1["#DLT overdose", "Mean"] <- mean(trials.dlt.over[1,], na.rm = TRUE)
    summary.mono.1["#DLT (all)", "Mean"] <- mean(trials.dlt.all[1,], na.rm = TRUE)
    summary.mono.1["% overdose", "Mean"] <- mean(trials.perc.over[1,], na.rm = TRUE)
    summary.mono.1["% DLT", "Mean"] <- mean(trials.perc.dlt[1,], na.rm = TRUE)

    summary.mono.1["#Pat underdose", "Min."] <- min(trials.n.under[1,], na.rm = TRUE)
    summary.mono.1["#Pat target dose", "Min."] <- min(trials.n.target[1,], na.rm = TRUE)
    summary.mono.1["#Pat overdose", "Min."] <- min(trials.n.over[1,], na.rm = TRUE)
    summary.mono.1["#Pat (all)", "Min."] <- min(trials.n.all[1,], na.rm = TRUE)
    summary.mono.1["#DLT underdose", "Min."] <- min(trials.dlt.under[1,], na.rm = TRUE)
    summary.mono.1["#DLT target dose", "Min."] <- min(trials.dlt.target[1,], na.rm = TRUE)
    summary.mono.1["#DLT overdose", "Min."] <- min(trials.dlt.over[1,], na.rm = TRUE)
    summary.mono.1["#DLT (all)", "Min."] <- min(trials.dlt.all[1,], na.rm = TRUE)
    summary.mono.1["% overdose", "Min."] <- min(trials.perc.over[1,], na.rm = TRUE)
    summary.mono.1["% DLT", "Min."] <- min(trials.perc.dlt[1,], na.rm = TRUE)

    summary.mono.1["#Pat underdose", "Max."] <- max(trials.n.under[1,], na.rm = TRUE)
    summary.mono.1["#Pat target dose", "Max."] <- max(trials.n.target[1,], na.rm = TRUE)
    summary.mono.1["#Pat overdose", "Max."] <- max(trials.n.over[1,], na.rm = TRUE)
    summary.mono.1["#Pat (all)", "Max."] <- max(trials.n.all[1,], na.rm = TRUE)
    summary.mono.1["#DLT underdose", "Max."] <- max(trials.dlt.under[1,], na.rm = TRUE)
    summary.mono.1["#DLT target dose", "Max."] <- max(trials.dlt.target[1,], na.rm = TRUE)
    summary.mono.1["#DLT overdose", "Max."] <- max(trials.dlt.over[1,], na.rm = TRUE)
    summary.mono.1["#DLT (all)", "Max."] <- max(trials.dlt.all[1,], na.rm = TRUE)
    summary.mono.1["% overdose", "Max."] <- max(trials.perc.over[1,], na.rm = TRUE)
    summary.mono.1["% DLT", "Max."] <- max(trials.perc.dlt[1,], na.rm = TRUE)


    summary.mono.1["#Pat underdose", "2.5%"] <- quantile(trials.n.under[1,], na.rm = TRUE, probs = 0.025)
    summary.mono.1["#Pat target dose", "2.5%"] <- quantile(trials.n.target[1,], na.rm = TRUE, probs = 0.025)
    summary.mono.1["#Pat overdose", "2.5%"] <- quantile(trials.n.over[1,], na.rm = TRUE, probs = 0.025)
    summary.mono.1["#Pat (all)", "2.5%"] <- quantile(trials.n.all[1,], na.rm = TRUE, probs = 0.025)
    summary.mono.1["#DLT underdose", "2.5%"] <- quantile(trials.dlt.under[1,], na.rm = TRUE, probs = 0.025)
    summary.mono.1["#DLT target dose", "2.5%"] <- quantile(trials.dlt.target[1,], na.rm = TRUE, probs = 0.025)
    summary.mono.1["#DLT overdose", "2.5%"] <- quantile(trials.dlt.over[1,], na.rm = TRUE, probs = 0.025)
    summary.mono.1["#DLT (all)", "2.5%"] <- quantile(trials.dlt.all[1,], na.rm = TRUE, probs = 0.025)
    summary.mono.1["% overdose", "2.5%"] <- quantile(trials.perc.over[1,], na.rm = TRUE, probs = 0.025)
    summary.mono.1["% DLT", "2.5%"] <- quantile(trials.perc.dlt[1,], na.rm = TRUE, probs = 0.025)


    summary.mono.1["#Pat underdose", "97.5%"] <- quantile(trials.n.under[1,], na.rm = TRUE, probs = 0.975)
    summary.mono.1["#Pat target dose", "97.5%"] <- quantile(trials.n.target[1,], na.rm = TRUE, probs = 0.975)
    summary.mono.1["#Pat overdose", "97.5%"] <- quantile(trials.n.over[1,], na.rm = TRUE, probs = 0.975)
    summary.mono.1["#Pat (all)", "97.5%"] <- quantile(trials.n.all[1,], na.rm = TRUE, probs = 0.975)
    summary.mono.1["#DLT underdose", "97.5%"] <- quantile(trials.dlt.under[1,], na.rm = TRUE, probs = 0.975)
    summary.mono.1["#DLT target dose", "97.5%"] <- quantile(trials.dlt.target[1,], na.rm = TRUE, probs = 0.975)
    summary.mono.1["#DLT overdose", "97.5%"] <- quantile(trials.dlt.over[1,], na.rm = TRUE, probs = 0.975)
    summary.mono.1["#DLT (all)", "97.5%"] <- quantile(trials.dlt.all[1,], na.rm = TRUE, probs = 0.975)
    summary.mono.1["% overdose", "97.5%"] <- quantile(trials.perc.over[1,], na.rm = TRUE, probs = 0.975)
    summary.mono.1["% DLT", "97.5%"] <- quantile(trials.perc.dlt[1,], na.rm = TRUE, probs = 0.975)


    for(d in 1:length(doses.mono1.a)){
      summ.ppd.m1a[2, d] <- mean(trials.p.dose.m1[d, ], na.rm=T)
      summ.ppd.m1a[3, d] <- median(trials.p.dose.m1[d, ], na.rm=T)
      summ.ppd.m1a[4, d] <- min(trials.p.dose.m1[d, ], na.rm=T)
      summ.ppd.m1a[5, d] <- max(trials.p.dose.m1[d, ], na.rm=T)
    }

    simulation.results[['results mono1.a']] <- results.mono.1
    simulation.results[['summary mono1.a']] <- summary.mono.1
    simulation.results[['MTDs mono1.a']] <- input.mono.1
    simulation.results[['#pat. mono1.a']] <- summ.ppd.m1a


  }


  #------------------------------------------------------------------------------------------------------------
  #same as for the other trials
  if(active.mono1.b){
    #fill the results matrix (straightforward)
    results.mono.4["Percentage", "underdose"] <- 100*results.mono.4["Number of trials", "underdose"]/n.studies
    results.mono.4["Percentage", "target dose"] <- 100*results.mono.4["Number of trials", "target dose"]/n.studies
    results.mono.4["Percentage", "overdose"] <- 100*results.mono.4["Number of trials", "overdose"]/n.studies
    results.mono.4["Percentage", "max n reached before MTD"] <- 100*results.mono.4["Number of trials", "max n reached before MTD"]/n.studies
    results.mono.4["Percentage", "all doses too toxic"] <- 100*results.mono.4["Number of trials", "all doses too toxic"]/n.studies

    n.studies.not.tox <- n.studies - results.mono.4["Number of trials", "all doses too toxic"]

    if(n.studies.not.tox == 0){
      results.mono.4["Percentage not all too toxic", "underdose"] <- ""
      results.mono.4["Percentage not all too toxic", "target dose"] <- ""
      results.mono.4["Percentage not all too toxic", "overdose"] <- ""
      results.mono.4["Percentage not all too toxic", "max n reached before MTD"] <- ""
      results.mono.4["Percentage not all too toxic", "all doses too toxic"] <- ""
    } else {
      results.mono.4["Percentage not all too toxic", "underdose"] <- 100*results.mono.4["Number of trials", "underdose"] / n.studies.not.tox
      results.mono.4["Percentage not all too toxic", "target dose"] <- 100*results.mono.4["Number of trials", "target dose"] / n.studies.not.tox
      results.mono.4["Percentage not all too toxic", "overdose"] <- 100*results.mono.4["Number of trials", "overdose"] / n.studies.not.tox
      results.mono.4["Percentage not all too toxic", "max n reached before MTD"] <- 100*results.mono.4["Number of trials", "max n reached before MTD"] / n.studies.not.tox
      results.mono.4["Percentage not all too toxic", "all doses too toxic"] <- NA
    }

    #fill the summary matrix

    summary.mono.4["#Pat underdose", "Median"] <- median(trials.n.under[4,], na.rm = TRUE)
    summary.mono.4["#Pat target dose", "Median"] <- median(trials.n.target[4,], na.rm = TRUE)
    summary.mono.4["#Pat overdose", "Median"] <- median(trials.n.over[4,], na.rm = TRUE)
    summary.mono.4["#Pat (all)", "Median"] <- median(trials.n.all[4,], na.rm = TRUE)
    summary.mono.4["#DLT underdose", "Median"] <- median(trials.dlt.under[4,], na.rm = TRUE)
    summary.mono.4["#DLT target dose", "Median"] <- median(trials.dlt.target[4,], na.rm = TRUE)
    summary.mono.4["#DLT overdose", "Median"] <- median(trials.dlt.over[4,], na.rm = TRUE)
    summary.mono.4["#DLT (all)", "Median"] <- median(trials.dlt.all[4,], na.rm = TRUE)
    summary.mono.4["% overdose", "Median"] <- median(trials.perc.over[4,], na.rm = TRUE)
    summary.mono.4["% DLT", "Median"] <- median(trials.perc.dlt[4,], na.rm = TRUE)

    summary.mono.4["#Pat underdose", "Mean"] <- mean(trials.n.under[4,], na.rm = TRUE)
    summary.mono.4["#Pat target dose", "Mean"] <- mean(trials.n.target[4,], na.rm = TRUE)
    summary.mono.4["#Pat overdose", "Mean"] <- mean(trials.n.over[4,], na.rm = TRUE)
    summary.mono.4["#Pat (all)", "Mean"] <- mean(trials.n.all[4,], na.rm = TRUE)
    summary.mono.4["#DLT underdose", "Mean"] <- mean(trials.dlt.under[4,], na.rm = TRUE)
    summary.mono.4["#DLT target dose", "Mean"] <- mean(trials.dlt.target[4,], na.rm = TRUE)
    summary.mono.4["#DLT overdose", "Mean"] <- mean(trials.dlt.over[4,], na.rm = TRUE)
    summary.mono.4["#DLT (all)", "Mean"] <- mean(trials.dlt.all[4,], na.rm = TRUE)
    summary.mono.4["% overdose", "Mean"] <- mean(trials.perc.over[4,], na.rm = TRUE)
    summary.mono.4["% DLT", "Mean"] <- mean(trials.perc.dlt[4,], na.rm = TRUE)

    summary.mono.4["#Pat underdose", "Min."] <- min(trials.n.under[4,], na.rm = TRUE)
    summary.mono.4["#Pat target dose", "Min."] <- min(trials.n.target[4,], na.rm = TRUE)
    summary.mono.4["#Pat overdose", "Min."] <- min(trials.n.over[4,], na.rm = TRUE)
    summary.mono.4["#Pat (all)", "Min."] <- min(trials.n.all[4,], na.rm = TRUE)
    summary.mono.4["#DLT underdose", "Min."] <- min(trials.dlt.under[4,], na.rm = TRUE)
    summary.mono.4["#DLT target dose", "Min."] <- min(trials.dlt.target[4,], na.rm = TRUE)
    summary.mono.4["#DLT overdose", "Min."] <- min(trials.dlt.over[4,], na.rm = TRUE)
    summary.mono.4["#DLT (all)", "Min."] <- min(trials.dlt.all[4,], na.rm = TRUE)
    summary.mono.4["% overdose", "Min."] <- min(trials.perc.over[4,], na.rm = TRUE)
    summary.mono.4["% DLT", "Min."] <- min(trials.perc.dlt[4,], na.rm = TRUE)

    summary.mono.4["#Pat underdose", "Max."] <- max(trials.n.under[4,], na.rm = TRUE)
    summary.mono.4["#Pat target dose", "Max."] <- max(trials.n.target[4,], na.rm = TRUE)
    summary.mono.4["#Pat overdose", "Max."] <- max(trials.n.over[4,], na.rm = TRUE)
    summary.mono.4["#Pat (all)", "Max."] <- max(trials.n.all[4,], na.rm = TRUE)
    summary.mono.4["#DLT underdose", "Max."] <- max(trials.dlt.under[4,], na.rm = TRUE)
    summary.mono.4["#DLT target dose", "Max."] <- max(trials.dlt.target[4,], na.rm = TRUE)
    summary.mono.4["#DLT overdose", "Max."] <- max(trials.dlt.over[4,], na.rm = TRUE)
    summary.mono.4["#DLT (all)", "Max."] <- max(trials.dlt.all[4,], na.rm = TRUE)
    summary.mono.4["% overdose", "Max."] <- max(trials.perc.over[4,], na.rm = TRUE)
    summary.mono.4["% DLT", "Max."] <- max(trials.perc.dlt[4,], na.rm = TRUE)


    summary.mono.4["#Pat underdose", "2.5%"] <- quantile(trials.n.under[4,], na.rm = TRUE, probs = 0.025)
    summary.mono.4["#Pat target dose", "2.5%"] <- quantile(trials.n.target[4,], na.rm = TRUE, probs = 0.025)
    summary.mono.4["#Pat overdose", "2.5%"] <- quantile(trials.n.over[4,], na.rm = TRUE, probs = 0.025)
    summary.mono.4["#Pat (all)", "2.5%"] <- quantile(trials.n.all[4,], na.rm = TRUE, probs = 0.025)
    summary.mono.4["#DLT underdose", "2.5%"] <- quantile(trials.dlt.under[4,], na.rm = TRUE, probs = 0.025)
    summary.mono.4["#DLT target dose", "2.5%"] <- quantile(trials.dlt.target[4,], na.rm = TRUE, probs = 0.025)
    summary.mono.4["#DLT overdose", "2.5%"] <- quantile(trials.dlt.over[4,], na.rm = TRUE, probs = 0.025)
    summary.mono.4["#DLT (all)", "2.5%"] <- quantile(trials.dlt.all[4,], na.rm = TRUE, probs = 0.025)
    summary.mono.4["% overdose", "2.5%"] <- quantile(trials.perc.over[4,], na.rm = TRUE, probs = 0.025)
    summary.mono.4["% DLT", "2.5%"] <- quantile(trials.perc.dlt[4,], na.rm = TRUE, probs = 0.025)


    summary.mono.4["#Pat underdose", "97.5%"] <- quantile(trials.n.under[4,], na.rm = TRUE, probs = 0.975)
    summary.mono.4["#Pat target dose", "97.5%"] <- quantile(trials.n.target[4,], na.rm = TRUE, probs = 0.975)
    summary.mono.4["#Pat overdose", "97.5%"] <- quantile(trials.n.over[4,], na.rm = TRUE, probs = 0.975)
    summary.mono.4["#Pat (all)", "97.5%"] <- quantile(trials.n.all[4,], na.rm = TRUE, probs = 0.975)
    summary.mono.4["#DLT underdose", "97.5%"] <- quantile(trials.dlt.under[4,], na.rm = TRUE, probs = 0.975)
    summary.mono.4["#DLT target dose", "97.5%"] <- quantile(trials.dlt.target[4,], na.rm = TRUE, probs = 0.975)
    summary.mono.4["#DLT overdose", "97.5%"] <- quantile(trials.dlt.over[4,], na.rm = TRUE, probs = 0.975)
    summary.mono.4["#DLT (all)", "97.5%"] <- quantile(trials.dlt.all[4,], na.rm = TRUE, probs = 0.975)
    summary.mono.4["% overdose", "97.5%"] <- quantile(trials.perc.over[4,], na.rm = TRUE, probs = 0.975)
    summary.mono.4["% DLT", "97.5%"] <- quantile(trials.perc.dlt[4,], na.rm = TRUE, probs = 0.975)

    for(d in 1:length(doses.mono1.b)){
      summ.ppd.m1b[2, d] <- mean(trials.p.dose.m4[d, ], na.rm=T)
      summ.ppd.m1b[3, d] <- median(trials.p.dose.m4[d, ], na.rm=T)
      summ.ppd.m1b[4, d] <- min(trials.p.dose.m4[d, ], na.rm=T)
      summ.ppd.m1b[5, d] <- max(trials.p.dose.m4[d, ], na.rm=T)
    }


    simulation.results[['results mono1.b']] <- results.mono.4
    simulation.results[['summary mono1.b']] <- summary.mono.4
    simulation.results[['MTDs mono1.b']] <- input.mono.4
    simulation.results[['#pat. mono1.b']] <- summ.ppd.m1b

  }


  #------------------------------------------------------------------------------------------------------------
  #same as for the other trials
  #-------------------
  #MONO compound 2
  #-------------------
  #prepare the results of the mono trial for compound 2
  if(active.mono2.a){
    results.mono.2["Percentage", "underdose"] <- 100*results.mono.2["Number of trials", "underdose"]/n.studies
    results.mono.2["Percentage", "target dose"] <- 100*results.mono.2["Number of trials", "target dose"]/n.studies
    results.mono.2["Percentage", "overdose"] <- 100*results.mono.2["Number of trials", "overdose"]/n.studies
    results.mono.2["Percentage", "max n reached before MTD"] <- 100*results.mono.2["Number of trials", "max n reached before MTD"]/n.studies
    results.mono.2["Percentage", "all doses too toxic"] <- 100*results.mono.2["Number of trials", "all doses too toxic"]/n.studies

    n.studies.not.tox <- n.studies - results.mono.2["Number of trials", "all doses too toxic"]

    if(n.studies.not.tox == 0){
      results.mono.2["Percentage not all too toxic", "underdose"] <- ""
      results.mono.2["Percentage not all too toxic", "target dose"] <- ""
      results.mono.2["Percentage not all too toxic", "overdose"] <- ""
      results.mono.2["Percentage not all too toxic", "max n reached before MTD"] <- ""
      results.mono.2["Percentage not all too toxic", "all doses too toxic"] <- ""
    } else {
      results.mono.2["Percentage not all too toxic", "underdose"] <- 100*results.mono.2["Number of trials", "underdose"] / n.studies.not.tox
      results.mono.2["Percentage not all too toxic", "target dose"] <- 100*results.mono.2["Number of trials", "target dose"] / n.studies.not.tox
      results.mono.2["Percentage not all too toxic", "overdose"] <- 100*results.mono.2["Number of trials", "overdose"] / n.studies.not.tox
      results.mono.2["Percentage not all too toxic", "max n reached before MTD"] <- 100*results.mono.2["Number of trials", "max n reached before MTD"] / n.studies.not.tox
      results.mono.2["Percentage not all too toxic", "all doses too toxic"] <- NA
    }

    #fill the summary matrix

    summary.mono.2["#Pat underdose", "Median"] <- median(trials.n.under[2,], na.rm = TRUE)
    summary.mono.2["#Pat target dose", "Median"] <- median(trials.n.target[2,], na.rm = TRUE)
    summary.mono.2["#Pat overdose", "Median"] <- median(trials.n.over[2,], na.rm = TRUE)
    summary.mono.2["#Pat (all)", "Median"] <- median(trials.n.all[2,], na.rm = TRUE)

    summary.mono.2["#DLT underdose", "Median"] <- median(trials.dlt.under[2,], na.rm = TRUE)
    summary.mono.2["#DLT target dose", "Median"] <- median(trials.dlt.target[2,], na.rm = TRUE)
    summary.mono.2["#DLT overdose", "Median"] <- median(trials.dlt.over[2,], na.rm = TRUE)
    summary.mono.2["#DLT (all)", "Median"] <- median(trials.dlt.all[2,], na.rm = TRUE)
    summary.mono.2["% overdose", "Median"] <- median(trials.perc.over[2,], na.rm = TRUE)
    summary.mono.2["% DLT", "Median"] <- median(trials.perc.dlt[2,], na.rm = TRUE)

    summary.mono.2["#Pat underdose", "Mean"] <- mean(trials.n.under[2,], na.rm = TRUE)
    summary.mono.2["#Pat target dose", "Mean"] <- mean(trials.n.target[2,], na.rm = TRUE)
    summary.mono.2["#Pat overdose", "Mean"] <- mean(trials.n.over[2,], na.rm = TRUE)
    summary.mono.2["#Pat (all)", "Mean"] <- mean(trials.n.all[2,], na.rm = TRUE)
    summary.mono.2["#DLT underdose", "Mean"] <- mean(trials.dlt.under[2,], na.rm = TRUE)
    summary.mono.2["#DLT target dose", "Mean"] <- mean(trials.dlt.target[2,], na.rm = TRUE)
    summary.mono.2["#DLT overdose", "Mean"] <- mean(trials.dlt.over[2,], na.rm = TRUE)
    summary.mono.2["#DLT (all)", "Mean"] <- mean(trials.dlt.all[2,], na.rm = TRUE)
    summary.mono.2["% overdose", "Mean"] <- mean(trials.perc.over[2,], na.rm = TRUE)
    summary.mono.2["% DLT", "Mean"] <- mean(trials.perc.dlt[2,], na.rm = TRUE)

    summary.mono.2["#Pat underdose", "Min."] <- min(trials.n.under[2,], na.rm = TRUE)
    summary.mono.2["#Pat target dose", "Min."] <- min(trials.n.target[2,], na.rm = TRUE)
    summary.mono.2["#Pat overdose", "Min."] <- min(trials.n.over[2,], na.rm = TRUE)
    summary.mono.2["#Pat (all)", "Min."] <- min(trials.n.all[2,], na.rm = TRUE)
    summary.mono.2["#DLT underdose", "Min."] <- min(trials.dlt.under[2,], na.rm = TRUE)
    summary.mono.2["#DLT target dose", "Min."] <- min(trials.dlt.target[2,], na.rm = TRUE)
    summary.mono.2["#DLT overdose", "Min."] <- min(trials.dlt.over[2,], na.rm = TRUE)
    summary.mono.2["#DLT (all)", "Min."] <- min(trials.dlt.all[2,], na.rm = TRUE)
    summary.mono.2["% overdose", "Min."] <- min(trials.perc.over[2,], na.rm = TRUE)
    summary.mono.2["% DLT", "Min."] <- min(trials.perc.dlt[2,], na.rm = TRUE)

    summary.mono.2["#Pat underdose", "Max."] <- max(trials.n.under[2,], na.rm = TRUE)
    summary.mono.2["#Pat target dose", "Max."] <- max(trials.n.target[2,], na.rm = TRUE)
    summary.mono.2["#Pat overdose", "Max."] <- max(trials.n.over[2,], na.rm = TRUE)
    summary.mono.2["#Pat (all)", "Max."] <- max(trials.n.all[2,], na.rm = TRUE)
    summary.mono.2["#DLT underdose", "Max."] <- max(trials.dlt.under[2,], na.rm = TRUE)
    summary.mono.2["#DLT target dose", "Max."] <- max(trials.dlt.target[2,], na.rm = TRUE)
    summary.mono.2["#DLT overdose", "Max."] <- max(trials.dlt.over[2,], na.rm = TRUE)
    summary.mono.2["#DLT (all)", "Max."] <- max(trials.dlt.all[2,], na.rm = TRUE)
    summary.mono.2["% overdose", "Max."] <- max(trials.perc.over[2,], na.rm = TRUE)
    summary.mono.2["% DLT", "Max."] <- max(trials.perc.dlt[2,], na.rm = TRUE)



    summary.mono.2["#Pat underdose", "2.5%"] <- quantile(trials.n.under[2,], na.rm = TRUE, probs = 0.025)
    summary.mono.2["#Pat target dose", "2.5%"] <- quantile(trials.n.target[2,], na.rm = TRUE, probs = 0.025)
    summary.mono.2["#Pat overdose", "2.5%"] <- quantile(trials.n.over[2,], na.rm = TRUE, probs = 0.025)
    summary.mono.2["#Pat (all)", "2.5%"] <- quantile(trials.n.all[2,], na.rm = TRUE, probs = 0.025)
    summary.mono.2["#DLT underdose", "2.5%"] <- quantile(trials.dlt.under[2,], na.rm = TRUE, probs = 0.025)
    summary.mono.2["#DLT target dose", "2.5%"] <- quantile(trials.dlt.target[2,], na.rm = TRUE, probs = 0.025)
    summary.mono.2["#DLT overdose", "2.5%"] <- quantile(trials.dlt.over[2,], na.rm = TRUE, probs = 0.025)
    summary.mono.2["#DLT (all)", "2.5%"] <- quantile(trials.dlt.all[2,], na.rm = TRUE, probs = 0.025)
    summary.mono.2["% overdose", "2.5%"] <- quantile(trials.perc.over[2,], na.rm = TRUE, probs = 0.025)
    summary.mono.2["% DLT", "2.5%"] <- quantile(trials.perc.dlt[2,], na.rm = TRUE, probs = 0.025)


    summary.mono.2["#Pat underdose", "97.5%"] <- quantile(trials.n.under[2,], na.rm = TRUE, probs = 0.975)
    summary.mono.2["#Pat target dose", "97.5%"] <- quantile(trials.n.target[2,], na.rm = TRUE, probs = 0.975)
    summary.mono.2["#Pat overdose", "97.5%"] <- quantile(trials.n.over[2,], na.rm = TRUE, probs = 0.975)
    summary.mono.2["#Pat (all)", "97.5%"] <- quantile(trials.n.all[2,], na.rm = TRUE, probs = 0.975)
    summary.mono.2["#DLT underdose", "97.5%"] <- quantile(trials.dlt.under[2,], na.rm = TRUE, probs = 0.975)
    summary.mono.2["#DLT target dose", "97.5%"] <- quantile(trials.dlt.target[2,], na.rm = TRUE, probs = 0.975)
    summary.mono.2["#DLT overdose", "97.5%"] <- quantile(trials.dlt.over[2,], na.rm = TRUE, probs = 0.975)
    summary.mono.2["#DLT (all)", "97.5%"] <- quantile(trials.dlt.all[2,], na.rm = TRUE, probs = 0.975)
    summary.mono.2["% overdose", "97.5%"] <- quantile(trials.perc.over[2,], na.rm = TRUE, probs = 0.975)
    summary.mono.2["% DLT", "97.5%"] <- quantile(trials.perc.dlt[2,], na.rm = TRUE, probs = 0.975)

    for(d in 1:length(doses.mono2.a)){
      summ.ppd.m2a[2, d] <- mean(trials.p.dose.m2[d, ], na.rm=T)
      summ.ppd.m2a[3, d] <- median(trials.p.dose.m2[d, ], na.rm=T)
      summ.ppd.m2a[4, d] <- min(trials.p.dose.m2[d, ], na.rm=T)
      summ.ppd.m2a[5, d] <- max(trials.p.dose.m2[d, ], na.rm=T)
    }

    simulation.results[['results mono2.a']] <- results.mono.2
    simulation.results[['summary mono2.a']] <- summary.mono.2
    simulation.results[['MTDs mono2.a']] <- input.mono.2
    simulation.results[['#pat. mono2.a']] <- summ.ppd.m2a
  }

  if(active.mono2.b){
    results.mono.5["Percentage", "underdose"] <- 100*results.mono.5["Number of trials", "underdose"]/n.studies
    results.mono.5["Percentage", "target dose"] <- 100*results.mono.5["Number of trials", "target dose"]/n.studies
    results.mono.5["Percentage", "overdose"] <- 100*results.mono.5["Number of trials", "overdose"]/n.studies
    results.mono.5["Percentage", "max n reached before MTD"] <- 100*results.mono.5["Number of trials", "max n reached before MTD"]/n.studies
    results.mono.5["Percentage", "all doses too toxic"] <- 100*results.mono.5["Number of trials", "all doses too toxic"]/n.studies

    n.studies.not.tox <- n.studies - results.mono.5["Number of trials", "all doses too toxic"]

    if(n.studies.not.tox == 0){
      results.mono.5["Percentage not all too toxic", "underdose"] <- ""
      results.mono.5["Percentage not all too toxic", "target dose"] <- ""
      results.mono.5["Percentage not all too toxic", "overdose"] <- ""
      results.mono.5["Percentage not all too toxic", "max n reached before MTD"] <- ""
      results.mono.5["Percentage not all too toxic", "all doses too toxic"] <- ""
    } else {
      results.mono.5["Percentage not all too toxic", "underdose"] <- 100*results.mono.5["Number of trials", "underdose"] / n.studies.not.tox
      results.mono.5["Percentage not all too toxic", "target dose"] <- 100*results.mono.5["Number of trials", "target dose"] / n.studies.not.tox
      results.mono.5["Percentage not all too toxic", "overdose"] <- 100*results.mono.5["Number of trials", "overdose"] / n.studies.not.tox
      results.mono.5["Percentage not all too toxic", "max n reached before MTD"] <- 100*results.mono.5["Number of trials", "max n reached before MTD"] / n.studies.not.tox
      results.mono.5["Percentage not all too toxic", "all doses too toxic"] <- NA
    }

    #fill the summary matrix

    summary.mono.5["#Pat underdose", "Median"] <- median(trials.n.under[5,], na.rm = TRUE)
    summary.mono.5["#Pat target dose", "Median"] <- median(trials.n.target[5,], na.rm = TRUE)
    summary.mono.5["#Pat overdose", "Median"] <- median(trials.n.over[5,], na.rm = TRUE)
    summary.mono.5["#Pat (all)", "Median"] <- median(trials.n.all[5,], na.rm = TRUE)

    summary.mono.5["#DLT underdose", "Median"] <- median(trials.dlt.under[5,], na.rm = TRUE)
    summary.mono.5["#DLT target dose", "Median"] <- median(trials.dlt.target[5,], na.rm = TRUE)
    summary.mono.5["#DLT overdose", "Median"] <- median(trials.dlt.over[5,], na.rm = TRUE)
    summary.mono.5["#DLT (all)", "Median"] <- median(trials.dlt.all[5,], na.rm = TRUE)
    summary.mono.5["% overdose", "Median"] <- median(trials.perc.over[5,], na.rm = TRUE)
    summary.mono.5["% DLT", "Median"] <- median(trials.perc.dlt[5,], na.rm = TRUE)

    summary.mono.5["#Pat underdose", "Mean"] <- mean(trials.n.under[5,], na.rm = TRUE)
    summary.mono.5["#Pat target dose", "Mean"] <- mean(trials.n.target[5,], na.rm = TRUE)
    summary.mono.5["#Pat overdose", "Mean"] <- mean(trials.n.over[5,], na.rm = TRUE)
    summary.mono.5["#Pat (all)", "Mean"] <- mean(trials.n.all[5,], na.rm = TRUE)
    summary.mono.5["#DLT underdose", "Mean"] <- mean(trials.dlt.under[5,], na.rm = TRUE)
    summary.mono.5["#DLT target dose", "Mean"] <- mean(trials.dlt.target[5,], na.rm = TRUE)
    summary.mono.5["#DLT overdose", "Mean"] <- mean(trials.dlt.over[5,], na.rm = TRUE)
    summary.mono.5["#DLT (all)", "Mean"] <- mean(trials.dlt.all[5,], na.rm = TRUE)
    summary.mono.5["% overdose", "Mean"] <- mean(trials.perc.over[5,], na.rm = TRUE)
    summary.mono.5["% DLT", "Mean"] <- mean(trials.perc.dlt[5,], na.rm = TRUE)

    summary.mono.5["#Pat underdose", "Min."] <- min(trials.n.under[5,], na.rm = TRUE)
    summary.mono.5["#Pat target dose", "Min."] <- min(trials.n.target[5,], na.rm = TRUE)
    summary.mono.5["#Pat overdose", "Min."] <- min(trials.n.over[5,], na.rm = TRUE)
    summary.mono.5["#Pat (all)", "Min."] <- min(trials.n.all[5,], na.rm = TRUE)
    summary.mono.5["#DLT underdose", "Min."] <- min(trials.dlt.under[5,], na.rm = TRUE)
    summary.mono.5["#DLT target dose", "Min."] <- min(trials.dlt.target[5,], na.rm = TRUE)
    summary.mono.5["#DLT overdose", "Min."] <- min(trials.dlt.over[5,], na.rm = TRUE)
    summary.mono.5["#DLT (all)", "Min."] <- min(trials.dlt.all[5,], na.rm = TRUE)
    summary.mono.5["% overdose", "Min."] <- min(trials.perc.over[5,], na.rm = TRUE)
    summary.mono.5["% DLT", "Min."] <- min(trials.perc.dlt[5,], na.rm = TRUE)

    summary.mono.5["#Pat underdose", "Max."] <- max(trials.n.under[5,], na.rm = TRUE)
    summary.mono.5["#Pat target dose", "Max."] <- max(trials.n.target[5,], na.rm = TRUE)
    summary.mono.5["#Pat overdose", "Max."] <- max(trials.n.over[5,], na.rm = TRUE)
    summary.mono.5["#Pat (all)", "Max."] <- max(trials.n.all[5,], na.rm = TRUE)
    summary.mono.5["#DLT underdose", "Max."] <- max(trials.dlt.under[5,], na.rm = TRUE)
    summary.mono.5["#DLT target dose", "Max."] <- max(trials.dlt.target[5,], na.rm = TRUE)
    summary.mono.5["#DLT overdose", "Max."] <- max(trials.dlt.over[5,], na.rm = TRUE)
    summary.mono.5["#DLT (all)", "Max."] <- max(trials.dlt.all[5,], na.rm = TRUE)
    summary.mono.5["% overdose", "Max."] <- max(trials.perc.over[5,], na.rm = TRUE)
    summary.mono.5["% DLT", "Max."] <- max(trials.perc.dlt[5,], na.rm = TRUE)



    summary.mono.5["#Pat underdose", "2.5%"] <- quantile(trials.n.under[5,], na.rm = TRUE, probs = 0.025)
    summary.mono.5["#Pat target dose", "2.5%"] <- quantile(trials.n.target[5,], na.rm = TRUE, probs = 0.025)
    summary.mono.5["#Pat overdose", "2.5%"] <- quantile(trials.n.over[5,], na.rm = TRUE, probs = 0.025)
    summary.mono.5["#Pat (all)", "2.5%"] <- quantile(trials.n.all[5,], na.rm = TRUE, probs = 0.025)
    summary.mono.5["#DLT underdose", "2.5%"] <- quantile(trials.dlt.under[5,], na.rm = TRUE, probs = 0.025)
    summary.mono.5["#DLT target dose", "2.5%"] <- quantile(trials.dlt.target[5,], na.rm = TRUE, probs = 0.025)
    summary.mono.5["#DLT overdose", "2.5%"] <- quantile(trials.dlt.over[5,], na.rm = TRUE, probs = 0.025)
    summary.mono.5["#DLT (all)", "2.5%"] <- quantile(trials.dlt.all[5,], na.rm = TRUE, probs = 0.025)
    summary.mono.5["% overdose", "2.5%"] <- quantile(trials.perc.over[5,], na.rm = TRUE, probs = 0.025)
    summary.mono.5["% DLT", "2.5%"] <- quantile(trials.perc.dlt[5,], na.rm = TRUE, probs = 0.025)


    summary.mono.5["#Pat underdose", "97.5%"] <- quantile(trials.n.under[5,], na.rm = TRUE, probs = 0.975)
    summary.mono.5["#Pat target dose", "97.5%"] <- quantile(trials.n.target[5,], na.rm = TRUE, probs = 0.975)
    summary.mono.5["#Pat overdose", "97.5%"] <- quantile(trials.n.over[5,], na.rm = TRUE, probs = 0.975)
    summary.mono.5["#Pat (all)", "97.5%"] <- quantile(trials.n.all[5,], na.rm = TRUE, probs = 0.975)
    summary.mono.5["#DLT underdose", "97.5%"] <- quantile(trials.dlt.under[5,], na.rm = TRUE, probs = 0.975)
    summary.mono.5["#DLT target dose", "97.5%"] <- quantile(trials.dlt.target[5,], na.rm = TRUE, probs = 0.975)
    summary.mono.5["#DLT overdose", "97.5%"] <- quantile(trials.dlt.over[5,], na.rm = TRUE, probs = 0.975)
    summary.mono.5["#DLT (all)", "97.5%"] <- quantile(trials.dlt.all[5,], na.rm = TRUE, probs = 0.975)
    summary.mono.5["% overdose", "97.5%"] <- quantile(trials.perc.over[5,], na.rm = TRUE, probs = 0.975)
    summary.mono.5["% DLT", "97.5%"] <- quantile(trials.perc.dlt[5,], na.rm = TRUE, probs = 0.975)

    for(d in 1:length(doses.mono2.b)){
      summ.ppd.m2b[2, d] <- mean(trials.p.dose.m5[d, ], na.rm=T)
      summ.ppd.m2b[3, d] <- median(trials.p.dose.m5[d, ], na.rm=T)
      summ.ppd.m2b[4, d] <- min(trials.p.dose.m5[d, ], na.rm=T)
      summ.ppd.m2b[5, d] <- max(trials.p.dose.m5[d, ], na.rm=T)
    }

    simulation.results[['results mono2.b']] <- results.mono.5
    simulation.results[['summary mono2.b']] <- summary.mono.5
    simulation.results[['MTDs mono2.b']] <- input.mono.5
    simulation.results[['#pat. mono2.b']] <- summ.ppd.m2b
  }


  if(output.sim.config){
    #----------------------------------------------------------------------------------------------
    #if given, write the historical data in the output list
    if(!is.null(historical.data)){
      simulation.results[['historical data']]<- data_matrix_jointBLRM(historical.data)
    }
    #write the parameters of the decision rules in the output list
    simulation.results[['prior']] <- prior_mat_out(m=prior.mu, t=prior.tau)

    simulation.results[["specifications"]] <- input_out_jointBLRM(
      active.mono1.a=active.mono1.a,
      active.mono1.b=active.mono1.b,
      active.mono2.a=active.mono2.a,
      active.mono2.b=active.mono2.b,
      active.combi.a=active.combi.a,
      active.combi.b=active.combi.b,

      dose.ref1 = dose.ref1,
      dose.ref2 = dose.ref2,

      seed = seed,
      dosing.intervals = dosing.intervals,
      saturating = saturating,
      ewoc = ewoc.threshold,
      loss.weights = loss.weights,
      dynamic.weights = dynamic.weights,

      esc.rule = esc.rule,
      esc.comp.max=esc.comp.max,

      cohort.size.mono1.a= cohort.size.mono1.a,
      cohort.prob.mono1.a= cohort.prob.mono1.a,
      cohort.size.mono2.a= cohort.size.mono2.a,
      cohort.prob.mono2.a= cohort.prob.mono2.a,
      cohort.size.mono1.b= cohort.size.mono1.b,
      cohort.prob.mono1.b= cohort.prob.mono1.b,
      cohort.size.mono2.b= cohort.size.mono2.b,
      cohort.prob.mono2.b= cohort.prob.mono2.b,
      cohort.size.combi.a= cohort.size.combi.a,
      cohort.prob.combi.a= cohort.prob.combi.a,
      cohort.size.combi.b= cohort.size.combi.b,
      cohort.prob.combi.b= cohort.prob.combi.b,

      esc.step.mono1.a=esc.step.mono1.a,
      esc.step.mono2.a=esc.step.mono2.a ,
      esc.step.mono1.b=esc.step.mono1.b,
      esc.step.mono2.b=esc.step.mono2.b,
      esc.step.combi.a1=esc.step.combi.a1,
      esc.step.combi.b1=esc.step.combi.b1,
      esc.step.combi.a2=esc.step.combi.a2,
      esc.step.combi.b2=esc.step.combi.b2,
      esc.constrain.mono1.a=esc.constrain.mono1.a,
      esc.constrain.mono2.a=esc.constrain.mono2.a,
      esc.constrain.mono1.b=esc.constrain.mono1.b,
      esc.constrain.mono2.b=esc.constrain.mono2.b,
      esc.constrain.combi.a1=esc.constrain.combi.a1,
      esc.constrain.combi.b1=esc.constrain.combi.b1,
      esc.constrain.combi.a2=esc.constrain.combi.a2,
      esc.constrain.combi.b2=esc.constrain.combi.b2,


      start.dose.mono1.a=start.dose.mono1.a,
      start.dose.mono2.a=start.dose.mono2.a,
      start.dose.mono1.b=start.dose.mono1.b,
      start.dose.mono2.b=start.dose.mono2.b,
      start.dose.combi.a1=start.dose.combi.a1,
      start.dose.combi.a2=start.dose.combi.a2,
      start.dose.combi.b1=start.dose.combi.b1,
      start.dose.combi.b2=start.dose.combi.b2,

      max.n.mono1.a = max.n.mono1.a,
      max.n.mono1.b = max.n.mono1.b,
      max.n.mono2.a = max.n.mono2.a,
      max.n.mono2.b = max.n.mono2.b,
      max.n.combi.a = max.n.combi.a,
      max.n.combi.b = max.n.combi.b,
      decision.combi.a = mtd.decision.combi.a,
      decision.combi.b = mtd.decision.combi.b,
      decision.mono1.a = mtd.decision.mono1.a,
      decision.mono1.b = mtd.decision.mono1.b,
      decision.mono2.a = mtd.decision.mono2.a,
      decision.mono2.b = mtd.decision.mono2.b,


      mtd.enforce.mono1.a = mtd.enforce.mono1.a,
      mtd.enforce.mono1.b = mtd.enforce.mono1.b,
      mtd.enforce.mono2.a = mtd.enforce.mono2.a,
      mtd.enforce.mono2.b = mtd.enforce.mono2.b,
      mtd.enforce.combi.a = mtd.enforce.combi.a,
      mtd.enforce.combi.b = mtd.enforce.combi.b,

      backfill.mono1.a = backfill.mono1.a,
      backfill.mono1.b = backfill.mono1.b,
      backfill.size.mono1.a = backfill.size.mono1.a,
      backfill.size.mono1.b = backfill.size.mono1.b,
      backfill.prob.mono1.a = backfill.prob.mono1.a,
      backfill.prob.mono1.b = backfill.prob.mono1.b,
      backfill.mono2.a = backfill.mono2.a,
      backfill.mono2.b = backfill.mono2.b,
      backfill.size.mono2.a = backfill.size.mono2.a,
      backfill.size.mono2.b = backfill.size.mono2.b,
      backfill.prob.mono2.a = backfill.prob.mono2.a,
      backfill.prob.mono2.b = backfill.prob.mono2.b,
      backfill.combi.a = backfill.combi.a,
      backfill.combi.b = backfill.combi.b,
      backfill.size.combi.a = backfill.size.combi.a,
      backfill.size.combi.b = backfill.size.combi.b,
      backfill.prob.combi.a = backfill.prob.combi.a,
      backfill.prob.combi.b = backfill.prob.combi.b,
      backfill.start.mono1.a = backfill.start.mono1.a,
      backfill.start.mono1.b = backfill.start.mono1.b,
      backfill.start.mono2.a = backfill.start.mono2.a,
      backfill.start.mono2.b = backfill.start.mono2.b,
      backfill.start.combi.a1 = backfill.start.combi.a1,
      backfill.start.combi.a2 = backfill.start.combi.a2,
      backfill.start.combi.b1 = backfill.start.combi.b1,
      backfill.start.combi.b2 = backfill.start.combi.b2
    )

    #write the input parameters and simulation options into the output list.
    options.stan <- as.matrix(c(chains, iter, warmup, adapt_delta, max_treedepth))
    rownames(options.stan) <- c('chains','iter', 'warmup','adapt_delta', 'max_treedepth')
    colnames(options.stan) <- "Value"
    simulation.results[['Stan options']] <-options.stan
  }
  #Working Dir
  if(!is.null(working.path) & !is.null(file.name)){
    if(clean.working.path){
      message("Cleaning up temporary data from working directory.")
      #cleanup of working directory
      file.remove(list.files(working.path, full.names = T, pattern = paste0("^", file.name, "_tmp.")))
    }
  }

  #-------------------------------------
  # Write the results to an excel-file
  #-------------------------------------
  if(!is.null(path) & !is.null(file.name)){
    if(dir.exists(file.path(path))){
      write.xlsx(simulation.results,
                 paste(file.path(path),"/",file.name,'.xlsx', sep=''),
                 colNames=TRUE, rowNames=TRUE, overwrite = TRUE)
    }
  }
  #output
  return(simulation.results)
}
#--------------------------------------------------------------------------------------------
#end of function
#--------------------------------------------------------------------------------------------

