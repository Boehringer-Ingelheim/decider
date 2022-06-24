#' The 'OncoBLRM' package.
#'
#' @description Performs Bayesian logistic regression for binary dose-toxicity
#' data from monotherapy and/or (two-component) combination therapy dose-finding trials.
#' The underlying model is referred to as Bayesian logistic regression model (BLRM)
#' and specified according to Neuenschwander et al. (2008, 2014, 2016).
#' All Bayesian models were implemented in the
#' in the \href{https://mc-stan.org}{Stan} modeling language (Stan Development Team (2020))
#' using the \code{\link[rstan:rstan]{rstan-package}}
#' and the \code{\link[rstantools:rstantools-package]{rstantools-package}}.
#'
#' Currently, only the so-called joint BLRM is included in the package. Different
#' methods for the derivation of dosing recommendations are supported, among others the
#' escalation with overdose control (EWOC) criterion that goes back to Babb et al. (1998).
#' The main functions are \code{\link[OncoBLRM:sim_jointBLRM]{sim_jointBLRM}()} and
#' \code{\link[OncoBLRM:scenario_jointBLRM]{scenario_jointBLRM}()}.
#' Refer to their documentation entries for a detailed description of the underlying methods.
#'
#' The methods implemented by this package were mainly developed in the context
#' of oncology dose-finding trials, but can also be applied to different therapeutic
#' areas provided that the methodology of estimating a monotonically increasing
#' dose-toxicity relationship based on binary safety data applies.
#' For the documentation and the argument names used by the package, the typical
#' terminology of oncology trials is used. That is, the dose to be determined
#' by the trials is referred to as the maximum tolerated dose (MTD), which is
#' determined based on binary safety data given by so-called dose-limiting
#' toxicities (DLT). These DLTs are a pre-specified set of adverse events that
#' are considered to be dose limiting. In this context, the goal of a trial is
#' to determine the MTD so that it has a true DLT rate in some pre-specified
#' dosing interval. Note that other therapeutic areas often use slightly
#' different notions instead.
#'
#' @docType package
#' @name OncoBLRM-package
#' @aliases OncoBLRM
#' @useDynLib OncoBLRM, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import ggplot2
#' @import openxlsx
#' @import flock
#' @import doRNG
#' @importFrom rstan sampling
#' @importFrom stats rbinom sd family quantile median
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach registerDoSEQ "%dopar%" getDoParRegistered
#' @importFrom doParallel registerDoParallel
#' @importFrom gridExtra grid.arrange
#' @importFrom scales breaks_extended
#' @importFrom grDevices cairo_pdf
#' @importFrom utils write.table
#' @importFrom openssl sha256
#'
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
#' Babb, J., Rogatko, A., & Zacks, S. (1998). Cancer phase I clinical trials: Efficient dose escalation with overdose control.
#' Statistics in medicine 17(10), 1103-1120.
#'
#' Zhou, H.,  Yuan, Y., & Nie, L. (2018). Accuracy, safety, and reliability of novel phase I designs.
#' Clinical Cancer Research, 24(21), 5483-5484 <doi: 10.1158/1078-0432.ccr-18-0168>.
#'
#' @seealso \code{\link[OncoBLRM:scenario_jointBLRM]{scenario_jointBLRM}()},
#' \code{\link[OncoBLRM:sim_jointBLRM]{sim_jointBLRM}()},\code{\link[OncoBLRM:fit_jointBLRM]{fit_jointBLRM}()},
#' \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling}},
#' \code{\link[rstan:rstan]{rstan-package}}.
NULL
