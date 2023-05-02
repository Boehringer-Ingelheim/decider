# decider R package -- development version
This package was created to allow the use of Bayesian logistic regression models 
(BLRMs) for decision making in oncology dose escalation trials with multiple 
arms. It focuses on the use of the meta-analytic ("joint") BLRM for mono- and
combination therapy, which uses a hierarchical prior structure. The package
comes with various functions allowing to evaluate the BLRM and a chosen prior
using simulations or hypothetical data scenarios.

Please note that the present version is an early development version of the 
package.

The main functionality is using a hierarchical (exchangeable) combination therapy BLRM, 
which is implemented using the \code{rstan} package. The package allows some variants 
of the general model:
* Inclusion of a binary covariate
* Use of two different interaction terms, namely, using the terminology introduced
in \ref{OncoBayes2}, the linear interaction term and the saturating interaction term.

The core functions are \code{scenario_jointBLRM} and \code{sim_jointBLRM}, for 
evaluation of hypothetical data scenarios, respectively trial simulations. The
latter function allows simulations of up to 6 parallel ongoing trial arms, with
arbitrarily many historical trials (mono and/or combination therapy), and provides
various options for customization -- for instance, cohort size, order of cohort
enrolment in different arms, dose escalation rules, and further things can be adjusted. 
Additionally, the simulations allow for the inclusion of so-called back-fill cohorts, 
where lower doses are backfilled while escalation is ongoing on higher doses.

Please refer to the vignette [...] for an introduction to how the functions of 
the \code{decider} can be used to conduct BLRM simulations with multiple arms
or evaluate hypothetical data scenarios.
