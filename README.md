<!-- badges: start -->
[![R-CMD-check](https://github.com/Boehringer-Ingelheim/decider/workflows/R-CMD-check/badge.svg)](https://github.com/Boehringer-Ingelheim/decider/actions)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# `decider` R package: decision making in multiple-arm oncology dose escalation trials with logistic regression

The `decider` package was created to allow the use of Bayesian logistic 
regression models (BLRMs) for decision making in oncology dose escalation trials 
with multiple arms. It focuses on the use of the meta-analytic ("joint") BLRM 
for mono- and combination therapy, which uses a hierarchical prior structure. 
The package comes with various functions allowing to evaluate the performance 
and operating characteristics of the BLRM and a chosen prior in a given trial
setting using trial simulations or hypothetical data scenarios.

## Background and getting started

The main model implemented in the package is the hierarchical (exchangeable) 
combination therapy BLRM following the model described e.g. by 
[Neuenschwander et al. (2016)](https://doi.org/10.1080/19466315.2016.1174149)
for oncology dose finding in combination therapy. The posterior computations for
the joint BLRM are implemented using the [`rstan`](https://mc-stan.org/) 
package. For a detailed model description, refer to the documentation of 
[`scenario_jointBLRM()`](https://boehringer-ingelheim.github.io/decider/reference/scenario_jointBLRM.html). The `decider `package
allows to use the basic hierarchical combination therapy model with exchangeable 
prior structure together with two variants of the general model:

* Inclusion of a binary covariate

* Use of two different interaction terms following the methodology described  
by [`OncoBayes2`](https://CRAN.R-project.org/package=OncoBayes2), namely, 
the linear interaction term and the saturating interaction term.

The core functions are 
[`scenario_jointBLRM()`](https://boehringer-ingelheim.github.io/decider/reference/scenario_jointBLRM.html) and 
[`sim_jointBLRM()`](https://boehringer-ingelheim.github.io/decider/reference/sim_jointBLRM.html), for evaluation of 
hypothetical data scenarios, respectively trial simulations. The latter function
allows simulations of up to 6 parallel ongoing trial arms, with arbitrarily many
historical trials (mono and/or combination therapy), and provides various 
options for customization -- for instance, cohort size, order of cohort 
enrolment in different arms, dose escalation rules, and further things can be 
adjusted. Additionally, the simulations allow for the inclusion of so-called 
back-fill cohorts, where lower doses are backfilled while escalation is ongoing 
on higher doses.

Please refer to the vignette 
[`Evaluating prior specifications of a joint BLRM for oncology dose finding`](https://boehringer-ingelheim.github.io/decider/articles/intro_jointBLRM.html)
for an introduction to how the functions of the `decider` can be used to conduct
BLRM simulations with multiple arms or evaluate hypothetical data scenarios.

## Installation


The development version can be installed as follows:

```{r}
if(!requireNamespace("remotes", quietly = TRUE)){
  install.packages("remotes")
}
remotes::install_github("https://github.com/Boehringer-Ingelheim/decider")
```

## Documentation

The package documentation is hosted [`here`](https://Boehringer-Ingelheim.github.io/decider/).
