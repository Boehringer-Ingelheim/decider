
#-------------------------------------------------------------------------------
#Internal functions that perform MCMC sampling in higher-level functions
#(used for calls to stan during scenario_jointBLRM and sim_jointBLRM)
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
#sampling_jointBLRM
#-------------------------------------------------------------------------------
#internal function used by scenario_jointBLRM:
#calls sampling method for stanmodel "JointBLRM" with given data,
#returns a matrix of the samples (does not return a stanmodel object).

#' @keywords internal
sampling_jointBLRM <- function(
  dose1,
  dose2,
  dose.ref1,
  dose.ref2,
  n.pat,
  n.dlt,
  n.study,
  MAP.prior = FALSE,
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
  iter=10000,
  warmup = floor(iter/2),
  refresh = floor(iter/10),
  adapt_delta = 0.95,
  max_treedepth = 15,
  chains = 4,
  seed=sample.int(.Machine$integer.max, 1)
){
  nobs <- length(dose1)
  nstd <- length(levels(factor(n.study)))

  return(as.matrix(sampling(
    object = stanmodels$jointBLRM,
    data = list(
      n_obs = nobs,
      n_studies = nstd,
      dose_1 = as.array(dose1),
      dose_2 = as.array(dose2),
      dose_c = c(dose.ref1, dose.ref2),
      n = as.array(n.pat),
      r = as.array(n.dlt),
      s = as.array(n.study),
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
    iter=iter, warmup=warmup, refresh=refresh, chains=chains, seed=seed,
    control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth),
    pars="log_ab")))
}


#' @keywords internal
sampling_covariate_jointBLRM <- function(
  dose1,
  dose2,
  dose.ref1,
  dose.ref2,
  n.pat,
  n.dlt,
  n.study,
  cov,
  MAP.prior = FALSE,
  saturating = FALSE,
  prior.mu = list(mu_a1 =  c(logit(0.33), 2),
                  mu_b1 =  c(0,          1),
                  mu_a2 =  c(logit(0.33), 2),
                  mu_b2 =  c(0,          1),
                  mu_eta = c(0,          1.121),

                  mu_c1 = c(0,  1),   #NEW
                  mu_c2 = c(0,  1)    #NEW
  ),
  prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
                   tau_b1 =  c(log(0.125), log(2)/1.96),
                   tau_a2 =  c(log(0.25),  log(2)/1.96),
                   tau_b2 =  c(log(0.125), log(2)/1.96),
                   tau_eta = c(log(0.125), log(2)/1.96),

                   tau_c1 = c(log(0.125), log(2)/1.96),   #NEW
                   tau_c2 = c(log(0.125), log(2)/1.96)    #NEW
  ),

  two_sided1 = TRUE,   #NEW
  two_sided2 = TRUE,   #NEW

  iter=10000,
  warmup = floor(iter/2),
  refresh = floor(iter/10),
  adapt_delta = 0.95,
  max_treedepth = 15,
  chains = 4,
  seed=sample.int(.Machine$integer.max, 1)
){
  nobs <- length(dose1)
  nstd <- length(levels(factor(n.study)))

  return(as.matrix(sampling(
    object = stanmodels$jointBLRMcov,
    data = list(
      n_obs = nobs,
      n_studies = nstd,
      dose_1 = as.array(dose1),
      dose_2 = as.array(dose2),
      dose_c = c(dose.ref1, dose.ref2),
      n = as.array(n.pat),
      r = as.array(n.dlt),
      c = as.array(cov),
      s = as.array(n.study),
      doMAP = ifelse(MAP.prior, 1, 0),
      saturating = ifelse(saturating, 1, 0),

      twoside1 = ifelse(two_sided1, 1, 0),
      twoside2 = ifelse(two_sided2, 1, 0),


      mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
                  prior.mu$mu_a2[1], prior.mu$mu_b2[1], prior.mu$mu_eta[1],
                  prior.mu$mu_c1[1], prior.mu$mu_c2[1]),
      sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
                prior.mu$mu_a2[2], prior.mu$mu_b2[2], prior.mu$mu_eta[2],
                prior.mu$mu_c1[2], prior.mu$mu_c2[2]),
      mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
                   prior.tau$tau_a2[1], prior.tau$tau_b2[1], prior.tau$tau_eta[1],
                   prior.tau$tau_c1[1], prior.tau$tau_c2[1]),
      sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
                 prior.tau$tau_a2[2], prior.tau$tau_b2[2], prior.tau$tau_eta[2],
                 prior.tau$tau_c1[2], prior.tau$tau_c2[2])
    ),
    iter=iter, warmup=warmup, refresh=refresh, chains=chains, seed=seed,
    control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth),
    pars="log_ab")))
}


# @keywords internal
# sampling_jointBLRMsat <- function(
#   dose1,
#   dose2,
#   dose.ref1,
#   dose.ref2,
#   n.pat,
#   n.dlt,
#   n.study,
#   MAP.prior = FALSE,
#   prior.mu = list(mu_a1 =  c(logit(0.33), 2),
#                   mu_b1 =  c(0,          1),
#                   mu_a2 =  c(logit(0.33), 2),
#                   mu_b2 =  c(0,          1),
#                   mu_eta = c(0,          1.121)
#   ),
#   prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
#                    tau_b1 =  c(log(0.125), log(2)/1.96),
#                    tau_a2 =  c(log(0.25),  log(2)/1.96),
#                    tau_b2 =  c(log(0.125), log(2)/1.96),
#                    tau_eta = c(log(0.125), log(2)/1.96)
#   ),
#   iter=10000,
#   warmup = floor(iter/2),
#   refresh = floor(iter/10),
#   adapt_delta = 0.95,
#   max_treedepth = 15,
#   chains = 4,
#   seed=sample.int(.Machine$integer.max, 1)
# ){
#   nobs <- length(dose1)
#   nstd <- length(levels(factor(n.study)))
#
#   return(as.matrix(sampling(
#     object = stanmodels$jointBLRMsaturating,
#     data = list(
#       n_obs = nobs,
#       n_studies = nstd,
#       dose_1 = as.array(dose1),
#       dose_2 = as.array(dose2),
#       dose_c = c(dose.ref1, dose.ref2),
#       n = as.array(n.pat),
#       r = as.array(n.dlt),
#       s = as.array(n.study),
#       doMAP = ifelse(MAP.prior, 1, 0),
#       mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
#                   prior.mu$mu_a2[1], prior.mu$mu_b2[1], prior.mu$mu_eta[1]),
#       sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
#                 prior.mu$mu_a2[2], prior.mu$mu_b2[2], prior.mu$mu_eta[2]),
#       mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
#                    prior.tau$tau_a2[1], prior.tau$tau_b2[1], prior.tau$tau_eta[1]),
#       sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
#                  prior.tau$tau_a2[2], prior.tau$tau_b2[2], prior.tau$tau_eta[2])
#     ),
#     iter=iter, warmup=warmup, refresh=refresh, chains=chains, seed=seed,
#     control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth),
#     pars="log_ab")))
# }

#-------------------------------------------------------------------------------
#post_tox_jointBLRM
#-------------------------------------------------------------------------------

#function that determines posterior distribution of DLT rate for multiple trials
#and multiple (potentially not tested) dose levels.
#Used by scenario_jointBLRM internally.
#Returns a matrix with interval probabilities and summary statistics of the
#posterior DLT rates for each given trial. Matrices are pre-formatted so that
#it matches the structure returned by scenario_jointBLRM.

#' @keywords internal
post_tox_jointBLRM <- function(
  study.interest,
  type.interest,
  dose1.interest,
  dose2.interest,
  dose1,
  dose2,
  dose.ref1,
  dose.ref2,
  n.pat,
  n.dlt,
  n.study,
  MAP.prior=FALSE,
  saturating = FALSE,
  dosing.intervals = c(0.16, 0.33),
  probs = c(0.025, 0.5, 0.975),
  prior.mu = list(mu_a1 =  c(-1.386294, 2),
                  mu_b1 =  c(0,         1),
                  mu_a2 =  c(-1.386294, 2),
                  mu_b2 =  c(0,         1),
                  mu_eta = c(0,         1.121)
  ),
  prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
                   tau_b1 =  c(log(0.125), log(2)/1.96),
                   tau_a2 =  c(log(0.25),  log(2)/1.96),
                   tau_b2 =  c(log(0.125), log(2)/1.96),
                   tau_eta = c(log(0.125), log(2)/1.96)
  ),
  iter = 10000,
  adapt_delta = 0.95,
  max_treedepth = 15,
  warmup = floor(iter/2),
  refresh = floor(iter/10),
  chains = 4,
  seed=sample.int(.Machine$integer.max, 1)
)
{

    samples <- sampling_jointBLRM(dose1=dose1,
                                     dose2=dose2,
                                     dose.ref1 = dose.ref1,
                                     dose.ref2 = dose.ref2,
                                     n.pat = n.pat,
                                     n.dlt = n.dlt,
                                     n.study = n.study,
                                     prior.mu = prior.mu,
                                     prior.tau = prior.tau,
                                     MAP.prior = MAP.prior,
                                     saturating = saturating,
                                     iter=iter,
                                     warmup=warmup,
                                     refresh=refresh,
                                     chains=chains,
                                     adapt_delta=adapt_delta,
                                     max_treedepth=max_treedepth,
                                     seed=seed
    )


  nstudies <- length(study.interest)
  ndose <- length(dose1.interest)
  nsamp <- length(samples[,1])
  dose1_ <- dose1.interest/dose.ref1
  dose2_ <- dose2.interest/dose.ref2

  res_tox <- list()
  resref_tox <- list()

  ninters <- length(dosing.intervals)+1
  nquan <- length(probs)
  name.inters <- rep(0, ninters)
  for(i in 1:ninters){
    if(i==1){
      name.inters[i] <- paste0("P([0,", dosing.intervals[1], "))")
    }else if(i==ninters){
      name.inters[i] <- paste0("P([",dosing.intervals[ninters-1], ",1])")
    }else{
      name.inters[i] <- paste0("P([",dosing.intervals[i-1],",",dosing.intervals[i],"))")
    }
  }

  for(studycurr in 1:nstudies){
    study <- study.interest[studycurr]
    type <- type.interest[studycurr]
    #determine type of interest and adjust doses accordingly
    if(type=="mono1"){
      inds <- which(dose1_>0 & dose2_==0)
      dose1_c <- dose1_[inds]
      dose2_c <- dose2_[inds]
      ndosec <- length(dose1_c)
      dose1.interestc <- dose1.interest[inds]
      dose2.interestc <- dose2.interest[inds]
    }else if(type=="mono2"){
      inds <- which(dose2_>0 & dose1_==0)
      dose1_c <- dose1_[inds]
      dose2_c <- dose2_[inds]
      ndosec <- length(dose1_c)
      dose1.interestc <- dose1.interest[inds]
      dose2.interestc <- dose2.interest[inds]
    }else if(type=="combi"){
      inds <- which(dose1_>0 & dose2_>0)
      dose1_c <- dose1_[inds]
      dose2_c <- dose2_[inds]
      ndosec <- length(dose1_c)
      dose1.interestc <- dose1.interest[inds]
      dose2.interestc <- dose2.interest[inds]
    }else{
      dose1_c <- dose1_
      dose2_c <- dose2_
      ndosec <- ndose
      dose1.interestc <- dose1.interest
      dose2.interestc <- dose2.interest
    }

    logdose1_c <- log(dose1_c)
    logdose2_c <- log(dose2_c)

    toxs <- matrix(NA, nrow = ndosec, ncol = nsamp)
    toxs1 <- rep(0, nsamp)
    toxs2 <- rep(0, nsamp)
    toxs0 <- rep(0, nsamp)
    #get study-specific posterior samples
    la1 <- samples[, paste0("log_ab[", study, ",1]")]
    lb1 <- samples[, paste0("log_ab[", study, ",2]")]
    la2 <- samples[, paste0("log_ab[", study, ",3]")]
    lb2 <- samples[, paste0("log_ab[", study, ",4]")]
    eta <- samples[, paste0("log_ab[", study, ",5]")]
    explb1 <- exp(lb1)
    explb2 <- exp(lb2)

    for(d in 1:ndosec){
      if(dose2_c[d]==0 & dose1_c[d]>0){
        #treatment mono1
        toxs[d, ] <- inv_logit(la1 + explb1*logdose1_c[d])
      }else if(dose1_c[d]==0 & dose2_c[d]>0){
        #treatment mono2
        toxs[d, ] <- inv_logit(la2 + explb2*logdose2_c[d])
      }else if(dose2_c[d]>0 & dose1_c[d]>0){
        #treatment combination
        toxs1 <- inv_logit(la1 + explb1*logdose1_c[d])
        toxs2 <- inv_logit(la2 + explb2*logdose2_c[d])
        toxs0 <- toxs1 + toxs2 - toxs1*toxs2
        if(saturating){
          toxs[d, ] <- inv_logit(logit(toxs0) +
                 2*eta*(dose1_c[d]*dose2_c[d])/(1 + dose1_c[d]*dose2_c[d]))
        }else{
          toxs[d, ] <- inv_logit(logit(toxs0) + eta*dose1_c[d]*dose2_c[d])
        }
      }else{
        #no treatment, put in 0
        toxs[d, ] <- rep(0, times = nsamp)
      }
    }

    n.entries <- 2+nquan+ninters
    if(nquan>0){
      colnames.res <- c("mean", "sd", paste0("q.", round(probs*100, digits=2), "%"), name.inters)
    }else{
      colnames.res <- c("mean", "sd", name.inters)
    }
    res_curr <- matrix(NA, nrow = ndosec, ncol = n.entries)
    rownames(res_curr) <- paste0(dose1.interestc, "+", dose2.interestc)
    colnames(res_curr) <-colnames.res

    for(d in 1:ndosec){
      res_curr[d, 1] <- mean(toxs[d, ])
      res_curr[d, 2] <- sd(toxs[d, ])
      if(nquan>0){
        res_curr[d, 3:(2+nquan)] <- quantile(toxs[d, ], probs = probs)
      }
      for(inter in 1:ninters){
        if(inter==1){
          res_curr[d, 2+nquan+1] <- mean(ifelse(toxs[d,]<dosing.intervals[1], 1, 0))
        }else if(inter==ninters){
          res_curr[d, 2+nquan+ninters] <- mean(ifelse(toxs[d,]>=dosing.intervals[ninters-1], 1, 0))
        }else{
          res_curr[d, 2+nquan+inter] <- mean(ifelse(toxs[d,]>=dosing.intervals[inter-1] &
                                                      toxs[d,]<dosing.intervals[inter],
                                                    1, 0))
        }
      }
    }
    res_tox[[paste0(studycurr)]] <- res_curr


    #compute samples from reference toxicities for study
    #(only required for dynamic loss escalation)
    ref_toxs <- matrix(NA, nrow = 3, ncol = nsamp)
    rtoxs0 <- rep(0, nsamp)
    ref_toxs[1, ] <- inv_logit(la1)
    ref_toxs[2, ] <- inv_logit(la2)
    logit_rtoxs0 <- logit(ref_toxs[1,]+ref_toxs[2,]-ref_toxs[1,]*ref_toxs[2,])
    ref_toxs[3, ] <- inv_logit(logit_rtoxs0 + eta)

    interprobs_ref <- matrix(NA, nrow=3, ncol=ninters)
    colnames(interprobs_ref) <- name.inters
    rownames(interprobs_ref) <- c("Ref.Mono1", "Ref.Mono2", "Ref.Combi")
    for(d in 1:3){
      for(inter in 1:ninters){
        if(inter==1){
          interprobs_ref[d, 1] <- mean(ifelse(ref_toxs[d,]<dosing.intervals[1], 1, 0))
        }else if(inter==ninters){
          interprobs_ref[d, ninters] <- mean(ifelse(ref_toxs[d,]>=dosing.intervals[ninters-1], 1, 0))
        }else{
          interprobs_ref[d, inter] <- mean(ifelse(ref_toxs[d,]>=dosing.intervals[inter-1] &
                                                    ref_toxs[d,]<dosing.intervals[inter],
                                                  1, 0))
        }
      }
    }

    resref_tox[[paste0(studycurr)]]<-interprobs_ref

  }

  return(list(
    "postsummaries"=res_tox,
    "refprobs"=resref_tox
  ))

}


#' @keywords internal
post_tox_covariate_jointBLRM <- function(
  study.interest,
  type.interest,

  cov.interest,                  #'NEW': 0, 1 (otherwise, do both)

  dose1.interest,
  dose2.interest,
  dose1,
  dose2,
  dose.ref1,
  dose.ref2,
  n.pat,
  n.dlt,
  n.study,
  cov,                                #'NEW'
  MAP.prior=FALSE,
  saturating = FALSE,
  dosing.intervals = c(0.16, 0.33),
  probs = c(0.025, 0.5, 0.975),
  prior.mu = list(mu_a1 =  c(-1.386294, 2),
                  mu_b1 =  c(0,         1),
                  mu_a2 =  c(-1.386294, 2),
                  mu_b2 =  c(0,         1),
                  mu_eta = c(0,         1.121),

                  mu_c1 = c(0,  1),   #NEW
                  mu_c2 = c(0,  1)    #NEW
  ),
  prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
                   tau_b1 =  c(log(0.125), log(2)/1.96),
                   tau_a2 =  c(log(0.25),  log(2)/1.96),
                   tau_b2 =  c(log(0.125), log(2)/1.96),
                   tau_eta = c(log(0.125), log(2)/1.96),

                   tau_c1 = c(log(0.125), log(2)/1.96),   #NEW
                   tau_c2 = c(log(0.125), log(2)/1.96)    #NEW
  ),

  two_sided1 = TRUE,                             #'NEW'
  two_sided2 = TRUE,                             #'NEW'

  iter = 10000,
  adapt_delta = 0.95,
  max_treedepth = 15,
  warmup = floor(iter/2),
  refresh = floor(iter/10),
  chains = 4,
  seed=sample.int(.Machine$integer.max, 1)
)
{

  samples <- sampling_covariate_jointBLRM(dose1=dose1,
                                dose2=dose2,
                                dose.ref1 = dose.ref1,
                                dose.ref2 = dose.ref2,
                                n.pat = n.pat,
                                n.dlt = n.dlt,
                                cov = cov,
                                n.study = n.study,
                                prior.mu = prior.mu,
                                prior.tau = prior.tau,
                                two_sided1 = two_sided1,
                                two_sided2 = two_sided2,
                                saturating = saturating,
                                MAP.prior = MAP.prior,
                                iter=iter,
                                warmup=warmup,
                                refresh=refresh,
                                chains=chains,
                                adapt_delta=adapt_delta,
                                max_treedepth=max_treedepth,
                                seed=seed
  )


  nstudies <- length(study.interest)
  ndose <- length(dose1.interest)
  nsamp <- length(samples[,1])
  dose1_ <- dose1.interest/dose.ref1
  dose2_ <- dose2.interest/dose.ref2

  res_tox <- list()
  resref_tox <- list()

  ninters <- length(dosing.intervals)+1
  nquan <- length(probs)
  name.inters <- rep(0, ninters)
  for(i in 1:ninters){
    if(i==1){
      name.inters[i] <- paste0("P([0,", dosing.intervals[1], "))")
    }else if(i==ninters){
      name.inters[i] <- paste0("P([",dosing.intervals[ninters-1], ",1])")
    }else{
      name.inters[i] <- paste0("P([",dosing.intervals[i-1],",",dosing.intervals[i],"))")
    }
  }

  for(studycurr in 1:nstudies){
    study <- study.interest[studycurr]
    type <- type.interest[studycurr]
    covcurr <- cov.interest[studycurr]

    if(covcurr == 0){
      #determine type of interest and adjust doses accordingly
      if(type=="mono1"){
        inds <- which(dose1_>0 & dose2_==0)
        dose1_c <- dose1_[inds]
        dose2_c <- dose2_[inds]
        ndosec <- length(dose1_c)
        dose1.interestc <- dose1.interest[inds]
        dose2.interestc <- dose2.interest[inds]
      }else if(type=="mono2"){
        inds <- which(dose2_>0 & dose1_==0)
        dose1_c <- dose1_[inds]
        dose2_c <- dose2_[inds]
        ndosec <- length(dose1_c)
        dose1.interestc <- dose1.interest[inds]
        dose2.interestc <- dose2.interest[inds]
      }else if(type=="combi"){
        inds <- which(dose1_>0 & dose2_>0)
        dose1_c <- dose1_[inds]
        dose2_c <- dose2_[inds]
        ndosec <- length(dose1_c)
        dose1.interestc <- dose1.interest[inds]
        dose2.interestc <- dose2.interest[inds]
      }else{
        dose1_c <- dose1_
        dose2_c <- dose2_
        ndosec <- ndose
        dose1.interestc <- dose1.interest
        dose2.interestc <- dose2.interest
      }
      cov_c <- rep(0, times = ndosec)

    }else if(covcurr == 1){
      #determine type of interest and adjust doses accordingly
      if(type=="mono1"){
        inds <- which(dose1_>0 & dose2_==0)
        dose1_c <- dose1_[inds]
        dose2_c <- dose2_[inds]
        ndosec <- length(dose1_c)
        dose1.interestc <- dose1.interest[inds]
        dose2.interestc <- dose2.interest[inds]
      }else if(type=="mono2"){
        inds <- which(dose2_>0 & dose1_==0)
        dose1_c <- dose1_[inds]
        dose2_c <- dose2_[inds]
        ndosec <- length(dose1_c)
        dose1.interestc <- dose1.interest[inds]
        dose2.interestc <- dose2.interest[inds]
      }else if(type=="combi"){
        inds <- which(dose1_>0 & dose2_>0)
        dose1_c <- dose1_[inds]
        dose2_c <- dose2_[inds]
        ndosec <- length(dose1_c)
        dose1.interestc <- dose1.interest[inds]
        dose2.interestc <- dose2.interest[inds]
      }else{
        dose1_c <- dose1_
        dose2_c <- dose2_
        ndosec <- ndose
        dose1.interestc <- dose1.interest
        dose2.interestc <- dose2.interest
      }
      cov_c <- rep(1, times = ndosec)

    }else{
      #determine type of interest and adjust doses accordingly
      if(type=="mono1"){
        inds <- which(dose1_>0 & dose2_==0)
        nunique <- length(inds)
        dose1_c <- c(dose1_[inds], dose1_[inds])
        dose2_c <- c(dose2_[inds], dose2_[inds])
        ndosec <- length(dose1_c)
        dose1.interestc <- c(dose1.interest[inds], dose1.interest[inds])
        dose2.interestc <- c(dose2.interest[inds], dose2.interest[inds])
      }else if(type=="mono2"){
        inds <- which(dose2_>0 & dose1_==0)
        nunique <- length(inds)
        dose1_c <- c(dose1_[inds], dose1_[inds])
        dose2_c <- c(dose2_[inds], dose2_[inds])
        ndosec <- length(dose1_c)
        dose1.interestc <- c(dose1.interest[inds], dose1.interest[inds])
        dose2.interestc <- c(dose2.interest[inds], dose2.interest[inds])
      }else if(type=="combi"){
        inds <- which(dose1_>0 & dose2_>0)
        nunique <- length(inds)
        dose1_c <- c(dose1_[inds], dose1_[inds])
        dose2_c <- c(dose2_[inds], dose2_[inds])
        ndosec <- length(dose1_c)
        dose1.interestc <- c(dose1.interest[inds], dose1.interest[inds])
        dose2.interestc <- c(dose2.interest[inds], dose2.interest[inds])
      }else{
        nunique <- length(dose1_)
        dose1_c <- c(dose1_, dose1_)
        dose2_c <- c(dose2_, dose2_)
        ndosec <- length(dose1_c)
        dose1.interestc <- c(dose1.interest, dose1.interest)
        dose2.interestc <- c(dose2.interest, dose2.interest)
      }
      cov_c <- c(rep(0, times = nunique), rep(1, times = nunique))

    }

    logdose1_c <- log(dose1_c)
    logdose2_c <- log(dose2_c)

    toxs <- matrix(NA, nrow = ndosec, ncol = nsamp)
    toxs1 <- rep(0, nsamp)
    toxs2 <- rep(0, nsamp)
    toxs0 <- rep(0, nsamp)
    #get study-specific posterior samples
    la1 <- samples[, paste0("log_ab[", study, ",1]")]
    lb1 <- samples[, paste0("log_ab[", study, ",2]")]
    la2 <- samples[, paste0("log_ab[", study, ",3]")]
    lb2 <- samples[, paste0("log_ab[", study, ",4]")]
    eta <- samples[, paste0("log_ab[", study, ",5]")]
    if(two_sided1){
      c1 <- samples[, paste0("log_ab[", study, ",6]")]
    }else{
      c1 <- exp(samples[, paste0("log_ab[", study, ",6]")])
    }

    if(two_sided2){
      c2 <- samples[, paste0("log_ab[", study, ",7]")]
    }else{
      c2 <- exp(samples[, paste0("log_ab[", study, ",7]")])
    }

    explb1 <- exp(lb1)
    explb2 <- exp(lb2)

    for(d in 1:ndosec){
      if(dose2_c[d]==0 & dose1_c[d]>0){
        #treatment mono1
        if(cov_c[d]==0){
          toxs[d, ] <- inv_logit(la1 + explb1*logdose1_c[d])
        }else{
          toxs[d, ] <- inv_logit(la1 + c1 + explb1*logdose1_c[d])
        }
      }else if(dose1_c[d]==0 & dose2_c[d]>0){
        #treatment mono2
        if(cov_c[d]==0){
          toxs[d, ] <- inv_logit(la2 + explb2*logdose2_c[d])
        }else{
          toxs[d, ] <- inv_logit(la2 + c2 + explb2*logdose2_c[d])
        }
      }else if(dose2_c[d]>0 & dose1_c[d]>0){
        #treatment combination
        if(cov_c[d]==0){
          toxs1 <- inv_logit(la1 + explb1*logdose1_c[d])
          toxs2 <- inv_logit(la2 + explb2*logdose2_c[d])
        }else{
          toxs1 <- inv_logit(la1 + c1 + explb1*logdose1_c[d])
          toxs2 <- inv_logit(la2 + c2 + explb2*logdose2_c[d])
        }
        toxs0 <- toxs1 + toxs2 - toxs1*toxs2
        if(saturating){
          toxs[d, ] <- inv_logit(logit(toxs0) +
                                   2*eta*(dose1_c[d]*dose2_c[d])/(1 + dose1_c[d]*dose2_c[d]))
        }else{
          toxs[d, ] <- inv_logit(logit(toxs0) + eta*dose1_c[d]*dose2_c[d])
        }
      }else{
        #no treatment, put in 0
        toxs[d, ] <- rep(0, times = nsamp)
      }
    }

    #additional entry for cov (0 or 1)
    n.entries <- 2+nquan+ninters + 1
    if(nquan>0){
      colnames.res <- c("mean", "sd", paste0("q.", round(probs*100, digits=2), "%"), name.inters, "covariate")
    }else{
      colnames.res <- c("mean", "sd", name.inters, "covariate")
    }
    res_curr <- matrix(NA, nrow = ndosec, ncol = n.entries)
    rownames(res_curr) <- paste0(dose1.interestc, "+", dose2.interestc)
    colnames(res_curr) <-colnames.res
    res_curr[, n.entries] <- cov_c

    for(d in 1:ndosec){
      res_curr[d, 1] <- mean(toxs[d, ])
      res_curr[d, 2] <- sd(toxs[d, ])
      if(nquan>0){
        res_curr[d, 3:(2+nquan)] <- quantile(toxs[d, ], probs = probs)
      }
      for(inter in 1:ninters){
        if(inter==1){
          res_curr[d, 2+nquan+1] <- mean(ifelse(toxs[d,]<dosing.intervals[1], 1, 0))
        }else if(inter==ninters){
          res_curr[d, 2+nquan+ninters] <- mean(ifelse(toxs[d,]>=dosing.intervals[ninters-1], 1, 0))
        }else{
          res_curr[d, 2+nquan+inter] <- mean(ifelse(toxs[d,]>=dosing.intervals[inter-1] &
                                                      toxs[d,]<dosing.intervals[inter],
                                                    1, 0))
        }
      }
    }
    res_tox[[paste0(studycurr)]] <- res_curr


    #compute samples from reference toxicities for study
    #(only required for dynamic loss escalation)
    ref_toxs <- matrix(NA, nrow = 6, ncol = nsamp)
    #rtoxs0 <- rep(0, nsamp)
    ref_toxs[1, ] <- inv_logit(la1)
    ref_toxs[2, ] <- inv_logit(la2)
    logit_rtoxs0 <- logit(ref_toxs[1,]+ref_toxs[2,]-ref_toxs[1,]*ref_toxs[2,])
    ref_toxs[3, ] <- inv_logit(logit_rtoxs0 + eta)

    #with covariate
    ref_toxs[4, ] <- inv_logit(la1 + c1)
    ref_toxs[5, ] <- inv_logit(la2 + c2)
    logit_rtoxs0c <- logit(ref_toxs[4,]+ref_toxs[5,]-ref_toxs[4,]*ref_toxs[5,])
    ref_toxs[6, ] <- inv_logit(logit_rtoxs0c + eta)

    interprobs_ref <- matrix(NA, nrow=6, ncol=ninters)
    colnames(interprobs_ref) <- name.inters
    rownames(interprobs_ref) <- c("Ref.Mono1", "Ref.Mono2", "Ref.Combi", "Ref.Mono1cov", "Ref.Mono2cov", "Ref.Combicov")
    for(d in 1:6){
      for(inter in 1:ninters){
        if(inter==1){
          interprobs_ref[d, 1] <- mean(ifelse(ref_toxs[d,]<dosing.intervals[1], 1, 0))
        }else if(inter==ninters){
          interprobs_ref[d, ninters] <- mean(ifelse(ref_toxs[d,]>=dosing.intervals[ninters-1], 1, 0))
        }else{
          interprobs_ref[d, inter] <- mean(ifelse(ref_toxs[d,]>=dosing.intervals[inter-1] &
                                                    ref_toxs[d,]<dosing.intervals[inter],
                                                  1, 0))
        }
      }
    }

    resref_tox[[paste0(studycurr)]]<-interprobs_ref

  }

  return(list(
    "postsummaries"=res_tox,
    "refprobs"=resref_tox
  ))

}



#Runtime optimization:
# generates hash value from compressed data set which is used
# to create file that stores MCMC results with unique (hashed)
# name for each data scenario. Using this, previously conducted
# MCMC results can be reloaded.
# The implementation uses the flock package to ensure that reloading
# can be carried out across different processes during parallelized
# simulation runs. Prior to saving or loading an MCMC result,
# the function will check a lock-file to see whether some
# other process is currently using the file. In this case,
# the process will wait until the file is unlocked.
# Unlocking commands will be carried out via the on.exit()-function,
# to ensure that the file does not retain locked in case of
# errors or manual stops by the user while a file is locked ("deadlock").
#

#' @keywords internal
post_tox_jointBLRM_sim <- function(
  study.interest,
  type.interest,
  dose1.interest,
  dose2.interest,
  BLRM,
  file.name,
  working.path,
  dose1,
  dose2,
  dose.ref1,
  dose.ref2,
  n.pat,
  n.dlt,
  n.study,
  d.loss = FALSE,
  dosing.intervals = c(0.16, 0.33),
  prior.mu = list(mu_a1 =  c(-1.386294, 2),
                  mu_b1 =  c(0,         1),
                  mu_a2 =  c(-1.386294, 2),
                  mu_b2 =  c(0,         1),
                  mu_eta = c(0,         1.121)
  ),
  prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
                   tau_b1 =  c(log(0.125), log(2)/1.96),
                   tau_a2 =  c(log(0.25),  log(2)/1.96),
                   tau_b2 =  c(log(0.125), log(2)/1.96),
                   tau_eta = c(log(0.125), log(2)/1.96)
  ),
  saturating = FALSE,
  iter = 10000,
  adapt_delta = 0.95,
  max_treedepth = 15,
  warmup = floor(iter/2),
  refresh = floor(iter/10),
  chains = 4,
  #n_tries = 10,
  seed=sample.int(.Machine$integer.max, 1)
)
{

  #determine how many observations remain after collapsing data
  procs <- rep(F, length(dose1))
  cnt <- 0
  for(i in 1:length(dose1)){
    if(!procs[i]){
      idx <- which(dose1[i]==dose1 &
                     dose2[i]==dose2 &
                     n.study[i] == n.study)
      procs[idx] <- rep(TRUE, times = length(idx))
      cnt <- cnt+1
    }
  }

  #initialize shorter data for MCMC
  n_new <- cnt
  nd1 <- rep(0, times = n_new)
  nd2 <- rep(0, times = n_new)
  nn <- rep(0, times = n_new)
  nr <- rep(0, times = n_new)
  ns <- rep(0, times = n_new)

  #merge cohorts with same dose in the same study to compress data
  procs2 <- rep(F, length(dose1))
  cnt2 <- 0
  for(i in 1:length(dose1)){
    if(!procs2[i]){
      idx <- which(dose1[i]==dose1 &
                     dose2[i]==dose2 &
                     n.study[i] == n.study)
      procs2[idx] <- rep(TRUE, times = length(idx))
      cnt2 <- cnt2+1
      nd1[cnt2] <- dose1[i]
      nd2[cnt2] <- dose2[i]
      ns[cnt2] <- n.study[i]
      nn[cnt2] <- sum(n.pat[idx])
      nr[cnt2] <- sum(n.dlt[idx])
    }
  }

  nstd <- length(levels(factor(ns)))


  #deprecated call to sampling_jointBLRM, code was inlined below
  # samples <- sampling_jointBLRM(dose1=nd1,
  #                               dose2=nd2,
  #                               dose.ref1 = dose.ref1,
  #                               dose.ref2 = dose.ref2,
  #                               n.pat = nn,
  #                               n.dlt = nr,
  #                               n.study = ns,
  #                               prior.mu = prior.mu,
  #                               prior.tau = prior.tau,
  #                               MAP.prior = FALSE,
  #                               iter=iter,
  #                               warmup=warmup,
  #                               refresh=refresh,
  #                               chains=chains,
  #                               adapt_delta=adapt_delta,
  #                               max_treedepth=max_treedepth,
  #                               seed=seed
  # )


  #if file.name and working.path are supplied, check whether MCMC
  #result can be loaded
  if(!is.null(file.name) & !is.null(working.path)){
    #generate unique string based on data for MCMC call
    name_code <- paste0(file.name, "_",study.interest, "_",
                        n_new, "_", nstd, "_",
                        paste0(paste0(nd1, "-", nd2, "-", nn, "-", nr ,"-", ns), collapse ="_")#,
                        #".RData"
    )
    #create hash value from string
    #hashing is used to avoid overly long file names (may cause issues)
    name_hash <- sha256(name_code, key = "jointBLRM")
    #add file.name to allow cleaning up temporary or old results (e.g. with different priors)
    name_tmp <- paste0(file.name, "_tmp_", name_hash)

    #Check if hashed file name was already created
    if(paste0(name_tmp, ".RData") %in% list.files(working.path)){
      #path to file with MCMC results
      data_path <- file.path(working.path, name_tmp)
      #create file for locking
      lock_read <- lock(paste0(data_path,".lock"), exclusive = TRUE)
      #make sure that is it unlocked when exiting
      on.exit(expr = unlock(lock_read))
      #try to read MCMC results
      cont_load <- tryCatch(load(paste0(data_path,".RData")), error=function(e) e)
      if("res_all" %in% cont_load){
        #return if loaded
        return(res_all)
      }
    }

  }

  #Note: the following is only executed if the result could not be loaded previously
  #MCMC run with stan
  samples <- as.matrix(sampling(
    object = stanmodels$jointBLRM,
    data = list(
      n_obs = n_new,
      n_studies = nstd,
      dose_1 = as.array(nd1),
      dose_2 = as.array(nd2),
      dose_c = c(dose.ref1, dose.ref2),
      n = as.array(nn),
      r = as.array(nr),
      s = as.array(ns),
      saturating = ifelse(saturating, 1, 0),
      doMAP = 0,
      mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
                  prior.mu$mu_a2[1], prior.mu$mu_b2[1], prior.mu$mu_eta[1]),
      sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
                prior.mu$mu_a2[2], prior.mu$mu_b2[2], prior.mu$mu_eta[2]),
      mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
                   prior.tau$tau_a2[1], prior.tau$tau_b2[1], prior.tau$tau_eta[1]),
      sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
                 prior.tau$tau_a2[2], prior.tau$tau_b2[2], prior.tau$tau_eta[2])
    ),
    iter=iter, warmup=warmup, refresh=refresh, chains=chains, seed=seed,
    control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth),
    pars="log_ab"))

  # }else{
  #   #we have loaded samples from previous run;
  #   #simply use those
  #   samples <- samples_load
  # }



  nstudies <- length(study.interest)
  ndose <- length(dose1.interest)
  nsamp <- length(samples[,1])
  dose1_ <- dose1.interest/dose.ref1
  dose2_ <- dose2.interest/dose.ref2

  #Interval names not needed (inferred from context later)
  #ninters <- length(dosing.intervals)+1
  # name.inters <- rep(0, ninters)
  # for(i in 1:ninters){
  #   if(i==1){
  #     name.inters[i] <- paste0("P([0,", dosing.intervals[1], "))")
  #   }else if(i==ninters){
  #     name.inters[i] <- paste0("P([",dosing.intervals[ninters-1], ",1])")
  #   }else{
  #     name.inters[i] <- paste0("P([",dosing.intervals[i-1],",",dosing.intervals[i],"))")
  #   }
  # }

  study <- study.interest
  type <- type.interest
  #determine type of interest and adjust doses accordingly
  if(type=="mono1"){
    inds <- which(dose1_>0 & dose2_==0)
    dose1_c <- dose1_[inds]
    dose2_c <- dose2_[inds]
    ndosec <- length(dose1_c)
    dose1.interestc <- dose1.interest[inds]
    dose2.interestc <- dose2.interest[inds]
  }else if(type=="mono2"){
    inds <- which(dose2_>0 & dose1_==0)
    dose1_c <- dose1_[inds]
    dose2_c <- dose2_[inds]
    ndosec <- length(dose1_c)
    dose1.interestc <- dose1.interest[inds]
    dose2.interestc <- dose2.interest[inds]
  }else if(type=="combi"){
    inds <- which(dose1_>0 & dose2_>0)
    dose1_c <- dose1_[inds]
    dose2_c <- dose2_[inds]
    ndosec <- length(dose1_c)
    dose1.interestc <- dose1.interest[inds]
    dose2.interestc <- dose2.interest[inds]
  }else{
    dose1_c <- dose1_
    dose2_c <- dose2_
    ndosec <- ndose
    dose1.interestc <- dose1.interest
    dose2.interestc <- dose2.interest
  }

  logdose1_c <- log(dose1_c)
  logdose2_c <- log(dose2_c)

  #initialize vectors with samples from DLT rate
  toxs <- matrix(NA, nrow = ndosec, ncol = nsamp)
  toxs1 <- rep(0, nsamp)
  toxs2 <- rep(0, nsamp)
  toxs0 <- rep(0, nsamp)

  #get study-specific posterior samples
  la1 <- samples[, paste0("log_ab[", study, ",1]")]
  #lb1 <- samples[, paste0("log_ab[", study, ",2]")]
  la2 <- samples[, paste0("log_ab[", study, ",3]")]
  #lb2 <- samples[, paste0("log_ab[", study, ",4]")]
  eta <- samples[, paste0("log_ab[", study, ",5]")]
  #explb1 <- exp(lb1)
  #explb2 <- exp(lb2)
  explb1 <- exp(samples[, paste0("log_ab[", study, ",2]")])
  explb2 <- exp(samples[, paste0("log_ab[", study, ",4]")])

  #for each dose: determine posterior DLT rate samples.
  for(d in 1:ndosec){
    if(dose2_c[d]==0 & dose1_c[d]>0){
      #treatment mono1
      #toxs[d, ] <- inv_logit(la1 + explb1*logdose1_c[d])
      toxs[d, ] <- inv_logit(samples[, paste0("log_ab[", study, ",1]")] + explb1*logdose1_c[d])
    }else if(dose1_c[d]==0 & dose2_c[d]>0){
      #treatment mono2
      #toxs[d, ] <- inv_logit(la2 + explb2*logdose2_c[d])
      toxs[d, ] <- inv_logit(samples[, paste0("log_ab[", study, ",3]")] + explb2*logdose2_c[d])
    }else if(dose2_c[d]>0 & dose1_c[d]>0){
      #treatment combination
      #toxs1 <- inv_logit(la1 + explb1*logdose1_c[d])
      #toxs2 <- inv_logit(la2 + explb2*logdose2_c[d])
      toxs1 <- inv_logit(samples[, paste0("log_ab[", study, ",1]")] + explb1*logdose1_c[d])
      toxs2 <- inv_logit(samples[, paste0("log_ab[", study, ",3]")] + explb2*logdose2_c[d])
      toxs0 <- toxs1 + toxs2 - toxs1*toxs2
      if(saturating){
        #toxs[d, ] <- inv_logit(logit(toxs0) + 2*eta*dose1_c[d]*dose2_c[d]/(1+dose1_c[d]*dose2_c[d]))
        toxs[d, ] <- inv_logit(logit(toxs0) + 2*samples[, paste0("log_ab[", study, ",5]")]*dose1_c[d]*dose2_c[d]/(1+dose1_c[d]*dose2_c[d]))
      }else{
        #toxs[d, ] <- inv_logit(logit(toxs0) + eta*dose1_c[d]*dose2_c[d])
        toxs[d, ] <- inv_logit(logit(toxs0) + samples[, paste0("log_ab[", study, ",5]")]*dose1_c[d]*dose2_c[d])
      }
    }else{
      #no treatment, put in 0 as dummy
      toxs[d, ] <- rep(0, times = nsamp)
    }
  }

  #number of dosing intervals and return structure
  n.entries <- length(dosing.intervals)+1
  res_curr <- matrix(NA, nrow = ndosec, ncol = n.entries)
  rownames(res_curr) <- paste0(dose1.interestc, "+", dose2.interestc)

  #determine interval probabilities
  for(d in 1:ndosec){
    for(inter in 1:n.entries){
      if(inter==1){
        res_curr[d, 1] <- mean(ifelse(toxs[d,]<dosing.intervals[1], 1, 0))
      }else if(inter==n.entries){
        res_curr[d, n.entries] <- mean(ifelse(toxs[d,]>=dosing.intervals[n.entries-1], 1, 0))
      }else{
        res_curr[d, inter] <- mean(ifelse(toxs[d,]>=dosing.intervals[inter-1] &
                                            toxs[d,]<dosing.intervals[inter],
                                          1, 0))
      }
    }
  }

  #only when dynamic loss is used
  if(d.loss){
    #compute samples from reference toxicities for study
    #(only required for dynamic loss escalation)
    #Note: reference dose not guaranteed to be one of the available doses in a
    #trial, so needs to be computed separately
    ref_toxs <- matrix(NA, nrow = 3, ncol = nsamp)
    rtoxs0 <- rep(0, nsamp)
    #Note: model for DLT rate simplifies for reference dose!
    ref_toxs[1, ] <- inv_logit(la1)
    ref_toxs[2, ] <- inv_logit(la2)
    logit_rtoxs0 <- logit(ref_toxs[1,]+ref_toxs[2,]-ref_toxs[1,]*ref_toxs[2,])
    ref_toxs[3, ] <- inv_logit(logit_rtoxs0 + eta)

    #determine interval probabilities
    interprobs_ref <- matrix(NA, nrow=3, ncol=n.entries)
    #colnames(interprobs_ref) <- name.inters
    rownames(interprobs_ref) <- c("Ref.Mono1", "Ref.Mono2", "Ref.Combi")
    for(d in 1:3){
      for(inter in 1:n.entries){
        if(inter==1){
          interprobs_ref[d, 1] <- mean(ifelse(ref_toxs[d,]<dosing.intervals[1], 1, 0))
        }else if(inter==n.entries){
          interprobs_ref[d, n.entries] <- mean(ifelse(ref_toxs[d,]>=dosing.intervals[n.entries-1], 1, 0))
        }else{
          interprobs_ref[d, inter] <- mean(ifelse(ref_toxs[d,]>=dosing.intervals[inter-1] &
                                                    ref_toxs[d,]<dosing.intervals[inter],
                                                  1, 0))
        }
      }
    }

    res_all <- list(
      "result"=res_curr,
      "refprobs"=interprobs_ref
    )

  }else{
    #Not dynamic loss: only return interval probabilities
    res_all <- res_curr
  }

  #save MCMC result for later use, and return
  if(!is.null(file.name) & !is.null(working.path)){
    #name_tmp is the previously created hash name for the result
    data_path <- file.path(working.path, name_tmp)
    #wait until potential other workers are finished processing this data
    lock_read <- lock(paste0(data_path,".lock"), exclusive = TRUE)
    #unlock for access.
    on.exit(unlock(lock_read))
    #save return value
    #browser()
    tryCatch(save(res_all, file = paste0(data_path,".RData")), error = function(e) e)
  }

  return(res_all)

  # #alternative case: simply return loaded result
  # }else{
  #
  # }

}


#' @keywords internal
post_tox_covariate_jointBLRM_sim <- function(
  study.interest,
  type.interest,
  cov.interest, ##NEW
  dose1.interest,
  dose2.interest,
  BLRM,
  file.name,
  working.path,
  dose1,
  dose2,
  dose.ref1,
  dose.ref2,
  n.pat,
  n.dlt,
  n.study,
  cov, ##NEW
  d.loss = FALSE,
  dosing.intervals = c(0.16, 0.33),
  prior.mu = list(mu_a1 =  c(-1.386294, 2),
                  mu_b1 =  c(0,         1),
                  mu_a2 =  c(-1.386294, 2),
                  mu_b2 =  c(0,         1),
                  mu_eta = c(0,         1.121),

                  mu_c1 = c(0, 1),   #NEW
                  mu_c2 = c(0, 1)
  ),
  prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
                   tau_b1 =  c(log(0.125), log(2)/1.96),
                   tau_a2 =  c(log(0.25),  log(2)/1.96),
                   tau_b2 =  c(log(0.125), log(2)/1.96),
                   tau_eta = c(log(0.125), log(2)/1.96),

                   tau_c1 = c(log(0.125), log(2)/1.96),   #NEW
                   tau_c2 = c(log(0.125), log(2)/1.96)    #NEW
  ),

  two_sided1 = TRUE,                             #'NEW'
  two_sided2 = TRUE,                             #'NEW'

  saturating = FALSE,
  iter = 10000,
  adapt_delta = 0.95,
  max_treedepth = 15,
  warmup = floor(iter/2),
  refresh = floor(iter/10),
  chains = 4,
  #n_tries = 10,
  seed=sample.int(.Machine$integer.max, 1)
)
{

  #determine how many observations remain after collapsing data
  procs <- rep(F, length(dose1))
  cnt <- 0
  for(i in 1:length(dose1)){
    if(!procs[i]){
      idx <- which(dose1[i]==dose1 &
                     dose2[i]==dose2 &
                     n.study[i] == n.study &
                     cov[i] == cov)
      procs[idx] <- rep(TRUE, times = length(idx))
      cnt <- cnt+1
    }
  }

  #initialize shorter data for MCMC
  n_new <- cnt
  nd1 <- rep(0, times = n_new)
  nd2 <- rep(0, times = n_new)
  nn <- rep(0, times = n_new)
  nr <- rep(0, times = n_new)
  ns <- rep(0, times = n_new)
  nc <- rep(0, times = n_new)

  #merge cohorts with same dose in the same study to compress data
  procs2 <- rep(F, length(dose1))
  cnt2 <- 0
  for(i in 1:length(dose1)){
    if(!procs2[i]){
      idx <- which(dose1[i]==dose1 &
                     dose2[i]==dose2 &
                     n.study[i] == n.study &
                     cov[i] == cov)
      procs2[idx] <- rep(TRUE, times = length(idx))
      cnt2 <- cnt2+1
      nd1[cnt2] <- dose1[i]
      nd2[cnt2] <- dose2[i]
      ns[cnt2] <- n.study[i]
      nn[cnt2] <- sum(n.pat[idx])
      nr[cnt2] <- sum(n.dlt[idx])
      nc[cnt2] <- cov[i]
    }
  }

  nstd <- length(levels(factor(ns)))


  #deprecated call to sampling_jointBLRM, code was inlined below
  # samples <- sampling_jointBLRM(dose1=nd1,
  #                               dose2=nd2,
  #                               dose.ref1 = dose.ref1,
  #                               dose.ref2 = dose.ref2,
  #                               n.pat = nn,
  #                               n.dlt = nr,
  #                               n.study = ns,
  #                               prior.mu = prior.mu,
  #                               prior.tau = prior.tau,
  #                               MAP.prior = FALSE,
  #                               iter=iter,
  #                               warmup=warmup,
  #                               refresh=refresh,
  #                               chains=chains,
  #                               adapt_delta=adapt_delta,
  #                               max_treedepth=max_treedepth,
  #                               seed=seed
  # )


  #if file.name and working.path are supplied, check whether MCMC
  #result can be loaded
  if(!is.null(file.name) & !is.null(working.path)){
    #generate unique string based on data for MCMC call
    name_code <- paste0(file.name, "_",study.interest, "_",
                        n_new, "_", nstd, "_",
                        paste0(paste0(nd1, "-", nd2, "-", nn, "-", nr ,"-", ns, "-", nc), collapse ="_")#,
                        #".RData"
    )
    #create hash value from string
    #hashing is used to avoid overly long file names (may cause issues)
    name_hash <- sha256(name_code, key = "jointBLRM")
    #add file.name to allow cleaning up temporary or old results (e.g. with different priors)
    name_tmp <- paste0(file.name, "_tmp_", name_hash)

    #Check if hashed file name was already created
    if(paste0(name_tmp, ".RData") %in% list.files(working.path)){
      #path to file with MCMC results
      data_path <- file.path(working.path, name_tmp)
      #create file for locking
      lock_read <- lock(paste0(data_path,".lock"), exclusive = TRUE)
      #make sure that is it unlocked when exiting
      on.exit(expr = unlock(lock_read))
      #try to read MCMC results
      cont_load <- tryCatch(load(paste0(data_path,".RData")), error=function(e) e)
      if("res_all" %in% cont_load){
        #return if loaded
        return(res_all)
      }
    }

  }

  #Note: the following is only executed if the result could not be loaded previously
  #MCMC run with stan
  samples <- as.matrix(sampling(
    object = stanmodels$jointBLRMcov,
    data = list(
      n_obs = n_new,
      n_studies = nstd,
      dose_1 = as.array(nd1),
      dose_2 = as.array(nd2),
      dose_c = c(dose.ref1, dose.ref2),
      n = as.array(nn),
      r = as.array(nr),
      s = as.array(ns),
      c = as.array(nc),
      saturating = ifelse(saturating, 1, 0),
      doMAP = 0,
      twoside1 = two_sided1,
      twoside2 = two_sided2,
      mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
                  prior.mu$mu_a2[1], prior.mu$mu_b2[1], prior.mu$mu_eta[1],
                  prior.mu$mu_c1[1], prior.mu$mu_c2[1]),
      sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
                prior.mu$mu_a2[2], prior.mu$mu_b2[2], prior.mu$mu_eta[2],
                prior.mu$mu_c1[2], prior.mu$mu_c2[2]),
      mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
                   prior.tau$tau_a2[1], prior.tau$tau_b2[1], prior.tau$tau_eta[1],
                   prior.tau$tau_c1[1], prior.tau$tau_c2[1]),
      sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
                 prior.tau$tau_a2[2], prior.tau$tau_b2[2], prior.tau$tau_eta[2],
                 prior.tau$tau_c1[2], prior.tau$tau_c2[2])
    ),
    iter=iter, warmup=warmup, refresh=refresh, chains=chains, seed=seed,
    control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth),
    pars="log_ab"))

  # }else{
  #   #we have loaded samples from previous run;
  #   #simply use those
  #   samples <- samples_load
  # }



  nstudies <- length(study.interest)
  ndose <- length(dose1.interest)
  nsamp <- length(samples[,1])
  dose1_ <- dose1.interest/dose.ref1
  dose2_ <- dose2.interest/dose.ref2

  #Interval names not needed (inferred from context later)
  #ninters <- length(dosing.intervals)+1
  # name.inters <- rep(0, ninters)
  # for(i in 1:ninters){
  #   if(i==1){
  #     name.inters[i] <- paste0("P([0,", dosing.intervals[1], "))")
  #   }else if(i==ninters){
  #     name.inters[i] <- paste0("P([",dosing.intervals[ninters-1], ",1])")
  #   }else{
  #     name.inters[i] <- paste0("P([",dosing.intervals[i-1],",",dosing.intervals[i],"))")
  #   }
  # }

  study <- study.interest
  type <- type.interest
  #determine type of interest and adjust doses accordingly
  if(type=="mono1"){
    inds <- which(dose1_>0 & dose2_==0)
    dose1_c <- dose1_[inds]
    dose2_c <- dose2_[inds]
    ndosec <- length(dose1_c)
    dose1.interestc <- dose1.interest[inds]
    dose2.interestc <- dose2.interest[inds]
  }else if(type=="mono2"){
    inds <- which(dose2_>0 & dose1_==0)
    dose1_c <- dose1_[inds]
    dose2_c <- dose2_[inds]
    ndosec <- length(dose1_c)
    dose1.interestc <- dose1.interest[inds]
    dose2.interestc <- dose2.interest[inds]
  }else if(type=="combi"){
    inds <- which(dose1_>0 & dose2_>0)
    dose1_c <- dose1_[inds]
    dose2_c <- dose2_[inds]
    ndosec <- length(dose1_c)
    dose1.interestc <- dose1.interest[inds]
    dose2.interestc <- dose2.interest[inds]
  }else{
    dose1_c <- dose1_
    dose2_c <- dose2_
    ndosec <- ndose
    dose1.interestc <- dose1.interest
    dose2.interestc <- dose2.interest
  }


  cov_ <- rep(0, times = ndosec)
  if(cov.interest == 1){
    cov_ <- rep(1, times = ndosec)
  }

  logdose1_c <- log(dose1_c)
  logdose2_c <- log(dose2_c)

  #initialize vectors with samples from DLT rate
  toxs <- matrix(NA, nrow = ndosec, ncol = nsamp)
  toxs1 <- rep(0, nsamp)
  toxs2 <- rep(0, nsamp)
  toxs0 <- rep(0, nsamp)

  #get study-specific posterior samples
  la1 <- samples[, paste0("log_ab[", study, ",1]")]
  #lb1 <- samples[, paste0("log_ab[", study, ",2]")]
  la2 <- samples[, paste0("log_ab[", study, ",3]")]
  #lb2 <- samples[, paste0("log_ab[", study, ",4]")]
  eta <- samples[, paste0("log_ab[", study, ",5]")]
  #explb1 <- exp(lb1)
  #explb2 <- exp(lb2)
  explb1 <- exp(samples[, paste0("log_ab[", study, ",2]")])
  explb2 <- exp(samples[, paste0("log_ab[", study, ",4]")])


  if(two_sided1){
    c1 <- samples[, paste0("log_ab[", study, ",6]")]
  }else{
    c1 <- exp(samples[, paste0("log_ab[", study, ",6]")])
  }

  if(two_sided2){
    c2 <- samples[, paste0("log_ab[", study, ",7]")]
  }else{
    c2 <- exp(samples[, paste0("log_ab[", study, ",7]")])
  }
  #for each dose: determine posterior DLT rate samples.
  for(d in 1:ndosec){
    if(dose2_c[d]==0 & dose1_c[d]>0){
      #treatment mono1
      #toxs[d, ] <- inv_logit(la1 + explb1*logdose1_c[d])
      toxs[d, ] <- inv_logit(samples[, paste0("log_ab[", study, ",1]")] + explb1*logdose1_c[d]
                             + cov_[d] * c1)
    }else if(dose1_c[d]==0 & dose2_c[d]>0){
      #treatment mono2
      #toxs[d, ] <- inv_logit(la2 + explb2*logdose2_c[d])
      toxs[d, ] <- inv_logit(samples[, paste0("log_ab[", study, ",3]")] + explb2*logdose2_c[d]
                             + cov_[d] * c2)
    }else if(dose2_c[d]>0 & dose1_c[d]>0){
      #treatment combination
      #toxs1 <- inv_logit(la1 + explb1*logdose1_c[d])
      #toxs2 <- inv_logit(la2 + explb2*logdose2_c[d])
      toxs1 <- inv_logit(samples[, paste0("log_ab[", study, ",1]")] + explb1*logdose1_c[d]
                         + cov_[d] * c1)
      toxs2 <- inv_logit(samples[, paste0("log_ab[", study, ",3]")] + explb2*logdose2_c[d]
                         + cov_[d] * c2)
      toxs0 <- toxs1 + toxs2 - toxs1*toxs2
      if(saturating){
        #toxs[d, ] <- inv_logit(logit(toxs0) + 2*eta*dose1_c[d]*dose2_c[d]/(1+dose1_c[d]*dose2_c[d]))
        toxs[d, ] <- inv_logit(logit(toxs0) + 2*samples[, paste0("log_ab[", study, ",5]")]*dose1_c[d]*dose2_c[d]/(1+dose1_c[d]*dose2_c[d]))
      }else{
        #toxs[d, ] <- inv_logit(logit(toxs0) + eta*dose1_c[d]*dose2_c[d])
        toxs[d, ] <- inv_logit(logit(toxs0) + samples[, paste0("log_ab[", study, ",5]")]*dose1_c[d]*dose2_c[d])
      }
    }else{
      #no treatment, put in 0 as dummy
      toxs[d, ] <- rep(0, times = nsamp)
    }
  }

  #number of dosing intervals and return structure
  n.entries <- length(dosing.intervals)+1
  res_curr <- matrix(NA, nrow = ndosec, ncol = n.entries)
  rownames(res_curr) <- paste0(dose1.interestc, "+", dose2.interestc)

  #determine interval probabilities
  for(d in 1:ndosec){
    for(inter in 1:n.entries){
      if(inter==1){
        res_curr[d, 1] <- mean(ifelse(toxs[d,]<dosing.intervals[1], 1, 0))
      }else if(inter==n.entries){
        res_curr[d, n.entries] <- mean(ifelse(toxs[d,]>=dosing.intervals[n.entries-1], 1, 0))
      }else{
        res_curr[d, inter] <- mean(ifelse(toxs[d,]>=dosing.intervals[inter-1] &
                                            toxs[d,]<dosing.intervals[inter],
                                          1, 0))
      }
    }
  }

  #only when dynamic loss is used
  if(d.loss){
    #compute samples from reference toxicities for study
    #(only required for dynamic loss escalation)
    #Note: reference dose not guaranteed to be one of the available doses in a
    #trial, so needs to be computed separately
    ref_toxs <- matrix(NA, nrow = 3, ncol = nsamp)
    rtoxs0 <- rep(0, nsamp)
    #Note: model for DLT rate simplifies for reference dose!
    ref_toxs[1, ] <- inv_logit(la1 + cov.interest*c1)
    ref_toxs[2, ] <- inv_logit(la2 + cov.interest*c2)
    logit_rtoxs0 <- logit(ref_toxs[1,]+ref_toxs[2,]-ref_toxs[1,]*ref_toxs[2,])
    ref_toxs[3, ] <- inv_logit(logit_rtoxs0 + eta)

    #determine interval probabilities
    interprobs_ref <- matrix(NA, nrow=3, ncol=n.entries)
    #colnames(interprobs_ref) <- name.inters
    rownames(interprobs_ref) <- c("Ref.Mono1", "Ref.Mono2", "Ref.Combi")
    for(d in 1:3){
      for(inter in 1:n.entries){
        if(inter==1){
          interprobs_ref[d, 1] <- mean(ifelse(ref_toxs[d,]<dosing.intervals[1], 1, 0))
        }else if(inter==n.entries){
          interprobs_ref[d, n.entries] <- mean(ifelse(ref_toxs[d,]>=dosing.intervals[n.entries-1], 1, 0))
        }else{
          interprobs_ref[d, inter] <- mean(ifelse(ref_toxs[d,]>=dosing.intervals[inter-1] &
                                                    ref_toxs[d,]<dosing.intervals[inter],
                                                  1, 0))
        }
      }
    }

    res_all <- list(
      "result"=res_curr,
      "refprobs"=interprobs_ref
    )

  }else{
    #Not dynamic loss: only return interval probabilities
    res_all <- res_curr
  }

  #save MCMC result for later use, and return
  if(!is.null(file.name) & !is.null(working.path)){
    #name_tmp is the previously created hash name for the result
    data_path <- file.path(working.path, name_tmp)
    #wait until potential other workers are finished processing this data
    lock_read <- lock(paste0(data_path,".lock"), exclusive = TRUE)
    #unlock for access.
    on.exit(unlock(lock_read))
    #save return value
    #browser()
    tryCatch(save(res_all, file = paste0(data_path,".RData")), error = function(e) e)
  }

  return(res_all)

  # #alternative case: simply return loaded result
  # }else{
  #
  # }

}


#-------------------------------------------------------------------------------
#Old Version of post_tox_jointBLRM_sim
#-------------------------------------------------------------------------------

#Function used by sim_jointBLRM to determine interval probabilities.
#Faster and less structured version of "post_tox_joitBLRM".
#Changes:
# -Only works for single trials of interest.
# -Inlines the function for sampling.
# -compresses data by merging observations at the same dose in the same study.
# -output has less structure as it only needs to work with sim_jointBLRM.


# @keywords internal
# post_tox_jointBLRM_sim <- function(
#   study.interest,
#   type.interest,
#   dose1.interest,
#   dose2.interest,
#   dose1,
#   dose2,
#   dose.ref1,
#   dose.ref2,
#   n.pat,
#   n.dlt,
#   n.study,
#   d.loss = FALSE,
#   dosing.intervals = c(0.16, 0.33),
#   prior.mu = list(mu_a1 =  c(-1.386294, 2),
#                   mu_b1 =  c(0,         1),
#                   mu_a2 =  c(-1.386294, 2),
#                   mu_b2 =  c(0,         1),
#                   mu_eta = c(0,         1.121)
#   ),
#   prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
#                    tau_b1 =  c(log(0.125), log(2)/1.96),
#                    tau_a2 =  c(log(0.25),  log(2)/1.96),
#                    tau_b2 =  c(log(0.125), log(2)/1.96),
#                    tau_eta = c(log(0.125), log(2)/1.96)
#   ),
#   saturating = FALSE,
#   iter = 10000,
#   adapt_delta = 0.95,
#   max_treedepth = 15,
#   warmup = floor(iter/2),
#   refresh = floor(iter/10),
#   chains = 4,
#   seed=sample.int(.Machine$integer.max, 1)
# )
# {
#
#   #determine how many observations remain after collapsing data
#   procs <- rep(F, length(dose1))
#   cnt <- 0
#   for(i in 1:length(dose1)){
#     if(!procs[i]){
#       idx <- which(dose1[i]==dose1 &
#                      dose2[i]==dose2 &
#                      n.study[i] == n.study)
#       procs[idx] <- rep(TRUE, times = length(idx))
#       cnt <- cnt+1
#     }
#   }
#
#   #initialize shorter data for MCMC
#   n_new <- cnt
#   nd1 <- rep(0, times = n_new)
#   nd2 <- rep(0, times = n_new)
#   nn <- rep(0, times = n_new)
#   nr <- rep(0, times = n_new)
#   ns <- rep(0, times = n_new)
#
#   #merge cohorts with same dose in the same study to compress data
#   procs2 <- rep(F, length(dose1))
#   cnt2 <- 0
#   for(i in 1:length(dose1)){
#     if(!procs2[i]){
#       idx <- which(dose1[i]==dose1 &
#                      dose2[i]==dose2 &
#                      n.study[i] == n.study)
#       procs2[idx] <- rep(TRUE, times = length(idx))
#       cnt2 <- cnt2+1
#       nd1[cnt2] <- dose1[i]
#       nd2[cnt2] <- dose2[i]
#       ns[cnt2] <- n.study[i]
#       nn[cnt2] <- sum(n.pat[idx])
#       nr[cnt2] <- sum(n.dlt[idx])
#     }
#   }
#
#   #deprecated call to sampling_jointBLRM, code was inlined below
#   # samples <- sampling_jointBLRM(dose1=nd1,
#   #                               dose2=nd2,
#   #                               dose.ref1 = dose.ref1,
#   #                               dose.ref2 = dose.ref2,
#   #                               n.pat = nn,
#   #                               n.dlt = nr,
#   #                               n.study = ns,
#   #                               prior.mu = prior.mu,
#   #                               prior.tau = prior.tau,
#   #                               MAP.prior = FALSE,
#   #                               iter=iter,
#   #                               warmup=warmup,
#   #                               refresh=refresh,
#   #                               chains=chains,
#   #                               adapt_delta=adapt_delta,
#   #                               max_treedepth=max_treedepth,
#   #                               seed=seed
#   # )
#
#   nstd <- length(levels(factor(ns)))
#
#   #MCMC run with stan
#     samples <- as.matrix(sampling(
#       object = stanmodels$jointBLRM,
#       data = list(
#         n_obs = n_new,
#         n_studies = nstd,
#         dose_1 = as.array(nd1),
#         dose_2 = as.array(nd2),
#         dose_c = c(dose.ref1, dose.ref2),
#         n = as.array(nn),
#         r = as.array(nr),
#         s = as.array(ns),
#         doMAP = 0,
#         saturating = ifelse(saturating, 1, 0),
#         mean_mu = c(prior.mu$mu_a1[1], prior.mu$mu_b1[1],
#                     prior.mu$mu_a2[1], prior.mu$mu_b2[1], prior.mu$mu_eta[1]),
#         sd_mu = c(prior.mu$mu_a1[2], prior.mu$mu_b1[2],
#                   prior.mu$mu_a2[2], prior.mu$mu_b2[2], prior.mu$mu_eta[2]),
#         mean_tau = c(prior.tau$tau_a1[1], prior.tau$tau_b1[1],
#                      prior.tau$tau_a2[1], prior.tau$tau_b2[1], prior.tau$tau_eta[1]),
#         sd_tau = c(prior.tau$tau_a1[2], prior.tau$tau_b1[2],
#                    prior.tau$tau_a2[2], prior.tau$tau_b2[2], prior.tau$tau_eta[2])
#       ),
#       iter=iter, warmup=warmup, refresh=refresh, chains=chains, seed=seed,
#       control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth),
#       pars="log_ab"))
#
#
#   nstudies <- length(study.interest)
#   ndose <- length(dose1.interest)
#   nsamp <- length(samples[,1])
#   dose1_ <- dose1.interest/dose.ref1
#   dose2_ <- dose2.interest/dose.ref2
#
#   #Interval names not needed (inferred from context later)
#   #ninters <- length(dosing.intervals)+1
#   # name.inters <- rep(0, ninters)
#   # for(i in 1:ninters){
#   #   if(i==1){
#   #     name.inters[i] <- paste0("P([0,", dosing.intervals[1], "))")
#   #   }else if(i==ninters){
#   #     name.inters[i] <- paste0("P([",dosing.intervals[ninters-1], ",1])")
#   #   }else{
#   #     name.inters[i] <- paste0("P([",dosing.intervals[i-1],",",dosing.intervals[i],"))")
#   #   }
#   # }
#
#   study <- study.interest
#   type <- type.interest
#   #determine type of interest and adjust doses accordingly
#   if(type=="mono1"){
#     inds <- which(dose1_>0 & dose2_==0)
#     dose1_c <- dose1_[inds]
#     dose2_c <- dose2_[inds]
#     ndosec <- length(dose1_c)
#     dose1.interestc <- dose1.interest[inds]
#     dose2.interestc <- dose2.interest[inds]
#   }else if(type=="mono2"){
#     inds <- which(dose2_>0 & dose1_==0)
#     dose1_c <- dose1_[inds]
#     dose2_c <- dose2_[inds]
#     ndosec <- length(dose1_c)
#     dose1.interestc <- dose1.interest[inds]
#     dose2.interestc <- dose2.interest[inds]
#   }else if(type=="combi"){
#     inds <- which(dose1_>0 & dose2_>0)
#     dose1_c <- dose1_[inds]
#     dose2_c <- dose2_[inds]
#     ndosec <- length(dose1_c)
#     dose1.interestc <- dose1.interest[inds]
#     dose2.interestc <- dose2.interest[inds]
#   }else{
#     dose1_c <- dose1_
#     dose2_c <- dose2_
#     ndosec <- ndose
#     dose1.interestc <- dose1.interest
#     dose2.interestc <- dose2.interest
#   }
#
#   logdose1_c <- log(dose1_c)
#   logdose2_c <- log(dose2_c)
#
#   #initialize vectors with samples from DLT rate
#   toxs <- matrix(NA, nrow = ndosec, ncol = nsamp)
#   toxs1 <- rep(0, nsamp)
#   toxs2 <- rep(0, nsamp)
#   toxs0 <- rep(0, nsamp)
#
#   #get study-specific posterior samples
#   la1 <- samples[, paste0("log_ab[", study, ",1]")]
#   lb1 <- samples[, paste0("log_ab[", study, ",2]")]
#   la2 <- samples[, paste0("log_ab[", study, ",3]")]
#   lb2 <- samples[, paste0("log_ab[", study, ",4]")]
#   eta <- samples[, paste0("log_ab[", study, ",5]")]
#   explb1 <- exp(lb1)
#   explb2 <- exp(lb2)
#
#   #for each dose: determine posterior DLT rate samples.
#   if(saturating){
#     for(d in 1:ndosec){
#       if(dose2_c[d]==0 & dose1_c[d]>0){
#         #treatment mono1
#         toxs[d, ] <- inv_logit(la1 + explb1*logdose1_c[d])
#       }else if(dose1_c[d]==0 & dose2_c[d]>0){
#         #treatment mono2
#         toxs[d, ] <- inv_logit(la2 + explb2*logdose2_c[d])
#       }else if(dose2_c[d]>0 & dose1_c[d]>0){
#         #treatment combination
#         toxs1 <- inv_logit(la1 + explb1*logdose1_c[d])
#         toxs2 <- inv_logit(la2 + explb2*logdose2_c[d])
#         toxs0 <- toxs1 + toxs2 - toxs1*toxs2
#         toxs[d, ] <- inv_logit(logit(toxs0) + 2*eta*dose1_c[d]*dose2_c[d]/(1+dose1_c[d]*dose2_c[d]))
#
#       }else{
#         #no treatment, put in 0 as dummy
#         toxs[d, ] <- rep(0, times = nsamp)
#       }
#     }
#   }else{
#     for(d in 1:ndosec){
#       if(dose2_c[d]==0 & dose1_c[d]>0){
#         #treatment mono1
#         toxs[d, ] <- inv_logit(la1 + explb1*logdose1_c[d])
#       }else if(dose1_c[d]==0 & dose2_c[d]>0){
#         #treatment mono2
#         toxs[d, ] <- inv_logit(la2 + explb2*logdose2_c[d])
#       }else if(dose2_c[d]>0 & dose1_c[d]>0){
#         #treatment combination
#         toxs1 <- inv_logit(la1 + explb1*logdose1_c[d])
#         toxs2 <- inv_logit(la2 + explb2*logdose2_c[d])
#         toxs0 <- toxs1 + toxs2 - toxs1*toxs2
#         toxs[d, ] <- inv_logit(logit(toxs0) + eta*dose1_c[d]*dose2_c[d])
#
#       }else{
#         #no treatment, put in 0 as dummy
#         toxs[d, ] <- rep(0, times = nsamp)
#       }
#     }
#   }
#   #number of dosing intervals and return structure
#   n.entries <- length(dosing.intervals)+1
#   res_curr <- matrix(NA, nrow = ndosec, ncol = n.entries)
#   rownames(res_curr) <- paste0(dose1.interestc, "+", dose2.interestc)
#
#   #determine interval probabilities
#   for(d in 1:ndosec){
#     for(inter in 1:n.entries){
#       if(inter==1){
#         res_curr[d, 1] <- mean(ifelse(toxs[d,]<dosing.intervals[1], 1, 0))
#       }else if(inter==n.entries){
#         res_curr[d, n.entries] <- mean(ifelse(toxs[d,]>=dosing.intervals[n.entries-1], 1, 0))
#       }else{
#         res_curr[d, inter] <- mean(ifelse(toxs[d,]>=dosing.intervals[inter-1] &
#                                             toxs[d,]<dosing.intervals[inter],
#                                           1, 0))
#       }
#     }
#   }
#
#   #only when dynamic loss is used
#   if(d.loss){
#     #compute samples from reference toxicities for study
#     #(only required for dynamic loss escalation)
#     #Note: reference dose not guaranteed to be one of the available doses in a
#     #trial, so needs to be computed separately
#     ref_toxs <- matrix(NA, nrow = 3, ncol = nsamp)
#     rtoxs0 <- rep(0, nsamp)
#     #Note: model for DLT rate simplifies for reference dose!
#     ref_toxs[1, ] <- inv_logit(la1)
#     ref_toxs[2, ] <- inv_logit(la2)
#     logit_rtoxs0 <- logit(ref_toxs[1,]+ref_toxs[2,]-ref_toxs[1,]*ref_toxs[2,])
#     ref_toxs[3, ] <- inv_logit(logit_rtoxs0 + eta)
#
#     #determine interval probabilities
#     interprobs_ref <- matrix(NA, nrow=3, ncol=n.entries)
#     #colnames(interprobs_ref) <- name.inters
#     rownames(interprobs_ref) <- c("Ref.Mono1", "Ref.Mono2", "Ref.Combi")
#     for(d in 1:3){
#       for(inter in 1:n.entries){
#         if(inter==1){
#           interprobs_ref[d, 1] <- mean(ifelse(ref_toxs[d,]<dosing.intervals[1], 1, 0))
#         }else if(inter==n.entries){
#           interprobs_ref[d, n.entries] <- mean(ifelse(ref_toxs[d,]>=dosing.intervals[n.entries-1], 1, 0))
#         }else{
#           interprobs_ref[d, inter] <- mean(ifelse(ref_toxs[d,]>=dosing.intervals[inter-1] &
#                                                     ref_toxs[d,]<dosing.intervals[inter],
#                                                   1, 0))
#         }
#       }
#     }
#     return(list(
#       "result"=res_curr,
#       "refprobs"=interprobs_ref
#     ))
#
#   }else{
#     #Not dynamic loss: only return interval probabilities
#     return(res_curr)
#   }
#
#
# }


# #-------------------------------------------------------------------------------
# #Deprecated: old version
# #-------------------------------------------------------------------------------
# # @keywords internal
#
# post_tox_jointBLRM_sim_old <- function(
#   study.interest,
#   type.interest,
#   dose1.interest,
#   dose2.interest,
#   dose1,
#   dose2,
#   dose.ref1,
#   dose.ref2,
#   n.pat,
#   n.dlt,
#   n.study,
#   d.loss = FALSE,
#   dosing.intervals = c(0.16, 0.33),
#   prior.mu = list(mu_a1 =  c(-1.386294, 2),
#                   mu_b1 =  c(0,         1),
#                   mu_a2 =  c(-1.386294, 2),
#                   mu_b2 =  c(0,         1),
#                   mu_eta = c(0,         1.121)
#   ),
#   prior.tau = list(tau_a1 =  c(log(0.25),  log(2)/1.96),
#                    tau_b1 =  c(log(0.125), log(2)/1.96),
#                    tau_a2 =  c(log(0.25),  log(2)/1.96),
#                    tau_b2 =  c(log(0.125), log(2)/1.96),
#                    tau_eta = c(log(0.125), log(2)/1.96)
#   ),
#   iter = 10000,
#   adapt_delta = 0.95,
#   max_treedepth = 15,
#   warmup = floor(iter/2),
#   refresh = floor(iter/10),
#   chains = 4,
#   seed=sample.int(.Machine$integer.max, 1)
# )
# {
#
#
#   samples <- sampling_jointBLRM(dose1=dose1,
#                                 dose2=dose2,
#                                 dose.ref1 = dose.ref1,
#                                 dose.ref2 = dose.ref2,
#                                 n.pat = n.pat,
#                                 n.dlt = n.dlt,
#                                 n.study = n.study,
#                                 prior.mu = prior.mu,
#                                 prior.tau = prior.tau,
#                                 MAP.prior = FALSE,
#                                 iter=iter,
#                                 warmup=warmup,
#                                 refresh=refresh,
#                                 chains=chains,
#                                 adapt_delta=adapt_delta,
#                                 max_treedepth=max_treedepth,
#                                 seed=seed
#   )
#
#   nstudies <- length(study.interest)
#   ndose <- length(dose1.interest)
#   nsamp <- length(samples[,1])
#   dose1_ <- dose1.interest/dose.ref1
#   dose2_ <- dose2.interest/dose.ref2
#
#   ninters <- length(dosing.intervals)+1
#   name.inters <- rep(0, ninters)
#   for(i in 1:ninters){
#     if(i==1){
#       name.inters[i] <- paste0("P([0,", dosing.intervals[1], "))")
#     }else if(i==ninters){
#       name.inters[i] <- paste0("P([",dosing.intervals[ninters-1], ",1])")
#     }else{
#       name.inters[i] <- paste0("P([",dosing.intervals[i-1],",",dosing.intervals[i],"))")
#     }
#   }
#
#     study <- study.interest
#     type <- type.interest
#     #determine type of interest and adjust doses accordingly
#     if(type=="mono1"){
#       inds <- which(dose1_>0 & dose2_==0)
#       dose1_c <- dose1_[inds]
#       dose2_c <- dose2_[inds]
#       ndosec <- length(dose1_c)
#       dose1.interestc <- dose1.interest[inds]
#       dose2.interestc <- dose2.interest[inds]
#     }else if(type=="mono2"){
#       inds <- which(dose2_>0 & dose1_==0)
#       dose1_c <- dose1_[inds]
#       dose2_c <- dose2_[inds]
#       ndosec <- length(dose1_c)
#       dose1.interestc <- dose1.interest[inds]
#       dose2.interestc <- dose2.interest[inds]
#     }else if(type=="combi"){
#       inds <- which(dose1_>0 & dose2_>0)
#       dose1_c <- dose1_[inds]
#       dose2_c <- dose2_[inds]
#       ndosec <- length(dose1_c)
#       dose1.interestc <- dose1.interest[inds]
#       dose2.interestc <- dose2.interest[inds]
#     }else{
#       dose1_c <- dose1_
#       dose2_c <- dose2_
#       ndosec <- ndose
#       dose1.interestc <- dose1.interest
#       dose2.interestc <- dose2.interest
#     }
#
#     logdose1_c <- log(dose1_c)
#     logdose2_c <- log(dose2_c)
#
#     toxs <- matrix(NA, nrow = ndosec, ncol = nsamp)
#     toxs1 <- rep(0, nsamp)
#     toxs2 <- rep(0, nsamp)
#     toxs0 <- rep(0, nsamp)
#     #get study-specific posterior samples
#     la1 <- samples[, paste0("log_ab[", study, ",1]")]
#     lb1 <- samples[, paste0("log_ab[", study, ",2]")]
#     la2 <- samples[, paste0("log_ab[", study, ",3]")]
#     lb2 <- samples[, paste0("log_ab[", study, ",4]")]
#     eta <- samples[, paste0("log_ab[", study, ",5]")]
#     explb1 <- exp(lb1)
#     explb2 <- exp(lb2)
#
#     for(d in 1:ndosec){
#       if(dose2_c[d]==0 & dose1_c[d]>0){
#         #treatment mono1
#         toxs[d, ] <- inv_logit(la1 + explb1*logdose1_c[d])
#       }else if(dose1_c[d]==0 & dose2_c[d]>0){
#         #treatment mono2
#         toxs[d, ] <- inv_logit(la2 + explb2*logdose2_c[d])
#       }else if(dose2_c[d]>0 & dose1_c[d]>0){
#         #treatment combination
#         toxs1 <- inv_logit(la1 + explb1*logdose1_c[d])
#         toxs2 <- inv_logit(la2 + explb2*logdose2_c[d])
#         toxs0 <- toxs1 + toxs2 - toxs1*toxs2
#         toxs[d, ] <- inv_logit(logit(toxs0) + eta*dose1_c[d]*dose2_c[d])
#       }else{
#         #no treatment, put in 0
#         toxs[d, ] <- rep(0, times = nsamp)
#       }
#     }
#
#     n.entries <- ninters
#     colnames.res <- c(name.inters)
#
#     res_curr <- matrix(NA, nrow = ndosec, ncol = n.entries)
#     rownames(res_curr) <- paste0(dose1.interestc, "+", dose2.interestc)
#     colnames(res_curr) <-colnames.res
#
#     for(d in 1:ndosec){
#       for(inter in 1:ninters){
#         if(inter==1){
#           res_curr[d, 1] <- mean(ifelse(toxs[d,]<dosing.intervals[1], 1, 0))
#         }else if(inter==ninters){
#           res_curr[d, ninters] <- mean(ifelse(toxs[d,]>=dosing.intervals[ninters-1], 1, 0))
#         }else{
#           res_curr[d, inter] <- mean(ifelse(toxs[d,]>=dosing.intervals[inter-1] &
#                                                       toxs[d,]<dosing.intervals[inter],
#                                                     1, 0))
#         }
#       }
#     }
#
#     if(d.loss){
#       #compute samples from reference toxicities for study
#       #(only required for dynamic loss escalation)
#       ref_toxs <- matrix(NA, nrow = 3, ncol = nsamp)
#       rtoxs0 <- rep(0, nsamp)
#       ref_toxs[1, ] <- inv_logit(la1)
#       ref_toxs[2, ] <- inv_logit(la2)
#       logit_rtoxs0 <- logit(ref_toxs[1,]+ref_toxs[2,]-ref_toxs[1,]*ref_toxs[2,])
#       ref_toxs[3, ] <- inv_logit(logit_rtoxs0 + eta)
#
#       interprobs_ref <- matrix(NA, nrow=3, ncol=ninters)
#       colnames(interprobs_ref) <- name.inters
#       rownames(interprobs_ref) <- c("Ref.Mono1", "Ref.Mono2", "Ref.Combi")
#       for(d in 1:3){
#         for(inter in 1:ninters){
#           if(inter==1){
#             interprobs_ref[d, 1] <- mean(ifelse(ref_toxs[d,]<dosing.intervals[1], 1, 0))
#           }else if(inter==ninters){
#             interprobs_ref[d, ninters] <- mean(ifelse(ref_toxs[d,]>=dosing.intervals[ninters-1], 1, 0))
#           }else{
#             interprobs_ref[d, inter] <- mean(ifelse(ref_toxs[d,]>=dosing.intervals[inter-1] &
#                                                       ref_toxs[d,]<dosing.intervals[inter],
#                                                     1, 0))
#           }
#         }
#       }
#       return(list(
#         "result"=res_curr,
#         "refprobs"=interprobs_ref
#       ))
#
#     }else{
#       return(res_curr)
#     }
#
#
# }

