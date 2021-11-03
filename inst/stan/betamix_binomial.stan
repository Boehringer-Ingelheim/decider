/*
Stanmodel to evaluate simple binomial beta-mixture models.

*/
data{
  //number of observations (pooled to the same dose,
  //only for convenience included)
  int nobs;
  //number of patients per observation
  int n[nobs];
  //number of responders/DLTs
  int r[nobs];
  //number of mixture comps
  int num_mix_comps;
  //mixture weights
  real mix_probs[num_mix_comps];
  //shape parameter 1
  real mix_a[num_mix_comps];
  //shape parameter 2
  real mix_b[num_mix_comps];
}

transformed data{
  //log mix weights for stan implementation
  vector[num_mix_comps] log_mix_probs;
  for(i in 1:num_mix_comps){
    log_mix_probs[i] = log(mix_probs[i]);
  }

}

parameters{
  //response rate (pooled across observations)
  real<lower=0, upper=1> p;
}

model{
  //for saving the individual beta-distributions
  vector[num_mix_comps] prior_comps;

  //Model for treatment per obs
  for(i in 1:nobs){
    r[i] ~ binomial(n[i], p);
  }

  //prior, n components
  for(i in 1:num_mix_comps){
    //ith component is beta(a[i], b[i])
    prior_comps[i] = beta_lpdf(p | mix_a[i], mix_b[i]);
  }
  //mixture on log scale
  //note that adding the log mix-weights translates
  //to component-wise multiplication on non-log scale.
  target += log_sum_exp(prior_comps + log_mix_probs);

}
//Note: posterior weights can be calculated analytically, no need to sample here.
//this stanmodel is mostly used a RNG from the posterior for convenience.




