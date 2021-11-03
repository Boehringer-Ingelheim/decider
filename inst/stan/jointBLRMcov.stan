/*Stan model for joint BLRM with covariates
--------------------------------------------------------------------------------
  Implements the joint BLRM as described in Neuenschwander et al., 2016,
  "On the use of co-data in clinical trials".
  A non-centered parametrization  is implemented by obtaining
  multivariate normals via multiplication with cholesky factors.
  The cholesky decomposition is implemented by hand, as it is
  available analytically in the required 2x2-case.
*/
functions{

  /*counts mono observations based on input dose levels
    Note: first input vector signals the component to be counted*/
  int count_n_mono(vector dose_1, vector dose_2, int n_obs){
    int res = 0;
    for(i in 1:n_obs){
      if(dose_1[i]>0 && dose_2[i]==0){
        res+=1;
      }
    }
    return res;
  }

  /*Computes permutation of input data, so that the first n_obs1 observations
    are mono1, the subsequent n_obs2 observations are mono2, and the remaining
    ones are combination therapy.
    Returns matrix with two rows, first row is the permutation for sorting, and
    second row contains the inverse permutation (to reverse sorted input to
    normal order).*/
  int[,] sort_idx(vector dose_1, vector dose_2,
                 int n_obs, int n_obs1, int n_obs2)
  {
    int res[2, n_obs] = rep_array(0, 2, n_obs);
    //n_obs1/n_obs2 allow to compute offsets for sorting by counting
    int cnt1 = 0;
    int cnt2 = 0;
    int cnt = 0;
    //loop over input and save correct placement
    for(i in 1:n_obs){
      if(dose_1[i]>0 && dose_2[i]==0){
        res[1, cnt1+1] = i;
        res[2, i] = cnt1+1;
        cnt1 += 1;
      }else if(dose_1[i]==0 && dose_2[i]>0){
        res[1, n_obs1 + 1 + cnt2] = i;
        res[2, i] = n_obs1 + 1 + cnt2;
        cnt2 += 1;
      }else if(dose_1[i]>0 && dose_2[i]>0){
        res[1, n_obs1 + n_obs2 + 1 + cnt] = i;
        res[2, i] = n_obs1 + n_obs2 + 1 + cnt;
        cnt += 1;
      }
    }
    return res;
  }

}


data{
  //number of observations/cohorts
  int<lower=0> n_obs;
  //number of studies
  int<lower=0> n_studies;
  //number of patients for each cohort
  int<lower=0> n[n_obs];
  //number of DLTs for each cohort
  int<lower=0> r[n_obs];
  //binary covariate per cohort
  int<lower=0, upper=1> c[n_obs];
  //study number for cohorts
  int<lower=1> s[n_obs];
  //indicates whether a MAP prior is computed
  int<lower=0, upper=1> doMAP;
  //indicates whether a saturating interaction term is used
  int<lower=0, upper=1> saturating;

  //indicates whether the covariate is one-sided or two-sided
  int<lower=0, upper=1> twoside1;
  int<lower=0, upper=1> twoside2;

  //reference doses
  vector<lower=0>[2] dose_c;
  //dose levels component 1 and 2 for each cohort
  vector<lower=0>[n_obs] dose_1;
  vector<lower=0>[n_obs] dose_2;

  /*hyper priors
    Notation and order of entries:
    mu =  (mu_alpha1,  mu_beta1,  mu_alpha2,  mu_beta2,  mu_eta, mu_gamma1, mu_gamma2)
    tau = (tau_alpha1, tau_beta1, tau_alpha2, tau_beta2, tau_eta, tau_gamma1, tau_gamma2)
  */
  //mean of hyper SD tau
  vector[7] mean_tau;
  //sd's of hyper SD tau
  vector<lower=0>[7] sd_tau;
  //mean of hyper mean mu
  vector[7] mean_mu;
  //mean of hyper sd mu
  vector<lower=0>[7] sd_mu;
}

transformed data{
  //internally generates a study without observations for MAP prior
  int<lower=1> num_s = doMAP? n_studies+1 : n_studies;
  //count number of mono observations
  int<lower=0, upper=n_obs> n_obs1 = count_n_mono(dose_1, dose_2, n_obs);
  int<lower=0, upper=n_obs> n_obs2 = count_n_mono(dose_2, dose_1, n_obs);
  //compute sort indices (only done once per call to stan for efficiency)
  int srt_idx[2, n_obs] = sort_idx(dose_1, dose_2, n_obs, n_obs1, n_obs2);
  //sort by applying computed sorting permutation
  int n_srt[n_obs] = n[srt_idx[1, 1:n_obs]];
  int r_srt[n_obs] = r[srt_idx[1, 1:n_obs]];
  int s_srt[n_obs] = s[srt_idx[1, 1:n_obs]];
  int c_srt[n_obs] = c[srt_idx[1, 1:n_obs]];
  vector[n_obs] vc_srt = to_vector(c_srt);
  //doses are also rescaled by reference dose after sorting
  vector[n_obs] dose_1_srt = dose_1[srt_idx[1, 1:n_obs]]/dose_c[1];
  vector[n_obs] ldose_1_srt = log(dose_1_srt);
  vector[n_obs] dose_2_srt = dose_2[srt_idx[1, 1:n_obs]]/dose_c[2];
  vector[n_obs] ldose_2_srt = log(dose_2_srt);
}

parameters{
  //hyper SDs
  real<lower=0> tau_1a;
  real<lower=0> tau_1b;
  real<lower=0> tau_2a;
  real<lower=0> tau_2b;
  real<lower=0> tau_eta;
  real<lower=0> tau_1c;
  real<lower=0> tau_2c;
  //correlation coefficients
  real<lower=-1, upper=1> rho12;
  real<lower=-1, upper=1> rho34;

  /*For non-centered parametrization:
    Sample only raw standard normal variables. These are later transformed to
    bivariate normals by multiplying with cholesky factor*/
  //matrix for log(alpha_ij), log(beta_ij) and eta_j (for comp i, study j)
  matrix[num_s, 7] log_ab_raw;
  //for hyper means
  real mu_raw[7];
}
transformed parameters{
  real mu_1a;
  real mu_1b;
  real mu_2a;
  real mu_2b;
  real mu_eta;
  real mu_1c;
  real mu_2c;

  matrix[num_s,7] log_ab;
  vector<lower=0, upper=1>[n_obs] p_srt;
  vector<lower=0, upper=1>[n_obs-n_obs1-n_obs2] p_2;
  vector<lower=0, upper=1>[n_obs-n_obs1-n_obs2] p_1;
  vector<lower=0, upper=1>[n_obs-n_obs1-n_obs2] p_0;

  //transform raw hyper means to correct distribution
  mu_1a = mean_mu[1] + sd_mu[1]*mu_raw[1];
  mu_1b = mean_mu[2] + sd_mu[2]*mu_raw[2];
  mu_2a = mean_mu[3] + sd_mu[3]*mu_raw[3];
  mu_2b = mean_mu[4] + sd_mu[4]*mu_raw[4];
  mu_eta = mean_mu[5] + sd_mu[5]*mu_raw[5];
  mu_1c = mean_mu[6] + sd_mu[6]*mu_raw[6];
  mu_2c = mean_mu[7] + sd_mu[7]*mu_raw[7];

  /*Hard-coded matrix multiplication with lower cholesky factor
    of covariance matrix. This can be done without saving the
    cholesky factor itself, as it is available analytically.
    The following means:
    log_ab = mu + L*log_ab_raw,
    where L is a lower triangular matrix with L*L^T=Sigma,
    for a covariance matrix Sigma.
    Note: For general
    Sigma = tau_1^2           rho*tau_1*tau_2
            rho*tau_1*tau_2   tau_2^2
    the lower cholesky factor is
    L =     tau_1         0
            tau_2*rho     tau_2*squareroot(1-rho^2)
    */
  log_ab[1:num_s,1] = mu_1a + tau_1a*log_ab_raw[1:num_s, 1];
  log_ab[1:num_s,2] = mu_1b + tau_1b*rho12*log_ab_raw[1:num_s, 1] +
                      tau_1b*sqrt(1-square(rho12))*log_ab_raw[1:num_s, 2];
  log_ab[1:num_s,3] = mu_2a + tau_2a*log_ab_raw[1:num_s, 3];
  log_ab[1:num_s,4] = mu_2b + tau_2b*rho34*log_ab_raw[1:num_s, 3] +
                      tau_2b*sqrt(1-square(rho34))*log_ab_raw[1:num_s, 4];
  log_ab[1:num_s,5] = mu_eta + tau_eta*log_ab_raw[1:num_s, 5];
  log_ab[1:num_s,6] = mu_1c + tau_1c*log_ab_raw[1:num_s, 6];
  log_ab[1:num_s,7] = mu_2c + tau_2c*log_ab_raw[1:num_s, 7];



  //toxicity models for mono and combination treatment are vectorized
  if(n_obs1>0){
    //treatments mono 1
    if(twoside1){
      p_srt[1:n_obs1] = inv_logit(log_ab[s_srt[1:n_obs1],1] +
          log_ab[s_srt[1:n_obs1], 6].*vc_srt[1:n_obs1]+
          exp(log_ab[s_srt[1:n_obs1],2]).*
          ldose_1_srt[1:n_obs1]);
    }else{
      //cov1 is one-sided
      p_srt[1:n_obs1] = inv_logit(log_ab[s_srt[1:n_obs1],1] +
          exp(log_ab[s_srt[1:n_obs1], 6]).*vc_srt[1:n_obs1]+
          exp(log_ab[s_srt[1:n_obs1],2]).*
          ldose_1_srt[1:n_obs1]);
    }
  }
  if(n_obs2>0){
    //treatments mono 2
    if(twoside2){
      p_srt[(n_obs1+1):(n_obs1+n_obs2)] =
          inv_logit(log_ab[s_srt[(n_obs1+1):(n_obs1 + n_obs2)],3] +
            log_ab[s_srt[(n_obs1+1):(n_obs1 + n_obs2)], 7].*vc_srt[(n_obs1+1):(n_obs1 + n_obs2)]+
            exp(log_ab[s_srt[(n_obs1+1): (n_obs1 + n_obs2)],4]).*
            ldose_2_srt[(n_obs1+1): (n_obs1 + n_obs2)]);
    }else{
      //cov one-sided
      p_srt[(n_obs1+1):(n_obs1+n_obs2)] =
          inv_logit(log_ab[s_srt[(n_obs1+1):(n_obs1 + n_obs2)],3] +
            exp(log_ab[s_srt[(n_obs1+1):(n_obs1 + n_obs2)], 7]).*vc_srt[(n_obs1+1):(n_obs1 + n_obs2)]+
            exp(log_ab[s_srt[(n_obs1+1): (n_obs1 + n_obs2)],4]).*
            ldose_2_srt[(n_obs1+1): (n_obs1 + n_obs2)]);
    }
  }
  if(n_obs-n_obs1-n_obs2>0){
    //treatments combination
    if(twoside2){
      p_2[1 : (n_obs-n_obs1-n_obs2)] =
          inv_logit(log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs],3] +
              log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs], 7].*
              vc_srt[(n_obs1 + n_obs2 + 1) : n_obs]+
              exp(log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs],4]).*
              ldose_2_srt[(n_obs1 + n_obs2 + 1) : n_obs]);
    }else{
      //cov2 one-sided
      p_2[1 : (n_obs-n_obs1-n_obs2)] =
          inv_logit(log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs],3] +
              exp(log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs], 7]).*
              vc_srt[(n_obs1 + n_obs2 + 1) : n_obs] +
              exp(log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs],4]).*
              ldose_2_srt[(n_obs1 + n_obs2 + 1) : n_obs]);
    }
    if(twoside1){
      p_1[1 : (n_obs-n_obs1-n_obs2)] =
          inv_logit(log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs],1] +
              log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs], 6].*
              vc_srt[(n_obs1 + n_obs2 + 1) : n_obs] +
              exp(log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs],2]).*
              ldose_1_srt[(n_obs1 + n_obs2 + 1) : n_obs]);
    }else{
      p_1[1 : (n_obs-n_obs1-n_obs2)] =
          inv_logit(log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs],1] +
              exp(log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs], 6]).*
              vc_srt[(n_obs1 + n_obs2 + 1) : n_obs] +
              exp(log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs],2]).*
              ldose_1_srt[(n_obs1 + n_obs2 + 1) : n_obs]);
    }
    p_0[1 :(n_obs-n_obs1-n_obs2)] = p_1[1 : (n_obs-n_obs1-n_obs2)] +
        p_2[1 : (n_obs-n_obs1-n_obs2)] -
        p_1[1 : (n_obs-n_obs1-n_obs2)] .*
        p_2[1 : (n_obs-n_obs1-n_obs2)];
    if(saturating){
      p_srt[(n_obs1 + n_obs2 + 1) : n_obs] =
          inv_logit(logit(p_0[1 : (n_obs-n_obs1-n_obs2)]) +
              (2*log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs],5].*
              (dose_1_srt[(n_obs1 + n_obs2 + 1) : n_obs].*
              dose_2_srt[(n_obs1 + n_obs2 + 1) : n_obs])./
              (1 + dose_1_srt[(n_obs1 + n_obs2 + 1) : n_obs].*
              dose_2_srt[(n_obs1 + n_obs2 + 1) : n_obs])
              ));
    }else{
      p_srt[(n_obs1 + n_obs2 + 1) : n_obs] =
          inv_logit(logit(p_0[1 : (n_obs-n_obs1-n_obs2)]) +
              log_ab[s_srt[(n_obs1 + n_obs2 + 1) : n_obs],5].*
              dose_1_srt[(n_obs1 + n_obs2 + 1) : n_obs].*
              dose_2_srt[(n_obs1 + n_obs2 + 1) : n_obs]);
    }
  }

}

model{
  //priors for hyper means (non-centered)
  mu_raw ~  std_normal();
  //priors for hyper SD
  tau_1a ~ lognormal(mean_tau[1], sd_tau[1]);
  tau_1b ~ lognormal(mean_tau[2], sd_tau[2]);
  tau_2a ~ lognormal(mean_tau[3], sd_tau[3]);
  tau_2b ~ lognormal(mean_tau[4], sd_tau[4]);
  tau_eta ~ lognormal(mean_tau[5], sd_tau[5]);
  tau_1c ~ lognormal(mean_tau[6], sd_tau[6]);
  tau_2c ~ lognormal(mean_tau[7], sd_tau[7]);
  //priors for correlation coefficients
  rho12 ~ uniform(-1,1);
  rho34 ~ uniform(-1,1);
  //priors for regression parameters (non-centered)
  for(k in 1:num_s){
    log_ab_raw[k, 1:7] ~ std_normal();
  }
  //binomial likelihood
  r_srt ~ binomial(n_srt, p_srt);
}

generated quantities{
  //just to provide the sorted toxicity parameters as output
  vector<lower=0, upper=1>[n_obs] p = p_srt[srt_idx[2,1:n_obs]];

}
