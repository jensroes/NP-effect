/* 
  Finite mixture model
This model is adapted from Vasishth et al. 2007 
Random intercepts for subj and items
*/
  
data {
  int<lower=1> N;                    // Number of observations
  real y[N];  		            //outcome

  int<lower=0> E;  // Number of experiments
  int<lower=0, upper=E> exps[N];  // Number of experiments
    
  int<lower=1> S;                  //number of subjects
  int<lower=1, upper=S> subj[N];   //subject id
    
  int<lower=1> I;                  //number of items
  int<lower=1, upper=I> items[N];   //items id
}


parameters {
  real<lower=0> delta;			// distribution + extra component
  real<lower=0> sigma;		// residual sd
  real<lower=0> eta; // hyper hyper variance for mixing proportion
  vector<lower=0>[E] sigma_diff;
  
  // Non-centre intercept
  real alpha_mu; // distribution
  real<lower=0> alpha_sigma;
  vector[E] alpha_raw;
  
  real theta; // Hyperhyper mixing proportion 
  vector[E] theta_e;   // Mixing proportion each experiment

  // For random effects
  vector[S] u; //subject intercepts
  real<lower=0> sigma_u;//subj sd
  
  vector[I] w; //items intercepts
  real<lower=0> sigma_w;//items sd
  
}


transformed parameters{
  vector[E] log_theta_e[2];
  vector[E] prob_e = inv_logit(theta_e); //probability of extreme values  
  real prob = inv_logit(theta); // probability of extreme values
  vector<lower=0>[E] sigmap_e = sigma + sigma_diff;
  vector<lower=0>[E] sigma_e = sigma - sigma_diff;
  vector[E] alpha = alpha_mu + alpha_sigma * alpha_raw;
  // inverse logit for prior on theta and log for mixing proprtion
  log_theta_e[1] = log_inv_logit(theta_e);
  log_theta_e[2] = log1m_inv_logit(theta_e);
}

model {
  // Priors
  alpha_mu ~ normal(7, 5);
  alpha_sigma ~ normal(0, 10);
  alpha_raw ~ normal(0, 1);
  delta ~ normal(0, 1);
  sigma_diff ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  
  // Experiment-level mixing proportion
  theta_e ~ normal(theta, eta);

  // Hyper-priors
  // NP level mixing propotion
  theta ~ normal(0, 1);
  eta ~ normal(1, 2); 

  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects
  
  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //items random effects
  
  // Likelihood	
  for(n in 1:N){
    real mu = alpha[exps[n]] + u[subj[n]] + w[items[n]];
    target += log_sum_exp(
      log_theta_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta, sigmap_e[exps[n]]), 
      log_theta_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
  }
}


generated quantities{
  real log_lik[N];
  real y_tilde[N];
  real<lower=0, upper=1> prob_tilde; 
  real alpha2 = alpha_mu + delta;
  vector[E] alpha2exp = alpha + delta;
  
  // likelihood: 
  for(n in 1:N){
    real mu = alpha[exps[n]] + u[subj[n]] + w[items[n]];
    log_lik[n] = log_sum_exp(
      log_theta_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta, sigmap_e[exps[n]]), 
      log_theta_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
    prob_tilde = bernoulli_rng(prob_e[exps[n]]);
    if(prob_tilde){
      y_tilde[n] = lognormal_rng(mu + delta, sigmap_e[exps[n]]);
    }
    else{
      y_tilde[n] = lognormal_rng(mu, sigma_e[exps[n]]);
    }
  }
}
