/* 
  Finite mixture model
  This model is adapted from Vasishth et al. 2007 
  Random intercepts for subj and items
  Mixture on simple and complex condition
  Latent mixture proportion based on indv mixture proportion for each experiment (and condition) 
  implemented as hyper parameter.
*/
  
data {
  int<lower=1> N;                    // Number of observations
  real y[N];  		            //outcome
  int<lower=1, upper=2> nounphrase[N];  //predictor
  
  int<lower=0> E;  // Number of experiments
  int<lower=0, upper=E> exps[N];  // Number of experiments

  int<lower=1> S;                  //number of subjects
  int<lower=1, upper=S> subj[N];   //subject id
    
  int<lower=1> I;                  //number of items
  int<lower=1, upper=I> items[N];   //items id
}


parameters {
  real<lower=0> sigma;		// residual sd
  real<lower=0> tau;      // variance between exps
  vector<lower=0>[E] sigma_diff;

  real theta_simple;         // Hyperparameter for mixing proportion (across experiments for each complex and simple)
  real theta_complex;
  vector[E] theta_simple_e;   // Mixing proportion by-condition and by-experiment
  vector[E] theta_complex_e;

  // Non-centre intercept
  real alpha_mu; // distribution
  real<lower=0> alpha_sigma;
  vector[E] alpha_raw;

  // Hyper parameters
  real<lower=0> delta_mu; // effect across all experiments
  real<lower=0> delta_tau; // variation between experiments

  vector<lower=0>[E] delta;
//  vector[E] delta_eta;

  // For random effects
  vector[S] u; //subject intercepts
  real<lower=0> sigma_u;//subj sd
  
  vector[I] w; //items intercepts
  real<lower=0> sigma_w;//items sd

}


transformed parameters{
  vector[E] log_theta_simple_e[2];
  vector[E] log_theta_complex_e[2];
  vector[E] prob_simple_e = inv_logit(theta_simple_e); //probability of extreme values by exp
  vector[E] prob_complex_e  = inv_logit(theta_complex_e); 
//  real prob = inv_logit(theta);
  real prob_simple = inv_logit(theta_simple); //probability of extreme values
  real prob_complex = inv_logit(theta_complex); 
  vector<lower=0>[E] sigmap_e = sigma + sigma_diff;
  vector<lower=0>[E] sigma_e = sigma - sigma_diff;
  vector[E] alpha = alpha_mu + alpha_sigma * alpha_raw;
//  vector<lower=0>[E] delta = delta_mu + delta_tau * delta_eta; // effect for each experiment

  // inverse logit for prior on theta and log for mixing proprtion
  log_theta_simple_e[1] = log_inv_logit(theta_simple_e);
  log_theta_simple_e[2] = log1m_inv_logit(theta_simple_e);
  log_theta_complex_e[1] = log_inv_logit(theta_complex_e);
  log_theta_complex_e[2] = log1m_inv_logit(theta_complex_e);
}

model {
  // Priors
	alpha_mu ~ normal(6, 2);
	alpha_sigma ~ normal(0, 2);
	alpha_raw ~ normal(0, 1);
	
  // Hyper priors	
	delta ~ normal(delta_mu, delta_tau);
  delta_tau ~ cauchy(0, 2.5);
  delta_mu ~ normal(0, 5);
//  delta_eta ~ normal(0, 1);

  sigma_diff ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  
  // Hyper-priors
  tau ~ normal(0, 1);
  theta_simple_e ~ normal(theta_simple, tau);
  theta_complex_e ~ normal(theta_complex, tau);
  
  theta_simple ~ normal(0, 1);
  theta_complex ~ normal(0, 1);

  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects
  
  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //items random effects
  
  // Likelihood	
  for(n in 1:N){
    real mu = alpha[exps[n]] + u[subj[n]] + w[items[n]];
    if(nounphrase[n]==1){
      target += log_sum_exp(
        log_theta_simple_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta[exps[n]], sigmap_e[exps[n]]), 
        log_theta_simple_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
    }
    if(nounphrase[n]==2){
      target += log_sum_exp(
        log_theta_complex_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta[exps[n]], sigmap_e[exps[n]]), 
        log_theta_complex_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
    }
  }
}


generated quantities{
  real log_lik[N];
  real y_tilde[N];
  real<lower=0,upper=1> prob_tilde; 
//  real alpha2 = alpha_mu + delta_mu;
  vector[E] alpha2exp = alpha + delta;
  real prob_diff = prob_complex - prob_simple;
  vector[E] prob_diff_e = prob_complex_e - prob_simple_e;

  // likelihood: 
  for(n in 1:N){
    real mu = alpha[exps[n]] + u[subj[n]] + w[items[n]];
    if(nounphrase[n] == 1){// simple NPs
      log_lik[n] = log_sum_exp(
        log_theta_simple_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta[exps[n]], sigmap_e[exps[n]]), 
        log_theta_simple_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
    	prob_tilde = bernoulli_rng(prob_simple_e[exps[n]]); 
      if(prob_tilde) { 
        y_tilde[n] = lognormal_rng(mu + delta[exps[n]], sigmap_e[exps[n]]);
      }
      else{
        y_tilde[n] = lognormal_rng(mu, sigma_e[exps[n]]);
      }
    }
    if(nounphrase[n] == 2){// complex NPs
      log_lik[n] = log_sum_exp(
        log_theta_complex_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta[exps[n]], sigmap_e[exps[n]]), 
        log_theta_complex_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
    	prob_tilde = bernoulli_rng(prob_complex_e[exps[n]]); 
      if(prob_tilde) { 
        y_tilde[n] = lognormal_rng(mu + delta[exps[n]], sigmap_e[exps[n]]);
      }
      else{
        y_tilde[n] = lognormal_rng(mu, sigma_e[exps[n]]);
      }
    }
  }
}
