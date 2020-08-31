/* 
  Finite mixture model
  This model is adapted from Vasishth et al. 2007 
  Random intercepts for subj and items
  Latent mixture proportion based on indv mixture proportion for each experiment (and condition) 
  implemented as hyper parameter.
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
  real<lower=0> delta;		// distribution + extra component
  real<lower=0> sigma;		// residual sd
  real theta;
  vector<lower=0>[E] sigma_diff;

  // Non-centre intercept
  real alpha_mu; // distribution
  real<lower=0> alpha_sigma;
  vector[E] alpha_raw;

  // For random effects
  vector[S] u; //subject intercepts
  real<lower=0> sigma_u;//subj sd
  
  vector[I] w; //items intercepts
  real<lower=0> sigma_w;//items sd
  
}


transformed parameters{
  real log_theta[2];
  real prob = inv_logit(theta);
  vector<lower=0>[E] sigmap_e = sigma + sigma_diff;
  vector<lower=0>[E] sigma_e = sigma - sigma_diff;
  vector[E] alpha = alpha_mu + alpha_sigma * alpha_raw;
  
  // inverse logit for prior on theta and log for mixing proprtion
  log_theta[1] = log_inv_logit(theta);
  log_theta[2] = log1m_inv_logit(theta);
}

model {
  // Priors
	alpha_mu ~ normal(7, 4);
	alpha_sigma ~ normal(0, 10);
	alpha_raw ~ normal(0, 2);
	
  delta ~ normal(0, 1);
  sigma_diff ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  
  // Hyper-priors
  theta ~ normal(-.5, 1);

  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects
  
  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //items random effects
  
  // Likelihood	
  for(n in 1:N){
    real mu = alpha[exps[n]] + u[subj[n]] + w[items[n]];
    target += log_sum_exp(
      log_theta[1] + lognormal_lpdf(y[n] | mu + delta, sigmap_e[exps[n]]), 
      log_theta[2] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
  }
}


generated quantities{
  real log_lik[N];
  real y_tilde[N];
  real<lower=0,upper=1> prob_tilde; 
  real alpha2 = alpha_mu + delta;
  vector[E] alpha2exp = alpha + delta;
 
  // likelihood: 
  for(n in 1:N){
    real mu = alpha[exps[n]] + u[subj[n]] + w[items[n]];
    log_lik[n] = log_sum_exp(
      log_theta[1] + lognormal_lpdf(y[n] | mu + delta, sigmap_e[exps[n]]), 
      log_theta[2] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
  	prob_tilde = bernoulli_rng(prob); 
    if(prob_tilde) { 
      y_tilde[n] = lognormal_rng(mu + delta, sigmap_e[exps[n]]);
    }
    else{
      y_tilde[n] = lognormal_rng(mu, sigma_e[exps[n]]);
    }
  }
}
