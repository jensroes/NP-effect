/* 
Finite mixture model
This model is adapted from Vasishth et al. 2007 
Random intercepts for subj and items
Mixture on simple and complex condition
*/

data {
	int<lower=1> N;                    // Number of observations
	real y[N];  		            //outcome
	int<lower=1, upper=2> nounphrase[N];  //predictor
	
	int<lower=1> S;                  //number of subjects
	int<lower=1, upper=S> subj[N];   //subject id

	int<lower=1> I;                  //number of items
	int<lower=1, upper=I> items[N];   //items id
}


parameters {
  // Non-centre intercept
  real alpha_mu; // distribution
  real<lower=0> alpha_sigma;
  real alpha_raw;

	real<lower=0> delta;			// distribution + extra component
	
	real<lower=0> sigma;		// residual sd
	real<lower=0> sigma_diff;
	
	real theta_simple; // mixing proportion
	real theta_complex;
//	real theta;
	
//	real<lower = 0> tau; // variance for mixing proportion

	// For random effects
	vector[S] u; //subject intercepts
	real<lower=0> sigma_u;//subj sd

	vector[I] w; //items intercepts
	real<lower=0> sigma_w;//items sd
}


transformed parameters{
	vector[2] log_prob_simple;
	vector[2] log_prob_complex;
	real alpha = alpha_mu + alpha_sigma*alpha_raw;
	real<lower=0> sigmap_e = sigma + sigma_diff;
	real<lower=0> sigma_e = sigma - sigma_diff;
	real prob_simple = inv_logit(theta_simple); //probability of extreme values
	real prob_complex = inv_logit(theta_complex); //probability of extreme values
//  real prob = inv_logit(theta);

  log_prob_simple[1] = log_inv_logit(theta_simple);
  log_prob_simple[2] = log1m_inv_logit(theta_simple);
  log_prob_complex[1] = log_inv_logit(theta_complex);
  log_prob_complex[2] = log1m_inv_logit(theta_complex);
}

model {
  // Priors
	alpha_mu ~ normal(6, 3);
	alpha_sigma ~ normal(0, 10);
	alpha_raw ~ normal(0, 1);

  delta ~ normal(0, 1);
  sigma_diff ~ normal(0,1);
  sigma ~ cauchy(0, 2.5);
  
  theta_simple ~ normal(0, 2);
  theta_complex ~ normal(0, 2);
  
//  theta ~ normal(0, 3);
//  tau ~ normal(0, 10);
  
    // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //items random effects


  // Likelihood	
	for(n in 1:N){
    real mu = alpha + u[subj[n]] + w[items[n]];
	  if(nounphrase[n]==1){
      target += log_sum_exp(
        log_prob_simple[1] + lognormal_lpdf(y[n] | mu + delta, sigmap_e), 
        log_prob_simple[2] + lognormal_lpdf(y[n] | mu, sigma_e));
    }
    if(nounphrase[n]==2){
     target += log_sum_exp(
       log_prob_complex[1] + lognormal_lpdf(y[n] | mu + delta, sigmap_e), 
       log_prob_complex[2] + lognormal_lpdf(y[n] | mu, sigma_e));
    }
	}
}


generated quantities{
  real log_lik[N];
  real y_tilde[N];
  real<lower=0,upper=1> prob_tilde; 
  real alpha2 = alpha + delta;
  real prob_diff = prob_complex - prob_simple;
  
  // likelihood: 
  for(n in 1:N){
    real mu = alpha + u[subj[n]] + w[items[n]];
    if(nounphrase[n] == 1){// simple NPs
    	log_lik[n] = log_sum_exp(
    	      		log_prob_simple[1] + lognormal_lpdf(y[n] | mu + delta, sigmap_e), 
    	    		  log_prob_simple[2] + lognormal_lpdf(y[n] | mu, sigma_e));
    	prob_tilde = bernoulli_rng(prob_simple); 
      if(prob_tilde) { 
        y_tilde[n] = lognormal_rng(mu + delta, sigmap_e);
      }
      else{
        y_tilde[n] = lognormal_rng(mu, sigma_e);
      }
    }
    if(nounphrase[n] == 2){// complex NPs
    	log_lik[n] = log_sum_exp(
    	      		log_prob_complex[1] + lognormal_lpdf(y[n] | mu + delta, sigmap_e), 
    	    		  log_prob_complex[2] + lognormal_lpdf(y[n] | mu, sigma_e));
    	prob_tilde = bernoulli_rng(prob_complex); 
      if(prob_tilde) { 
        y_tilde[n] = lognormal_rng(mu + delta, sigmap_e);
      }
      else{
        y_tilde[n] = lognormal_rng(mu, sigma_e);
      }
    }
  }
}
