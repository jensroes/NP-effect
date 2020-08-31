/* 
Finite mixture model
This model is adapted from Vasishth et al. 2007 
Random intercepts for subj and items
Mixture on simple and complex condition
*/

data {
	int<lower=1> N;                    // Number of observations
	real y[N];  		            //outcome

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
	
	real theta; // mixing proportion

	// For random effects
	vector[S] u; //subject intercepts
	real<lower=0> sigma_u;//subj sd

	vector[I] w; //items intercepts
	real<lower=0> sigma_w;//items sd
}


transformed parameters{
	vector[2] log_prob;
	real alpha = alpha_mu + alpha_sigma*alpha_raw;
	real<lower=0> sigmap_e = sigma + sigma_diff;
	real<lower=0> sigma_e = sigma - sigma_diff;
  real prob = inv_logit(theta);

  log_prob[1] = log_inv_logit(theta);
  log_prob[2] = log1m_inv_logit(theta);
}

model {
  // Priors
	alpha_mu ~ normal(6, 3);
	alpha_sigma ~ normal(0, 10);
	alpha_raw ~ normal(0, 1);

  delta ~ normal(0, 1);
  sigma_diff ~ normal(0,1);
  sigma ~ cauchy(0, 2.5);
  
  theta ~ normal(0, 2);

  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //items random effects


  // Likelihood	
	for(n in 1:N){
    real mu = alpha + u[subj[n]] + w[items[n]];
    target += log_sum_exp(
        log_prob[1] + lognormal_lpdf(y[n] | mu + delta, sigmap_e), 
        log_prob[2] + lognormal_lpdf(y[n] | mu, sigma_e));
	}
}


generated quantities{
  real log_lik[N];
  real y_tilde[N];
  real<lower=0,upper=1> prob_tilde; 
  real alpha2 = alpha + delta;
  
  // likelihood: 
  for(n in 1:N){
    real mu = alpha + u[subj[n]] + w[items[n]];
  	log_lik[n] = log_sum_exp(
  	      		log_prob[1] + lognormal_lpdf(y[n] | mu + delta, sigmap_e), 
  	    		  log_prob[2] + lognormal_lpdf(y[n] | mu, sigma_e));
    prob_tilde = bernoulli_rng(prob); 
    if(prob_tilde) { 
      y_tilde[n] = lognormal_rng(mu + delta, sigmap_e);
    }
    else{
      y_tilde[n] = lognormal_rng(mu, sigma_e);
    }
  }
}
