/*
lognormal distribution 
by-subject and by-item intercepts (no random slopes; see Vasishth et al. 2017)
model for NP complexity effect as typically assessed in the literature
*/

data {
	int<lower=1> N;                    // Number of observations
	real<lower=0> y[N];  		            //outcome
	int<lower=0, upper=1> nounphrase[N];  //predictor
	
	int<lower=1> S;                  //number of subjects
	int<lower=1, upper=S> subj[N];   //subject id

	int<lower=1> I;                  //number of items
	int<lower=1, upper=I> items[N];   //items id
}

parameters {
  real beta;
	real<lower=0> sigma;		// residual sd
	real<lower=0> sigma_diff;

  // Non-centre intercept
  real alpha_mu; // distribution
  real<lower=0> alpha_sigma;
  real alpha_raw;

	// For random effects
	vector[S] u; //subject intercepts
	real<lower=0> sigma_u;//subj sd

	vector[I] w; //items intercepts
	real<lower=0> sigma_w;//items sd
}

transformed parameters{
  real<lower=0> sigmap_e = sigma + sigma_diff;
	real<lower=0> sigma_e = sigma - sigma_diff;
	real alpha = alpha_mu + alpha_sigma*alpha_raw;
}

model {
	// Priors
	alpha_mu ~ normal(6, 2);
	alpha_sigma ~ normal(0, 10);
	alpha_raw ~ normal(0, 1);
	beta ~ normal(0, 1);
	sigma ~ cauchy(0, 2.5); 
  sigma_diff ~ normal(0, 1);

	// REs priors
	sigma_u ~ normal(0,2.5);
	u ~ normal(0, sigma_u); //subj random effects

	sigma_w ~ normal(0,2.5);
	w ~ normal(0, sigma_w); //items random effects


	// Likelihood
	for(n in 1:N){
	  real mu = alpha + u[subj[n]] + w[items[n]];
    if(nounphrase[n]==0){
  	  y[n] ~ lognormal(mu, sigma_e);
    }
    if(nounphrase[n]==1){
  	  y[n] ~ lognormal(mu + beta, sigmap_e);
    }
	}
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_tilde;
  real alpha2 = alpha + beta;

  for (n in 1:N) {
	  real mu = alpha + u[subj[n]] + w[items[n]];
    if(nounphrase[n]==0){
      log_lik[n] = lognormal_lpdf(y[n] | mu, sigma_e);
      y_tilde[n] = lognormal_rng(mu, sigma_e);
    }
    if(nounphrase[n]==1){
      log_lik[n] = lognormal_lpdf(y[n] | mu + beta, sigmap_e);
      y_tilde[n] = lognormal_rng(mu + beta, sigmap_e);
    }
  }
}


