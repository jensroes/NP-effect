/*
lognormal distribution 
by-subject and by-item intercepts (no random slopes; see Vasishth et al. 2017)
model for NP complexity effect as typically assessed in the literature
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
	real<lower=0> sigma;		// residual sd

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
  vector[N] mu;
	real alpha = alpha_mu + alpha_sigma*alpha_raw;

  for(n in 1:N){
    mu[n] = alpha + u[subj[n]] + w[items[n]];
  }
}


model {
	// Priors
	alpha_mu ~ normal(6, 3);
	alpha_sigma ~ normal(0, 10);
	alpha_raw ~ normal(0, 1);
	
	// REs priors
	sigma_u ~ normal(0,2.5);
	u ~ normal(0, sigma_u); //subj random effects

	sigma_w ~ normal(0,2.5);
	w ~ normal(0, sigma_w); //items random effects

	sigma ~ cauchy(0, 2.5); 

	// Likelihood
	y ~ lognormal(mu, sigma);

}

generated quantities {
  vector[N] log_lik;
  vector[N] y_tilde;

  for (n in 1:N) {
    log_lik[n] = lognormal_lpdf(y[n] | mu[n], sigma);
    y_tilde[n] = lognormal_rng(mu[n], sigma);
  }
}
