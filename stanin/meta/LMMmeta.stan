/*
lognormal distribution 
by-subject and by-item intercepts (no random slopes; see Vasishth et al. 2017)
model for NP complexity effect as typically assessed in the literature
*/

data {
	int<lower=1> N;                    // Number of observations
	real<lower=0> y[N];  		            //outcome
	real nounphrase[N];  //predictor
	
	int<lower=0> E;  // Number of experiments
	int<lower=0, upper=E> exps[N];  // Number of experiments

	int<lower=1> S;                  //number of subjects
	int<lower=1, upper=S> subj[N];   //subject id

	int<lower=1> I;                  //number of items
	int<lower=1, upper=I> items[N];   //items id
}


parameters {
	vector<lower=0>[E] sigma;		// residual sd

  // Non-centre intercept
  real alpha_mu; // distribution
  real<lower=0> alpha_sigma;
  vector[E] alpha_raw;

  // Hyper parameters
  real beta; // effect across all experiments
  real<lower=0> tau; // variation between experiments
  vector[E] eta;

	// For random effects
	vector[S] u; //subject intercepts
	real<lower=0> sigma_u;//subj sd

	vector[I] w; //items intercepts
	real<lower=0> sigma_w;//items sd

}


transformed parameters{
  vector[N] mu;
  vector[N] sigma_exp;
  vector[E] betae = beta + tau * eta; // effect for each experiment
	vector[E] alpha = alpha_mu + alpha_sigma*alpha_raw;
  
  for(n in 1:N){
    mu[n] = alpha[exps[n]] + nounphrase[n]*betae[exps[n]] + u[subj[n]] + w[items[n]];
    sigma_exp[n] = sigma[exps[n]];
  }
}


model {
	// Priors
	alpha_mu ~ normal(7, 3);
	alpha_sigma ~ cauchy(0, 1);
	alpha_raw ~ normal(0, 1);
	sigma ~ cauchy(0, 2.5); 

  // Hyper priors	
	beta ~ normal(0, 10);
  tau ~ cauchy(0, 1);
  eta ~ normal(0, 1);
  
	// REs priors
	sigma_u ~ normal(0,2.5);
	u ~ normal(0, sigma_u); //subj random effects

	sigma_w ~ normal(0,2.5);
	w ~ normal(0, sigma_w); //items random effects
	
	// Likelihood
	y ~ lognormal(mu, sigma_exp);

}

generated quantities {
  vector[N] log_lik;
  vector[N] y_tilde;
  real alpha2 = alpha_mu + beta;
  vector[E] alpha2e = alpha + betae;
  
  for (n in 1:N) {
    log_lik[n] = lognormal_lpdf(y[n] | mu[n], sigma_exp[n]);
    y_tilde[n] = lognormal_rng(mu[n], sigma_exp[n]);
  }
}


