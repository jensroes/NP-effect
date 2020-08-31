/*
lognormal distribution 
by-subject and by-item intercepts (no random slopes; see Vasishth et al. 2017)
model for NP complexity effect as typically assessed in the literature
analysis for each modality
*/

data {
	int<lower=1> N;                    // Number of observations
	real<lower=0> y[N];  		            //outcome
	real<lower=0, upper=1> nounphrase[N];  //predictor
	int<lower=1, upper=2> modality[N];
	
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
  real beta_writing; // effect across all experiments by modality
  real beta_speech; 
  
  real<lower=0> tau; // variation between experiments
  vector[E] eta;

  // Hyper hyper parameter
  real<lower=0> psi; // variance
  real beta; // effect

	// For random effects
	vector[S] u; //subject intercepts
	real<lower=0> sigma_u;//subj sd

	vector[I] w; //items intercepts
	real<lower=0> sigma_w;//items sd
}


transformed parameters{
	vector[E] alpha = alpha_mu + alpha_sigma * alpha_raw;
  vector[E] beta_writing_e = beta_writing + tau * eta;
  vector[E] beta_speech_e = beta_speech + tau * eta;
}


model {
	// Priors
	alpha_mu ~ normal(7, 3);
	alpha_sigma ~ normal(0, 10);
	alpha_raw ~ normal(0, 1);

	sigma ~ cauchy(0, 2.5); 
	
	// Hyper priors	
	beta_speech ~ normal(beta, psi);
	beta_writing ~ normal(beta, psi);
	
  // Hyper hyper priors
  beta ~ normal(0, 10);
  tau ~ cauchy(0, 1);
  eta ~ normal(0, 1);
  psi ~ cauchy(0, 1);

	// REs priors
	sigma_u ~ normal(0,2.5);
	u ~ normal(0, sigma_u); //subj random effects

	sigma_w ~ normal(0,2.5);
	w ~ normal(0, sigma_w); //items random effects

	// Likelihood
  for(n in 1:N){
    real mu = alpha[exps[n]] + u[subj[n]] + w[items[n]];
    if(modality[n]==1){// speech
      y[n] ~ lognormal(mu + beta_speech_e[exps[n]]*nounphrase[n], sigma[exps[n]]);    
    }
    if(modality[n]==2){// writing
      y[n] ~ lognormal(mu + beta_writing_e[exps[n]]*nounphrase[n], sigma[exps[n]]);    
    }
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_tilde;
  vector[E] alpha2_speech_e = alpha + beta_speech_e;
  vector[E] alpha2_writing_e = alpha + beta_writing_e;
  real alpha2_speech = alpha_mu + beta_speech;
  real alpha2_writing = alpha_mu + beta_writing;
  real alpha2 = alpha_mu + beta;
  
  for (n in 1:N) {
    real mu = alpha[exps[n]] + u[subj[n]] + w[items[n]];
    if(modality[n]==1){// speech
      log_lik[n] = lognormal_lpdf(y[n] | mu + beta_speech_e[exps[n]]*nounphrase[n], sigma[exps[n]]);
      y_tilde[n] = lognormal_rng(mu + beta_speech_e[exps[n]]*nounphrase[n], sigma[exps[n]]);
    }
    if(modality[n]==2){// writing
      log_lik[n] = lognormal_lpdf(y[n] | mu + beta_writing_e[exps[n]]*nounphrase[n], sigma[exps[n]]);
      y_tilde[n] = lognormal_rng(mu + beta_writing_e[exps[n]]*nounphrase[n], sigma[exps[n]]);
    }
  }
}
