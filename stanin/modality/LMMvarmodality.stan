/*
  lognormal distribution 
  by-subject and by-item intercepts (no random slopes; see Vasishth et al. 2017)
  model for NP complexity effect 
  heteroscedasticity for np conditions
*/
  
data {
  int<lower=1> N;                    // Number of observations
  real<lower=0> y[N];  		            //outcome
  int<lower=0, upper=1> nounphrase[N];  //predictor
	int<lower=1, upper=2> modality[N];

  int<lower=0> E;  // Number of experiments
  int<lower=0, upper=E> exps[N];  // Number of experiments
  
  int<lower=1> S;                  //number of subjects
  int<lower=1, upper=S> subj[N];   //subject id
  
  int<lower=1> I;                  //number of items
  int<lower=1, upper=I> items[N];   //items id
}


parameters {
  real beta; // effect across all experiments
  real<lower=0> tau; // variation between experiments
  vector[E] eta;
  
  real beta_writing; // effect across all experiments by modality
  real beta_speech; 
  real psi;

  real<lower=0> sigma;		// residual sd
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
  vector[E] beta_speech_e = beta_speech + tau * eta; // effect for each experiment
  vector[E] beta_writing_e = beta_writing + tau * eta; // effect for each experiment
  vector[E] alpha = alpha_mu + alpha_sigma * alpha_raw;
  vector<lower=0>[E] sigmap_e = sigma + sigma_diff;
  vector<lower=0>[E] sigma_e = sigma - sigma_diff;
}



model {
  // Priors
  alpha_mu ~ normal(7, 4);
  alpha_sigma ~ normal(0, 5);
  alpha_raw ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5); 
  sigma_diff ~ normal(0,1);
  
	// Hyper priors	
	beta_speech ~ normal(beta, psi);
	beta_writing ~ normal(beta, psi);
	
  // Hyper hyper priors
  beta ~ normal(0, 5);
  tau ~ normal(0, 5);
  eta ~ normal(0, 1);
  psi ~ normal(0, 5);
  
  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects
  
  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //items random effects
  
  // Likelihood
  for(n in 1:N){
    real mu = alpha[exps[n]] + u[subj[n]] + w[items[n]];
    if(nounphrase[n]==0){
      y[n] ~ lognormal(mu, sigma_e[exps[n]]);    
    }
    if(nounphrase[n]==1){
      if(modality[n]==1){// speech
        y[n] ~ lognormal(mu + beta_speech_e[exps[n]], sigmap_e[exps[n]]);    
      }
      if(modality[n]==2){// writing
        y[n] ~ lognormal(mu + beta_writing_e[exps[n]], sigmap_e[exps[n]]);    
      }
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
    if(nounphrase[n]==0){
      log_lik[n] = lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]);
      y_tilde[n] = lognormal_rng(mu, sigma_e[exps[n]]);
    }
    if(nounphrase[n]==1){
      if(modality[n]==1){// speech
        log_lik[n] = lognormal_lpdf(y[n] | mu + beta_speech_e[exps[n]], sigmap_e[exps[n]]);
        y_tilde[n] = lognormal_rng(mu + beta_speech_e[exps[n]], sigmap_e[exps[n]]);
      }
      if(modality[n]==2){// writing
        log_lik[n] = lognormal_lpdf(y[n] | mu + beta_writing_e[exps[n]], sigmap_e[exps[n]]);
        y_tilde[n] = lognormal_rng(mu + beta_writing_e[exps[n]], sigmap_e[exps[n]]);
      }
    }
  }
}
