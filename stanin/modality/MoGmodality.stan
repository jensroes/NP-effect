/* 
Finite mixture model
This model is adapted from Vasishth et al. 2007 
Random intercepts for subj and items
Mixture on simple and complex condition
Hyper parameter to differentiate between modality
*/

data {
	int<lower=1> N;                    // Number of observations
	real y[N];  		            //outcome
	int<lower=1, upper=2> nounphrase[N];  // 1 simple, 2 complex
	int<lower=1, upper=2> modality[N]; // 1 speech, 2 writing
	
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
  real<lower=0> psi; // variance for experiment level mixing proportion
	vector<lower=0>[E] sigma_diff;

  // Non-centre intercept
  real alpha_mu; // distribution
  real<lower=0> alpha_sigma;
  vector[E] alpha_raw;

  real theta_simple;         // Hyperparameter for mixing proportion (across modality)
  real theta_complex;         // Hyperparameter for mixing proportion (across modality)
  real theta_simple_speech;   // Mixing proportion by modality
  real theta_complex_speech;   // Mixing proportion by modality
  real theta_simple_writing;   // Mixing proportion by modality
  real theta_complex_writing;   // Mixing proportion by modality

  vector[E] theta_simple_writing_e;   // Mixing proportion by modality for each experiment
  vector[E] theta_complex_writing_e;   
  vector[E] theta_simple_speech_e;   
  vector[E] theta_complex_speech_e;   

	// For random effects
	vector[S] u; //subject intercepts
	real<lower=0> sigma_u;//subj sd

	vector[I] w; //items intercepts
	real<lower=0> sigma_w;//items sd
}


transformed parameters{
  vector[E] log_theta_simple_writing_e[2];
  vector[E] log_theta_simple_speech_e[2];
  vector[E] log_theta_complex_writing_e[2];
  vector[E] log_theta_complex_speech_e[2];
  vector[E] prob_simple_writing_e = inv_logit(theta_simple_writing_e); //probability of extreme values by modality
  vector[E] prob_complex_writing_e = inv_logit(theta_complex_writing_e); 
  vector[E] prob_simple_speech_e = inv_logit(theta_simple_speech_e); 
  vector[E] prob_complex_speech_e = inv_logit(theta_complex_speech_e); 
  real prob_simple_writing = inv_logit(theta_simple_writing); //probability of extreme values by modality
  real prob_complex_writing = inv_logit(theta_complex_writing); 
  real prob_simple_speech = inv_logit(theta_simple_speech); 
  real prob_complex_speech = inv_logit(theta_complex_speech); 
  real prob_simple = inv_logit(theta_simple); //probability of extreme values
  real prob_complex = inv_logit(theta_complex); 
	vector<lower=0>[E] sigmap_e = sigma + sigma_diff;
	vector<lower=0>[E] sigma_e = sigma - sigma_diff;
	vector[E] alpha = alpha_mu + alpha_sigma * alpha_raw;

  // inverse logit for prior on theta and log for mixing proprtion
  log_theta_simple_writing_e[1] = log_inv_logit(theta_simple_writing_e);
  log_theta_simple_writing_e[2] = log1m_inv_logit(theta_simple_writing_e);
  log_theta_simple_speech_e[1] = log_inv_logit(theta_simple_speech_e);
  log_theta_simple_speech_e[2] = log1m_inv_logit(theta_simple_speech_e);
  log_theta_complex_writing_e[1] = log_inv_logit(theta_complex_writing_e);
  log_theta_complex_writing_e[2] = log1m_inv_logit(theta_complex_writing_e);
  log_theta_complex_speech_e[1] = log_inv_logit(theta_complex_speech_e);
  log_theta_complex_speech_e[2] = log1m_inv_logit(theta_complex_speech_e);
}

model {
  // Priors
	alpha_mu ~ normal(7, 5);
	alpha_sigma ~ normal(0, 10);
	alpha_raw ~ normal(0, 1);
	
  delta ~ normal(0, 1);
  sigma_diff ~ normal(0,1);
  sigma ~ cauchy(0, 2.5);
  
  // Hyper(hyper)-priors
  eta ~ normal(0, 1); 
  psi ~ normal(0, 1);

  // Experiment-level mixing proportion
  theta_simple_writing_e ~ normal(theta_simple_writing, psi);
  theta_simple_speech_e ~ normal(theta_simple_speech, psi);
  theta_complex_writing_e ~ normal(theta_complex_writing, psi);
  theta_complex_speech_e ~ normal(theta_complex_speech, psi);

  // Modality level mixing proportion
  theta_simple_speech ~ normal(theta_simple, eta);
  theta_simple_writing ~ normal(theta_simple, eta);
  theta_complex_speech ~ normal(theta_complex, eta);
  theta_complex_writing ~ normal(theta_complex, eta);

  // NP level mixing propotion
  theta_simple ~ normal(0, 2);
  theta_complex ~ normal(0, 2);
  
  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //items random effects
  
  // Likelihood	
	for(n in 1:N){
    real mu = alpha[exps[n]] + u[subj[n]] + w[items[n]];
	  if(nounphrase[n]==1){ // simple NPs 
	    if(modality[n]==1){ // speech
        target += log_sum_exp(
    	    log_theta_simple_speech_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta, sigmap_e[exps[n]]), 
    	    log_theta_simple_speech_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
  	  }
 	    if(modality[n]==2){ // writing
    	  target += log_sum_exp(
    	    log_theta_simple_writing_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta, sigmap_e[exps[n]]), 
    	    log_theta_simple_writing_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
  	  }
	   }
     if(nounphrase[n]==2){ // complex NPs 
	    if(modality[n]==1){ // speech
        target += log_sum_exp(
          log_theta_complex_speech_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta, sigmap_e[exps[n]]), 
          log_theta_complex_speech_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
  	  }
 	    if(modality[n]==2){ // writing
        target += log_sum_exp(
          log_theta_complex_writing_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta, sigmap_e[exps[n]]), 
          log_theta_complex_writing_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
  	  }
	  }
	}
}


generated quantities{
  real log_lik[N];
  real y_tilde[N];
  real<lower=0, upper=1> prob_tilde; 
  real prob_diff = prob_complex - prob_simple;
  real prob_diff_speech = prob_complex_speech - prob_simple_speech;
  real prob_diff_writing = prob_complex_writing - prob_simple_writing;
  vector[E] prob_diff_speech_e = prob_complex_speech_e - prob_simple_speech_e;
  vector[E] prob_diff_writing_e = prob_complex_writing_e - prob_simple_writing_e;
  real alpha2 = alpha_mu + delta;
  vector[E] alpha2exp = alpha + delta;

  // likelihood: 
  for(n in 1:N){
    real mu = alpha[exps[n]] + u[subj[n]] + w[items[n]];
    if(nounphrase[n]==1){ // simple NPs 
      if(modality[n]==1){ // speech
        log_lik[n] = log_sum_exp(
  	      log_theta_simple_speech_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta, sigmap_e[exps[n]]), 
  	      log_theta_simple_speech_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
  	    prob_tilde = bernoulli_rng(prob_simple_speech_e[exps[n]]);
  	    if(prob_tilde){
          y_tilde[n] = lognormal_rng(mu + delta, sigmap_e[exps[n]]);
        }
        else{
          y_tilde[n] = lognormal_rng(mu, sigma_e[exps[n]]);
        }
  	  }
  	  if(modality[n]==2){ // writing
        log_lik[n] = log_sum_exp(
  	      log_theta_simple_writing_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta, sigmap_e[exps[n]]), 
  	      log_theta_simple_writing_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
  	    prob_tilde = bernoulli_rng(prob_simple_writing_e[exps[n]]);
  	    if(prob_tilde){
          y_tilde[n] = lognormal_rng(mu + delta, sigmap_e[exps[n]]);
        }
        else{
          y_tilde[n] = lognormal_rng(mu, sigma_e[exps[n]]);
        }	    
	    }
    }
	  if(nounphrase[n]==2){ // complex NPs 
      if(modality[n]==1){ // speech
        log_lik[n] = log_sum_exp(
  	      log_theta_complex_speech_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta, sigmap_e[exps[n]]), 
  	      log_theta_complex_speech_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
  	    prob_tilde = bernoulli_rng(prob_complex_speech_e[exps[n]]);
  	    if(prob_tilde){
          y_tilde[n] = lognormal_rng(mu + delta, sigmap_e[exps[n]]);
        }
        else{
          y_tilde[n] = lognormal_rng(mu, sigma_e[exps[n]]);
        }	    
  	  }
  	  if(modality[n]==2){ // writing
        log_lik[n] = log_sum_exp(
  	      log_theta_complex_writing_e[1, exps[n]] + lognormal_lpdf(y[n] | mu + delta, sigmap_e[exps[n]]), 
  	      log_theta_complex_writing_e[2, exps[n]] + lognormal_lpdf(y[n] | mu, sigma_e[exps[n]]));
  	    prob_tilde = bernoulli_rng(prob_complex_writing_e[exps[n]]);
  	    if(prob_tilde){
          y_tilde[n] = lognormal_rng(mu + delta, sigmap_e[exps[n]]);
        }
        else{
          y_tilde[n] = lognormal_rng(mu, sigma_e[exps[n]]);
        }
	    }
    }
  }
}
