# Load packages
rm(list=ls());gc()
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
n_chain = n_core = 3 # number of cores/chains
iterations = 30000

# Load data
data <- read_csv("dfs/RoeserEtAl2019.csv")  
cond <- factor(data$structure, levels = c("simple", "complex"))
nounphrase <- as.numeric(cond)# - 1
mod <- factor(data$modality, levels = c("spoken", "written"))
modality <- as.numeric(mod)# - 1

# Data as list
dat <- 
  within(list(), {
    N <- nrow(data)
    y <- data$onset
    nounphrase <- nounphrase 
    modality <- modality
    subj <- as.integer(data$exp_subj)
    S <- max(as.integer(data$exp_subj))
    items <- as.integer(data$item)
    I <- max( as.integer(data$item))
    E <- max(data$Exp)
    exps <- data$Exp
  }
  );str(dat)

# Starting values
start <- 
  function(chain_id = 1){
    list(alpha_mu = 7,
         alpha_sigma = .1,
         alpha_raw = rep(0,dat$E),
         delta = .1,
         sigma = .1,
         eta = .1,
         psi = .1,
         sigma_diff = rep(.01, dat$E),
         theta_complex = 0,
         theta_simple = 0,
         theta_simple_writing = 0,
         theta_complex_writing = 0,
         theta_simple_speech = 0,
         theta_complex_speech = 0,
         theta_simple_writing_e = rep(0, dat$E),
         theta_complex_writing_e = rep(0, dat$E),
         theta_simple_speech_e = rep(0, dat$E),
         theta_complex_speech_e = rep(0, dat$E),
         u = rep(0, max(data$exp_subj)),
         w = rep(0, length(unique(dat$item))),
         sigma_u = .1,
         sigma_w = .1
    )
  }
start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# Check compiling
mog <- stan(file = "stanin/modality/MoGmodality.stan", data=dat, chains=0)

# Fit model
m <- stan(fit = mog, 
         data = dat,
         init = start_ll,
         iter = iterations,
         warmup = iterations/2,
         chains = n_chain, 
         cores = n_core, 
         include = FALSE,
         pars = c("mu", "u", "w", "theta_simple", "theta_complex",
                  "theta_simple_writing", "theta_complex_writing",
                  "theta_simple_speech", "theta_complex_speech",
                  "theta_simple_writing_e", "theta_complex_writing_e",
                  "theta_simple_speech_e", "theta_complex_speech_e",
                  "log_theta_simple_writing", "log_theta_complex_writing",
                  "log_theta_simple_speech", "log_theta_complex_speech"), 
         save_warmup = FALSE, # Don't save the warmup
         refresh = 200,
         thin = 1,
         seed = 81,
         control = list(max_treedepth = 16,
                         adapt_delta = 0.99,
                         stepsize = 2))

# Save posterior samples
saveRDS(m,
        file="stanout/modality/MoGmodality.rda",
        compress="xz")

(param <- names(m)[!grepl("log_|tilde|lp__|sigma_u|sigma_w|_e|alpha2", names(m))])
traceplot(m, param)

