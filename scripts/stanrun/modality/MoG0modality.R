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
mod <- factor(data$modality, levels = c("spoken", "written"))
modality <- as.numeric(mod)# - 1

# Data as list
dat <- 
  within(list(), {
    N <- nrow(data)
    y <- data$onset
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
         tau = .1,
         sigma_diff = rep(.01, dat$E),
         theta = 0,
         theta_writing = 0,
         theta_speech = 0,
         theta_writing_e = rep(0, dat$E),
         theta_speech_e = rep(0, dat$E),
         u = rep(0, max(data$exp_subj)),
         w = rep(0, length(unique(dat$item))),
         sigma_u = .1,
         sigma_w = .1
    )
  }
start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# Check compiling
mog <- stan(file = "stanin/modality/MoG0modality.stan", data=dat, chains=0)

# Fit model
m <- stan(fit = mog, 
         data = dat,
         init = start_ll,
         iter = iterations,
         warmup = iterations/2,
         chains = n_chain, 
         cores = n_core, 
         include = FALSE,
         pars = c("mu", "u", "w", "theta_writing","theta_speech", "theta",
                   "theta_writing_e", "theta_speech_e",
                  "log_theta_writing_e", "log_theta_speech_e", "prob_tilde"), 
         save_warmup = FALSE, # Don't save the warmup
         refresh = 200,
         thin = 1,
         seed = 81,
         control = list(max_treedepth = 16,
                         adapt_delta = 0.99,
                         stepsize = 2))

# Save posterior samples
saveRDS(m,
        file="stanout/modality/MoG0modality.rda",
        compress="xz")

(param <- names(m)[!grepl("log_|tilde|lp__|sigma_u|sigma_w|_e|alpha2", names(m))])
traceplot(m, param)

