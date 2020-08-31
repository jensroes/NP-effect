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
modality <- as.numeric(mod)
cond <- factor(data$structure, levels = c("simple", "complex"))
nounphrase <- as.numeric(cond) - 1

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
         beta_writing = 0,
         beta_speech = 0,
         beta = 0,
         psi = .1,
         tau = .1,
         eta = rep(0, dat$E),
         sigma = rep(.01, dat$E),
         u = rep(0, max(data$exp_subj)),
         w = rep(0, length(unique(dat$item))),
         sigma_u = .1,
         sigma_w = .1
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# Check compiling
lmm <- stan(file = "stanin/modality/LMMmodality.stan", data=dat, chains=0)

# Fit model
m <- stan(fit = lmm, 
         data = dat,
         init = start_ll,
         iter = iterations,
         warmup = iterations/2,
         chains = n_chain, 
         cores = n_core, 
         include = FALSE,
         pars = c("mu", "u", "w", "sigma_exp"), 
         save_warmup = FALSE, # Don't save the warmup
         refresh = 2000,
         thin = 1,
         seed = 81,
         control = list(max_treedepth = 16,
                         adapt_delta = 0.99,
                         stepsize = 2))

# Save posterior samples
saveRDS(m,
        file="stanout/modality/LMMmodality.rda",
        compress="xz")

(param <- names(m)[!grepl("log_|tilde|lp__|sigma_u|sigma_w", names(m))])
traceplot(m, param)

