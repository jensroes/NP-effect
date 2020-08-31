# Load packages
source("functions/get_data.R")
library(tidyverse)
library(magrittr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
n_chain = n_core = 3 # number of cores/chains
iterations = 20000

# Load data
file <- "dfs/exp2.csv"
data <- get_data(file, print_summary = TRUE) %>%
  filter(onset <= 7500)

# Data as list
dat <- 
  within(list(), {
    N <- nrow(data)
    y <- data$onset
    subj <- as.integer(data$subj)
    S <- max(data$subj)
    items <- as.integer(data$item)
    I <- max(data$item)
  }
  )
str(dat)

# Starting values
start <- 
  function(chain_id = 1){
    list(alpha_mu = 6,
         alpha_raw = 0,
         alpha_sigma = .1,
         sigma = .1,
         sigma_diff = .01,
         delta = 0.1,
         theta = 0,
         u = rep(0, max(data$subj)),
         w = rep(0, length(unique(data$item))),
         sigma_u = .1,
         sigma_w = .1
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )


# Check compiling
mog <- stan(file = "stanin/experiments/MoG0.stan", data=dat, chains=0)

# Fit model
m = stan(fit = mog, 
         data = dat,
         init = start_ll,
         iter = iterations,
         warmup = iterations/2,
         chains = n_chain, 
         cores = n_core, 
         include = FALSE,
         pars = c("mu", "log_prob", "theta"), 
         save_warmup = FALSE, # Don't save the warmup
         refresh = 1000,
         thin = 1,
         seed = 81,
         control = list(max_treedepth = 16,
                        adapt_delta = 0.99,
                        stepsize = 2))


# Save posterior samples
saveRDS(m,
        file="stanout/experiments/exp2/MoG0.rda",
        compress="xz")

param <- c("alpha", "delta", "prob", "sigmap_e", "sigma_e")
traceplot(m, param)
