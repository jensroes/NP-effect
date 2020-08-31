# Load packages
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
n_chain = n_core = 3 # number of cores/chains
iterations = 20000

# Load data
data <- read_csv("dfs/PooledData.csv") %>%
  filter(Exp != "Martin et al. (2010, Exp. 4b)") %>% # remove Martin et al. (2010, 4b)
  mutate(ExpID = as.integer(factor(ExpID)),
         subj = as.integer(factor(subj)),
         structure = factor(structure, levels = c("simple", "complex"), ordered = T),
         nounphrase = as.integer(structure)) 

# Data as list
dat <- 
  within(list(), {
    N <- nrow(data)
    y <- data$onset
    nounphrase <- data$nounphrase 
    subj <- data$subj
    S <- max(data$subj)
    items <- data$item
    I <- max(data$item)
    exps <- data$ExpID
    E <- max(data$ExpID)
  }
  );str(dat)

# Starting values
start <- 
  function(chain_id = 1){
    list(
      alpha_mu = 6,
      alpha_sigma = .1,
      alpha_raw = rep(0, dat$E),
      delta = .1,
      theta = 0,
      sigma = 1,
      sigma_diff = rep(.01, dat$E),
      u = rep(0, max(dat$subj)),
      w = rep(0, max(dat$item)),
      sigma_u = .1,
      sigma_w = .1
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# Check compiling
mog <- stan(file = "stanin/meta/MoG0.stan", data=dat, chains=0)

# Fit model
m <- stan(fit = mog, 
         data = dat,
         init = start_ll,
         iter = iterations,
         warmup= iterations/2,
         chains = n_chain, 
         cores = n_core, 
         refresh = 200,
         include = FALSE, # Don't include the following parameters in the output
         pars = c("mu",  "log_theta", "u", "w", "theta"),
         save_warmup = FALSE, # Don't save the warmup
         thin = 1,
         seed = 81,
         control = list(max_treedepth = 16,
                        adapt_delta = 0.99,
                        stepsize = 2));gc()

# Save posterior samples
saveRDS(m,
        file="stanout/meta/MoG0.rda",
        compress="xz")

(param <- names(m)[!grepl("log_|tilde|lp__|alpha2|sigma_e|sigmap", names(m))])
traceplot(m, param)
