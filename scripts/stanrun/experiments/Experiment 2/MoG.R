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

cond <- factor(data$structure, levels = c("simple-complex", "complex-simple"))
nounphrase <- as.numeric(cond)# - 1

#data %>% ggplot(aes(x=log(onset), color = structure)) + geom_density()

# Stan model
# Data as list
dat <- 
  within(list(), {
    N <- nrow(data)
    y <- data$onset
    nounphrase <- nounphrase 
    K <- length(unique(nounphrase))
    subj <- as.integer(data$subj)
    S <- max(as.integer(data$subj))
    items <- as.integer(data$item)
    I <- max( as.integer(data$item))
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
         theta_complex = 0,
         theta_simple = 0,
         u = rep(0, max(data$subj)),
         w = rep(0, length(unique(data$item))),
         sigma_u = .1,
         sigma_w = .1
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )


# Check compiling
mog <- stan(file = "stanin/experiments/MoG.stan", data=dat, chains=0)

# Fit model
m = stan(fit = mog, 
         data = dat,
         init = start_ll,
         iter = iterations,
         warmup = iterations/2,
         chains = n_chain, 
         cores = n_core, 
         include = FALSE,
         pars = c("mu", "log_prob_simple", "log_prob_complex",
                  "theta_simple", "theta_complex", "prob_tilde"), 
         save_warmup = FALSE, # Don't save the warmup
         refresh = 250,
         thin = 1,
         seed = 81,
         control = list(max_treedepth = 16,
                         adapt_delta = 0.99,
                         stepsize = 2))


# Save posterior samples
saveRDS(m,
        file="stanout/experiments/exp2/MoG.rda",
        compress="xz")

param <- c("alpha", "delta", "prob_simple", "prob_complex", "sigmap_e", "sigma_e")
traceplot(m, param)
