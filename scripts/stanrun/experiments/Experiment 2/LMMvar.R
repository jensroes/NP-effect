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
nounphrase <- as.numeric(cond) - 1

# Stan model
# Data as list
dat <- 
  within(list(), {
    N <- nrow(data)
    y <- data$onset
    nounphrase <- nounphrase
    subj <- as.integer(data$subj)
    S <- max(as.integer(data$subj))
    items <- as.integer(data$item)
    I <- max( as.integer(data$item))
  }
  )
str(dat)

# Initialise start values
start <- function(chain_id = 1){
  list( alpha_mu = 7
        , alpha_sigma = .1
        , alpha_raw = 0
        , beta = 0
        , sigma = 1
        , sigma_diff = .1
        , u = rep(0, dat$S)
        , w = rep(0, dat$I)
        , sigma_u = 0.1
        , sigma_w = 0.1)}

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# Check compiling
lmm <- stan(file = "stanin/experiments/LMMvar.stan", data=dat, chains=0)

# Fit model
m = stan(fit = lmm, 
         data = dat,
         init = start_ll,
         iter = iterations,
         warmup = iterations/2,
         chains = n_chain, 
         cores = n_core, 
         save_warmup = FALSE, # Don't save the warmup
         refresh = 1000,
         thin = 1,
         seed = 81,
         control = list(max_treedepth = 16,
                        adapt_delta = 0.99,
                        stepsize = 2))


# Save posterior samples
saveRDS(m,
        file="stanout/experiments/exp2/LMMvar.rda",
        compress="xz")

param <- c("alpha", "beta", "sigmap_e", "sigma_e")
traceplot(m, param)

