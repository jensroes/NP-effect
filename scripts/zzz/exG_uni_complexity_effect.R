# Ex Gaussian and uniform mixture model [03/08/2018]
# Load packages
source("functions/functions.R")
library(tidyverse)
library(rstan)
library(rethinking)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
n_chain = n_core = 3 # number of cores/chains
iterations = 4e3

# Load data
data <- read_delim("dfs/onsetdata.txt", delim = "\t")  
cond <- factor(data$structure, levels = c("simple", "complex"))
nounphrase <- as.integer(cond) - 1

#data %>% ggplot(aes(x=log(onset), color = structure)) + geom_density()

# Stan model
# Data as list
dat <- 
  within(list(), {
    N <- nrow(data)
    y <- data$onset
    nounphrase <- nounphrase 
    subj <- as.integer(data$exp_subj)
    S <- max(as.integer(data$exp_subj))
    items <- as.integer(data$item)
    I <- max( as.integer(data$item))
  }
  )
str(dat)


# Fit ex Gaussian model
# ----

# Check compiling
exMoG <- stan(file = "stanin/v2/exMoG.stan", data=dat, chains=0)

# Fit model
m = stan(fit = exMoG,
         data = dat,
         iter = iterations,
         warmup= iterations/2,
         chains = n_chain, 
         cores = n_core, 
         refresh = 100,
         seed = 35,
         control = list(max_treedepth = 16,
                        adapt_delta = 0.99))

# Save posterior samples
saveRDS(m,
        file="stanout/exMoG.rda",
        compress="xz")

param <- c( "beta", "tau", "lambda", "sigma")
traceplot(m, param)

as.data.frame(m, pars = param) %>%
  gather(Parameter, value) %>%
  dplyr::group_by(Parameter) %>%
  dplyr::summarise_all(funs( 
    M = dmode(.),
    Lower = HPDI(., prob = .95)[1],
    Upper = HPDI(., prob = .95)[2])
  )


