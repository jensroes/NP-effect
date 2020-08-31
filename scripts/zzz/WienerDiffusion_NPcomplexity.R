# MOG [02/08/2018]
# Load packages
library(bayesplot)
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
nounphrase <- as.numeric(cond)# - 1


library(RWiener)

set.seed(0)
dat <- rwiener(n = 5000, alpha = 1, tau = 0.2, beta = 0.75, delta = 0.5)

#data %>% ggplot(aes(x=log(onset), color = structure)) + geom_density()

# Stan model
# Data as list
dat <- 
  within(list(), {
    N <- nrow(data)
    y <- data$onset
    nounphrase <- nounphrase 
    K <- length(unique(nounphrase))
    subj <- as.integer(data$exp_subj)
    S <- max(as.integer(data$exp_subj))
    items <- as.integer(data$item)
    I <- max( as.integer(data$item))
  }
  )
str(dat)



# Check compiling
#ig <- stan(file = "stanin/v2/inverseGaussian.stan", data=dat, chains=0)

# Fit model
#m = stan(fit = ig, 
#         data = dat,
#         iter = iterations,
#         warmup= iterations/2,
#         chains = n_chain, 
#         cores = n_core, 
#         refresh = 100,
#         seed = 365)

# Save posterior samples
#saveRDS(m,
#        file="stanout/ig.rda",
#        compress="xz")

#param <- c("mu", "lambda")
#traceplot(m, param)



as.data.frame(m, pars = param) %>%
  gather(Parameter, value) %>%
  dplyr::group_by(Parameter) %>%
  dplyr::summarise(M = mean(value),
                   Lower = HPDI(value, prob = .95)[1],
                   Upper = HPDI(value, prob = .95)[2])
