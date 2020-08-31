library(brms)

library(bayesplot)
source("functions/functions.R")
library(tidyverse)
library(rethinking)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
n_chain = n_core = 3 # number of cores/chains
iterations = 4e3

# Load data
data <- read_delim("dfs/onsetdata.txt", delim = "\t")  
data$cond <- factor(data$structure, levels = c("simple", "complex"))
data$nounphrase <- as.integer(cond) - 1
data$subj <- as.integer(data$exp_subj)
data$items <- as.integer(data$item)

m0 <- brm(onset ~ 1 + (1|subj) + (1|items), data = data,
         family = lognormal(),
         seed = 1, chains = 3, iter = 3000)

m <- brm(onset ~ nounphrase + (1|subj) + (1|items), data = data,
           family = lognormal(),
           seed = 1, chains = 3, iter = 3000)

m1 <- brm(onset ~ nounphrase + (nounphrase|subj) + (1|items), data = data,
         family = lognormal(),
         seed = 1, chains = 3, iter = 3000)

loo.m0 <- loo(m0, reloo = TRUE)
loo.mlmm <- loo(m, reloo = TRUE)
loo.mlmm1 <- loo(m1, reloo = TRUE)
#loo.m0$elpd_loo
loo.m0$p_loo
loo::compare(loo.m0,loo.mlmm)
#print(loo.m0)
