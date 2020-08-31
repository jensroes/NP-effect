# Ex Gaussian [03/08/2018]
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
exg <- stan(file = "stanin/v2/exgaussian_v3.stan", data=dat, chains=0)

# Fit model
m = stan(fit = exg,
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
        file="stanout/exg_v3.rda",
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


exg <- stan(file = "stanin/v2/exgaussian_v4.stan", data=dat, chains=0)

# Fit model
m = stan(fit = exg,
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
        file="stanout/exg_v4.rda",
        compress="xz")


# with two taus but one beta
exg <- stan(file = "stanin/v2/exgaussian_v2.stan", data=dat, chains=0)

# Fit model
m = stan(fit = exg,
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
        file="stanout/exg_v2.rda",
        compress="xz")


param <- c( "beta", "lambda_simple", "lambda_complex","tau_simple", "tau_complex", "sigma")
traceplot(m, param)

as.data.frame(m, pars = param) %>%
  gather(Parameter, value) %>%
  dplyr::group_by(Parameter) %>%
  dplyr::summarise_all(funs( 
    M = dmode(.),
    Lower = HPDI(., prob = .95)[1],
    Upper = HPDI(., prob = .95)[2])
  )


m2 <- readRDS(file="stanout/exg_v3.rda")
m3 <- readRDS(file="stanout/exg.rda")

library(loo)
log_lik <- extract_log_lik(m)
loo.m <- loo(log_lik)
log_lik <- extract_log_lik(m2)
loo.m2 <- loo(log_lik)
log_lik <- extract_log_lik(m3)
loo.m3 <- loo(log_lik)

loo::compare(loo.m3,loo.m2,loo.m)

plot(loo.m2)
summary(m, param)$summary
#library(gamlss.dist)
#exGAUS()   # 
#simple <- rexGAUS(1000, mu=872, nu=0.00281, sigma=177)
#hist(simple)
#complex<- rexGAUS(1000, mu=872, nu=0.00207, sigma=177)
#hist(complex)

# library(gamlss)
# m1<-gamlss(y~1, family=exGAUS) 
# plot(m1)
#curve(dexGAUS(x, mu=846 ,sigma=199,nu=381), 100, 600, 
#      main = "The ex-GAUS  density mu=300 ,sigma=35,nu=100")
#plot(function(x) pexGAUS(x, mu=300,sigma=35,nu=100), 100, 600, 
#     main = "The ex-GAUS  cdf mu=300, sigma=35, nu=100")


