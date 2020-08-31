# log normal model (standard analysis)
# Fixed effects only

# Load packages
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(bayesplot)
source("functions/functions.R")

# Load data
data <- fread("dfs/onsetdata.txt", head = T, stringsAsFactors = FALSE, sep = "\t")  
cond <- factor(data$structure, levels = c("simple", "complex"))
contrast <- as.numeric(cond) - 1

# Stan model
# Data as list
dat <- 
  within(list(), {
    N <- nrow(data)
    y <- data$onset
    cond <- contrast 
  }
  )
str(dat)

# Check compiling
lognorm <- stan(file = "stanin/lognormal.stan", data=dat, chains=0)

# Fit model
n_chain = n_core = 3 # number of cores/chains
iterations = 2e3
m = stan(fit = lognorm, 
         data = dat,
         iter = iterations,
         warmup= iterations/2,
         chains = n_chain, 
         cores = n_core, 
         refresh = -1,
         seed = 365)

# Save posterior samples
saveRDS(m,
        file="stanout/lognormal.rda",
        compress="xz")

# print summary
print(m,pars=c( "simple", "complex", "simple_sigma", "complex_sigma"),probs=c(0.025,0.975))

# summary
samps <- as.data.frame(m)
(param <- names(samps)[c(1:4)])

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 95% intervals")
y_pred <- as.matrix(m, pars = c("y_tilde"))

y_samp = y_pred[sample(nrow(y_pred), 1000, replace=TRUE),  ]

color_scheme_set("red")
ppc_dens_overlay(y = data$onset, 
                 yrep = y_samp)

# Model checks
samples = get_samples(
  stan_out = m
  , pars_to_keep = c("simple", "complex", "simple_sigma", "complex_sigma")
)


#samples = ggmcmc::ggs(m
#                     , family = "simple"
#                      , inc_warmup = FALSE
#                      , stan_include_auxiliar = FALSE
#)

# traceplot
ggmcmc::ggmcmc(
  D = samples
  ,file = NULL
  ,plot = 'ggs_traceplot'
  ,param_page = 4
)

# full vs partial density histogramm of mu values
ggmcmc::ggmcmc(
  D = samples,
  file = NULL,
  plot = 'ggs_compare_partial'
)

# histogram of mu
ggmcmc::ggmcmc(
  D = samples
  ,file = NULL
  ,plot = 'ggs_histogram'
  ,param_page = 8
)

# auto-correlation
ggmcmc::ggmcmc(
  D = samples,
  file = NULL,
  plot = 'ggs_autocorrelation'
)



