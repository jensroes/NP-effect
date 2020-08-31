# Compare to conditions of which one consists of a mixture of two normals and one only has one normal
library(data.table)
library(rstan)
library(parallel) # parallelise chains in rstan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data <- fread("dfs/onsetdata.txt", head = T, stringsAsFactors = FALSE, sep = "\t")  
cond <- factor(data$structure, levels = c("simple", "complex"))
#contrast <- (as.numeric(cond) - 1.5)*2
contrast <- as.numeric(cond) - 1

N <- nrow(data)
y <- data$onset

# Stan model
# Data as list
dat <- 
  within(list(), {
    N <- N
    y <- y
    contrast <- contrast # 0 comes from single normal distribution; 1 comes from mixture of K normals
  }
  )
str(dat)

# Check compiling
rstan_options(auto_write = TRUE)
model <- stan(file = "stanin/lognormalbycond.stan", data=dat, chains=0)

n_chain = n_core = 3 # number of cores/chains
iterations = 500
m = stan(fit = model, 
           data = dat,
           iter = iterations,
           warmup= iterations/2,
           chains = n_chain, 
           cores = n_core, 
           refresh = -1,
           seed = 365)

saveRDS(m,
     file="stanout/lognormalbycond.rda",
     compress="xz")

m <- readRDS(file="stanout/lognormalbycond.rda")
print(m,pars=c("alpha","beta", "alpha_sigma", "beta_sigma"),probs=c(0.025,0.975))
# alpha and betas are about identical but the variance of the first mixture component beta[1] is very large
# indv differences? > varying intercepts model!!!
# few extreme values?

# Model checks
samples = ggmcmc::ggs(m
                      , family = "beta"
                      , inc_warmup = FALSE
                      , stan_include_auxiliar = FALSE
)

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

