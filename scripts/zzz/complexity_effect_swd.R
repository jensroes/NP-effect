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


# Check compiling
swd <- stan(file = "stanin/v2/swd.stan", data=dat, chains=0)

# Fit model
m = stan(fit = swd, 
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
        file="stanout/swd.rda",
        compress="xz")


(param <- names(m)[1:8])

traceplot(m, param)

library(loo)
log_lik <- extract_log_lik(m)
m.loo <- loo(log_lik)
pareto_k_table(m.loo)
plot(m.loo)

as.data.frame(m, pars = param) %>%
  gather(Parameter, value) %>%
  dplyr::group_by(Parameter) %>%
  dplyr::summarise(M = mean(value),
                   Lower = HPDI(value, prob = .95)[1],
                   Upper = HPDI(value, prob = .95)[2])



#plot
gg_rt<-ggplot(data,aes(x=onset))+
  geom_histogram()+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  geom_text(mapping = aes(x=5000,y=1500,label = paste0("Min. ", min(data$onset))))+
  geom_text(mapping = aes(x=5000,y=1000,label = paste0("Mean. ", round(mean(data$onset),2))))+
  geom_text(mapping = aes(x=5000,y=500,label = paste("SD: ", round(sd(data$onset),3),sep="")))

#plot
gg_rt

gg_res <- stan_hist(m,pars = param)

library(Rmisc)
multiplot(gg_rt,gg_res)
