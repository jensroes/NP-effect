library(rethinking)
library(tidyverse)
library(rstan)
library(gridExtra)
source("functions/functions.R")

# Data
data <- read_delim("dfs/onsetdata.txt", delim = "\t")  
data$cond <- factor(data$structure, levels = c("simple", "complex"))
data$nounphrase <- as.numeric(data$cond)# - 1

# Model
m <- readRDS(file="stanout/mog.rda")

#(param <- names(m)[1:6])

#names(m)[!grepl(pattern = "log_lik", names(m))]
param <- c("NP", "NP2", "NP0" , 
           "sigma_e",  # NP1
           "sigmap_e", # 1 - NP2; 2 - NP0
           "theta_simple", # 1 - NP0; 2 - NP; 3 - NP2
           "theta_complex" # 1 - NP0; 2 - NP; 3 - NP2
           )

summary(m, param)$summary
traceplot(m, param)

samps <- as.data.frame(m, pars = param) %>%
   transmute(NP2 = NP2,
         NP = NP,
         NP0 = NP0,
         psimple = `theta_simple[2]`,
         psimple0 = `theta_simple[1]`,
         psimple2 = `theta_simple[3]`,
         
         pcomplex = `theta_complex[2]`,
         pcomplex0 = `theta_complex[1]`,
         pcomplex2 = `theta_complex[3]`,
         
         sigma = `sigma_e`,
         sigma0 = `sigmap_e[2]`,
         sigma2 = `sigmap_e[1]`
         ) %>%
    gather(Parameter, value) %>%
    group_by(Parameter) %>%
    dplyr::summarise( 
      mu = dmode(value),
      Lower = HPDI(value, prob = .95)[1],
      Upper = HPDI(value, prob = .95)[2]
      )
    

xmax = 5000
binwidth = 60
mu2 <- samps[samps$Parameter == "NP2",]$mu
mu1 <- samps[samps$Parameter == "NP",]$mu
mu0 <- samps[samps$Parameter == "NP0",]$mu

sigma2 <- samps[samps$Parameter == "sigma2",]$mu
sigma1 <- samps[samps$Parameter == "sigma",]$mu
sigma0 <- samps[samps$Parameter == "sigma0",]$mu

psimple2 <- samps[samps$Parameter == "psimple2",]$mu
psimple1 <- samps[samps$Parameter == "psimple",]$mu
psimple0 <- samps[samps$Parameter == "psimple0",]$mu

pcomplex2 <- samps[samps$Parameter == "pcomplex2",]$mu
pcomplex1 <- samps[samps$Parameter == "pcomplex",]$mu
pcomplex0 <- samps[samps$Parameter == "pcomplex0",]$mu

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dlnorm(x, mu, sigma, log= FALSE)
}

psimple <- data %>% 
  filter(nounphrase == 1) %>%
  mutate(x = onset) %>%
  ggplot() +
  geom_histogram(aes(x, ..density.., color = "data"), binwidth = binwidth, 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mu = mu0, sigma = sigma0, lam = psimple0 ),
                aes(colour = "K1"), lwd = .25) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mu = mu1, sigma = sigma1, lam = psimple1 ),
                aes(colour = "K2"), lwd = .25) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mu = mu2, sigma = sigma2, lam = psimple2 ),
                aes(colour = "K3"), lwd = .25) +
  scale_colour_manual("Mixture components:", 
                      values = c("grey", "darkred", "darkblue", "forestgreen"),
                      breaks = c("K1", "K2", "K3", "data"),
                      labels = c(bquote(italic("K")[1]),
                                 bquote(italic("K")[2]),
                                 bquote(italic("K")[3]),
                                 "data") 
                      ) +
  ylab("") +
  xlab("") +
  ggtitle("Simple NP") +
  scale_x_continuous(breaks = seq(0, 6000, 1000), limits = c(0, xmax)) + 
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        legend.position = c(.8,.65),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(color = "white"),
        plot.title = element_text(size = 14),
        text=element_text(family="ArialMT"))

pcomplex <- data %>% 
  filter(nounphrase == 2) %>%
  mutate(x = onset) %>%
  ggplot() +
  geom_histogram(aes(x, ..density.., colour = "data"), binwidth = binwidth,  
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mu = mu0, sigma = sigma0, lam = pcomplex0 ),
                aes(color = "K1"), lwd = .25) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mu = mu1, sigma = sigma1, lam = pcomplex1),
                aes(color = "K2"), lwd = .25) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mu = mu2, sigma = sigma2, lam = pcomplex2),
                aes(color = "K3"), lwd = .25) +
  scale_colour_manual("Mixture components:", 
                      values = c("grey", "darkred", "darkblue", "forestgreen"),
                      breaks = c("K1", "K2", "K3", "data"),
                      labels = c(bquote(italic("K")[1]),
                                 bquote(italic("K")[2]),
                                 bquote(italic("K")[3]),
                                 "data") 
  ) +
  ylab("") +
  xlab("Onset latency (in ms)") +
  ggtitle("Conjoined NP") +
  scale_x_continuous(breaks = seq(0, 6000, 1000), limits = c(0, xmax)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(size =12),
        legend.position = "none",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 14),
        panel.grid.minor = element_blank(),
        text=element_text(family="ArialMT"))

p <- grid.arrange(psimple,pcomplex) 

ggsave("/home/jens/Documents/conferences/abstracts/cuny2019/gfx/results_thom.pdf", p, width = 8, height = 6,  bg = "transparent")

