library(rethinking)
library(tidyverse)
library(rstan)
library(gridExtra)
source("functions/functions.R")

# Data
data <- read_delim("dfs/onsetdata.txt", delim = "\t")  

alpha = .50
binwidth = 70

psimple <- data %>% 
  filter(structure == "simple") %>%
  mutate(x = onset) %>%
  ggplot() +
  geom_histogram(aes(x, ..density.., fill = "data"), binwidth = binwidth, color = "white", alpha = alpha) +
  scale_fill_manual("Mixture\ncomponents:", 
                    breaks = c("data"),
                    values = c(
                      "data" = alpha("grey30", alpha)),
                    labels = c("data") 
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5))) +
  ylab("") +
  xlab("") +
  ggtitle("Simple NP") +
  scale_x_continuous(breaks = seq(0, 6000, 1000), limits = c(0, 6000)) + 
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        legend.position = c(1,1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_text(face = "bold", size = 12, vjust = .5),
        legend.text = element_text(size = 11),
        legend.key.height =  unit(.75,"line"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "white", size = 12),
        plot.title = element_text(size = 13, face="bold"),
        text=element_text(family="ArialMT"));psimple

pcomplex <- data %>% 
  filter(structure == "complex") %>%
  mutate(x = onset) %>%
  ggplot() +
  geom_histogram(aes(x, ..density.., fill = "data"), binwidth = binwidth, color = "white", alpha = alpha) +
  scale_fill_manual("", 
                    breaks = c("data"),
                    values = c(
                      "data" = alpha("grey30", alpha)),
                    labels = c("data") ) +
  ylab("") +
  xlab("onset latency (in ms)") +
  ggtitle("Conjoined NP") +
  scale_x_continuous(breaks = seq(0, 6000, 1000), limits = c(0, 6000)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(size =12),
        legend.position = "none",
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(),
        text=element_text(family="ArialMT"));pcomplex

grid.arrange(psimple,pcomplex) 

ggsave("slides-nut-bps/gfx/resultsdata.pdf", width = 5, height = 4,  bg = "transparent")



# Model
lmm <- readRDS(file="stanout/lmm.rda")

names(lmm)
param <- c("alpha", "beta", "sigma" )

samps <- as.data.frame(lmm, pars = param) %>%
  transmute(conjoined = exp(alpha + beta),
            simple = exp(alpha),
            beta = exp(alpha + beta) -exp(alpha) ) %>% as.tibble();samps

samps %>% 
  ggplot(aes(x=beta)) +
  geom_histogram( color = "white", alpha = alpha, fill = "black", binwidth = 5) +
  ylab("Posterior probability") +
  xlab(bquote(hat(beta)~"(in ms)")) +
  ggtitle(bquote("NP effect ("*hat(beta)==~.(round(dmode(samps$beta)))*"ms,"~"HPDI ["*.(round(HPDI(samps$beta, prob = .95)[1]))*"ms,"~.(round(HPDI(samps$beta, prob = .95)[2]))*"ms])"), 
          subtitle = "conjoined NP - simple NP") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text( size = 12),
        axis.title = element_text(size=12), 
        plot.title = element_text(size = 13, face="bold"),
        plot.subtitle = element_text(size = 10),
        text=element_text(family="ArialMT"))

ggsave("slides-ntu-bps/gfx/resultsdatalmm.pdf", width = 5, height = 5,  bg = "transparent")


# Model
m <- readRDS(file="stanout/mogK3.rda")

param <- c("NP", "NP2", "NP0" , 
           "sigma_e",  # NP1
           "sigmap_e", # 1 - NP2; 2 - NP0
           "theta_simple", # 1 - NP0; 2 - NP; 3 - NP2
           "theta_complex" # 1 - NP0; 2 - NP; 3 - NP2
)


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
  filter(structure == "simple") %>%
  mutate(x = onset) %>%
  ggplot() +
  geom_histogram(aes(x, ..density.., fill = "data"), binwidth = binwidth, color = "white", alpha = alpha) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = mu0, sigma = sigma0, lam = psimple0 ), aes(fill= "K1"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = mu1, sigma = sigma1, lam = psimple1 ), aes(fill= "K2"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = mu2, sigma = sigma2, lam = psimple2 ), aes(fill = "K3"), lwd = .25) +
  scale_fill_manual("Mixture\ncomponents:", 
                    breaks = c("K1", "K2", "K3", "data"),
                    values = c(
                      "data" = alpha("grey70", alpha),
                      "K1" = alpha("darkred", alpha),
                      "K2" = alpha("darkblue", alpha),
                      "K3" = alpha("forestgreen", alpha)),
                    labels = c(bquote(italic("K")[1]),
                               bquote(italic("K")[2]),
                               bquote(italic("K")[3]),
                               "data") 
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.2))) +
  ylab("") +
  xlab("") +
  ggtitle("Simple NP") +
  scale_x_continuous(breaks = seq(0, 6000, 1000), limits = c(0, 6000)) + 
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        legend.position = c(1,1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_text(face = "bold", size = 12, vjust = .5),
        legend.text = element_text(size = 11),
        legend.key.height =  unit(.75,"line"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "white", size = 12),
        plot.title = element_text(size = 13, face="bold"),
        text=element_text(family="ArialMT"));psimple

pcomplex <- data %>% 
  filter(structure == "complex") %>%
  mutate(x = onset) %>%
  ggplot() +
  geom_histogram(aes(x, ..density.., fill = "data"), binwidth = binwidth, color = "white", alpha = alpha) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = mu0, sigma = sigma0, lam = pcomplex0 ), aes(fill= "K1"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = mu1, sigma = sigma1, lam = pcomplex1), aes(fill = "K2"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = mu2, sigma = sigma2, lam = pcomplex2), aes(fill = "K3"), lwd = .25) +
  scale_fill_manual("Mixture\ncomponents:", 
                    breaks = c("K1", "K2", "K3", "data"),
                    values = c(
                      "data" = alpha("grey70", alpha),
                      "K1" = alpha("darkred", alpha),
                      "K2" = alpha("darkblue", alpha),
                      "K3" = alpha("forestgreen", alpha)),
                    labels = c(bquote(italic("K")[1]),
                               bquote(italic("K")[2]),
                               bquote(italic("K")[3]),
                               "data") 
  ) +
  ylab("") +
  xlab("onset latency (in ms)") +
  ggtitle("Conjoined NP") +
  scale_x_continuous(breaks = seq(0, 6000, 1000), limits = c(0, 6000)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(size =12),
        legend.position = "none",
        plot.title = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(),
        text=element_text(family="ArialMT"));pcomplex

p <- grid.arrange(psimple,pcomplex) 

ggsave("slides-ntu-bps/gfx/resultsmog.pdf", p, width = 6, height = 5,  bg = "transparent")


