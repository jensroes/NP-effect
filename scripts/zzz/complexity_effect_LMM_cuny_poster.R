library(rethinking)
library(tidyverse)
library(rstan)
source("functions/functions.R")

# Load LMM model
mLMM <- readRDS(file="stanout/lmm.rda")

param <- c("alpha", "beta", "sigma", "sigma_u", "sigma_w")

summary(mLMM, param)$summary
traceplot(mLMM, param)

# Calculate the complexity effect from LMM
as.data.frame(mLMM, pars = c("alpha", "beta", "sigma", "sigma_u", "sigma_w")) %>%
#  mutate(diff = exp(alpha + beta + ((sigma + sigma_u + sigma_w)/2)) - exp(alpha + ((sigma + sigma_u + sigma_w)/2))) %>%
  mutate(diff = exp(alpha + beta) - exp(alpha)) %>%
  dplyr::select(diff) %>%
  dplyr::summarise_all(funs( 
    M = dmode(.),
    Lower = HPDI(., prob = .95)[1],
    Upper = HPDI(., prob = .95)[2])
  ) -> sumsta

as.data.frame(mLMM, pars = c("alpha", "beta")) %>% as.tibble() %>%
  mutate(diff = exp(alpha + beta) - exp(alpha)) %>%
  dplyr::select(diff) %>%
  ggplot(aes(x = diff)) +
  geom_histogram(alpha =.35) +
  labs(x = bquote(hat(beta)~"in ms"), y = " ") +
  theme_minimal() +
  ggtitle(bquote("NP complexity effect"~hat(beta)~"(LMM)"),
          subtitle = "with MAP and 95% HPDI") +
  theme(axis.text.y = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size = 10), 
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size =10)) +
  geom_errorbarh(aes(y =90, xmin=as.vector(sumsta$Lower), xmax =as.vector(sumsta$Upper)),
                 size =.5, height = 30, linetype = "dashed") +
  geom_point(aes(y=90, x= sumsta$M), size =1, shape =21)

ggsave("~/Documents/poster_cuny2019/gfx/lmm.pdf", width = 3.5, height = 3.5,  bg = "transparent")

