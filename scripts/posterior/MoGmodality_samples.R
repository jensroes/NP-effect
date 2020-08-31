library(magrittr)
library(tidyverse)
library(rstan)

dir("stanout/modality")

m <- readRDS("stanout/modality/MoGmodality.rda")
(params <- names(m)[!grepl("sigma|mu\\[|log_|w\\[|u\\[|_tilde|lp__|eta|tau|delta",names(m))])
traceplot(m, params)

ps <- rstan::extract(m, pars = params)
ps %<>% as_tibble()

write_csv(ps, "stanout/modality/MoGmodality_posterior.csv")
