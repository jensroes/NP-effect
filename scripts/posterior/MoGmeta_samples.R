library(magrittr)
library(tidyverse)
library(rstan)

dir("stanout/meta")

m <- readRDS("stanout/meta/MoGmeta.rda")
(params <- names(m)[!grepl("sigma|mu\\[|log_|w\\[|u\\[|_tilde|lp__|eta|tau|delta",names(m))])
traceplot(m, params)

ps <- rstan::extract(m, pars = params)
ps %<>% as_tibble()

write_csv(ps, "stanout/posterior/MoGmeta_posterior.csv")
