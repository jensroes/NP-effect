library(magrittr)
library(tidyverse)
library(rstan)

dir("stanout")

m <- readRDS("stanout/meta/LMMmeta.rda")
(params <- names(m)[!grepl("sigma_exp|mu\\[|log_|w\\[|u\\[|_tilde|lp__|sigma|eta|tau",names(m))])
traceplot(m, params)

ps <- rstan::extract(m, pars = params)
ps %<>% as_tibble()

write_csv(ps, "stanout/posterior/LMMmeta.csv")
