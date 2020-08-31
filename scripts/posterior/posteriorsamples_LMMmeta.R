library(magrittr)
library(tidyverse)
library(rstan)

path <- "stanout"
dir(path)

m <- readRDS("stanout/LMMmeta.rda")
params <- names(m)[grepl("sigma|beta|alpha|tau|eta",names(m))]
#traceplot(m, keep)

ps <- extract(m, pars = params)
ps %<>% as_tibble()

write_csv(ps, "stanout/LMMmeta_posterior.csv")
