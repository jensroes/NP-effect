library(tidyverse)
library(stringi)
library(xtable)

dir("results")

exp1 <- read_csv("results/loos_exp1.csv") %>%
  mutate(exp = "Experiment 1")
  
exp2 <- read_csv("results/loos_exp2.csv") %>%
  mutate(exp = "Experiment 2")


bind_rows(exp1, exp2) %>%
  mutate_if(is.numeric, round, 0) %>%
  unite(col = "elpd_diff", elpd_diff, se_diff, sep = " (") %>%
  unite(col = "elpd_loo", elpd_loo, se_elpd_loo, sep = " (") %>%
  mutate(elpd_diff = paste0(elpd_diff, ")")) %>%
  mutate(elpd_loo = paste0(elpd_loo, ")")) %>%
  mutate(Model = recode(Model, MoG = "MoG-1",
                               MoG0 = "MoG-0",
                               LMM = "LMM-1",
                               LMM0 = "LMM-0",
                               LMMvar = "LMM-2")) %>% 
  select(exp, everything()) %>%
  mutate(exp = c("Expeirment 1", rep("", 4), "Expeirment 2", rep("", 4))) -> looc;looc

stri_sub(looc$elpd_loo, 4, 1) <- ","

#looc %>% pivot_wider(names_from = exp, values_from = c(elpd_diff, elpd_loo))

names(looc)[c(3,4)] <- c("$\\Delta\\widehat{elpd}$", "$\\widehat{elpd}$")


print(xtable(looc), include.rownames=FALSE)
