library(magrittr)
library(tidyverse)
library(stringi)


exp1 <- read_csv("results/loos_exp1_diff.csv") %>%
  mutate(Experiment = "1")

exp2 <- read_csv("results/loos_exp2_diff.csv") %>%
  mutate(Experiment = "2")

bind_rows(exp1, exp2) %>%
  mutate_if(is.numeric, round, 0) %>%
  #  filter(elpd_diff != 0) %>%
  unite(col = "elpd_diff", elpd_diff, se_diff, sep = " (") %>%
  unite(col = "elpd_loo", elpd_loo, se_elpd_loo, sep = " (") %>%
  mutate(elpd_diff = paste0(elpd_diff, ")"),
         elpd_loo = paste0(elpd_loo, ")"),
         Comparison = recode(Comparison, "LMM vs LMM0" = "LMM-1 -- LMM-0",
                             "MOG vs LMM" = "MoG-1 -- LMM-1"),
         Model = recode(Model, LMM0 = "LMM-0",
                        LMM = "LMM-1",
                        MoG = "MoG-1"),
         elpd_diff = ifelse(Comparison == "LMM-1 -- LMM-0" & 
                              Model == "LMM-1", gsub("[-]", "", elpd_diff), elpd_diff)) -> looc;looc

stri_sub(looc$elpd_loo, 4, 1) <- ","


model <- looc %>% #filter(elpd_diff == "0 (0)") %>%
  filter(Model %in% c("LMM-1", "MoG-1")) %>%
  mutate(elpd_loo = paste0(Model, ": $\\widehat{elpd}$=", elpd_loo )) %>%
  select(Model, Experiment, elpd_loo) %>% unique()

loos <- mutate(looc, result = paste0(Comparison, ": ", elpd_diff )) %>%
  filter(elpd_diff != "0 (0)") %>%
  mutate(elpd_diff = paste0(Comparison, ": $\\Delta\\widehat(elpd)$=", elpd_diff)) %>%
  select(Comparison, Experiment, elpd_diff) %>%
  separate(Comparison, into = c("Model", "remove"), sep = " -- ", remove = FALSE) %>%
  full_join(model) %>%
  select(-remove) 

# For slide
result <- loos %>% select(Experiment, elpd_diff, elpd_loo) ;result
result %>% pull(elpd_diff)
result %>% pull(elpd_loo)
