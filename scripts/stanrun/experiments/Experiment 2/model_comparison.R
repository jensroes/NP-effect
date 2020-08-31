library(loo)
library(tidyverse)
library(magrittr)

path <- "stanout/experiments/exp2/"

(files <- dir(path, pattern = ".rda"))

ms <- gsub(files, pattern = ".rda", replacement = "")

for(i in 1:length(files)){
  m <- readRDS(paste0(path, files[i] ))
  log_lik <- extract_log_lik(m, merge_chains = F) 
  r_eff <- relative_eff(exp(log_lik)) 
  assign(paste0("loo_",ms[i]), loo(log_lik, r_eff = r_eff, cores = 2))  
  print(ms[i]); rm(list = "m"); gc()
}

(loos <- ls(pattern = "loo_"))
mc <- do.call(what = loo_compare, args = lapply(loos, as.name))
tibble(name = loos, model = paste0("model",1:length(loos))) -> names

mc %<>% as.data.frame() %>%
  round(2) %>%
  mutate(model=row.names(.)) %>%
  select(model, elpd_diff:se_elpd_loo) %>%
  left_join(names) %>%
  select(-model) %>%
  mutate(name = gsub("loo_", "", name)) %>%
  select(name, elpd_diff:se_elpd_loo) %>%
  rename(Model = name);mc

file_out <- paste0("results/experiments/loos_exp2.csv")
write_csv(mc, file_out)

loo_compare(loo_LMM, loo_LMM0) %>%
  as_tibble() %>%
  filter(elpd_diff != 0) %>%
  mutate(Model = "LMM vs LMM0") %>%
  select(Model, elpd_diff, se_diff) %>%
  mutate(elpd_diff = abs(elpd_diff)) -> lmm

loo_compare(loo_MoG, loo_LMM) %>%
  as_tibble() %>%
  filter(elpd_diff != 0) %>%
  mutate(Model = "MoG vs LMM") %>%
  select(Model, elpd_diff, se_diff)  -> moglmm

models <- bind_rows(lmm, moglmm)  
file_out <- paste0("results/experiments/loos_exp2_diff.csv")
write_csv(models, file_out)




