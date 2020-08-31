library(tidyverse)
library(stringi)
library(xtable)

read_csv("results/meta.csv") %>%
  mutate_if(is.numeric, round, 0) %>%
  unite(col = "elpd_diff", elpd_diff, se_diff, sep = " (") %>%
  unite(col = "elpd_loo", elpd_loo, se_elpd_loo, sep = " (") %>%
  mutate(elpd_diff = paste0(elpd_diff, ")")) %>%
  mutate(elpd_loo = paste0(elpd_loo, ")")) %>%
  mutate(Type = recode(Model, MoGmeta = "MoG-1",
                               MoG0meta = "MoG-0",
                               LMMmeta = "LMM-1",
                               LMM0meta = "LMM-0",
                               LMMvarmeta = "LMM-2")) %>%
  mutate(Description = recode(Model, MoGmeta = "Mixing proportions by NP type",
                       MoG0meta = "Null model",
                       LMMmeta = "Latency difference",
                       LMM0meta = "Null model",
                       LMMvarmeta = "Larger variance for conjoined NP")) %>%
  select(Type, elpd_diff, elpd_loo, Description) %>%
  rename(Model = Type) -> looc;looc

stri_sub(looc$elpd_loo, 5, 1) <- ","
stri_sub(looc[nchar(looc$elpd_diff) >= 10,]$elpd_diff, 3, 1) <- ","

names(looc)[c(2,3)] <- c("$\\Delta\\widehat{elpd}$", "$\\widehat{elpd}$")

print(xtable(looc), include.rownames=FALSE)
