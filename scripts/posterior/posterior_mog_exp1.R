# Load packages
source("functions/get_data.R")
library(tidyverse)
library(magrittr)
library(rstan)
library(ggExtra)

# Load data
#file <- "dfs/exp1.csv"
#data <- get_data(file, print_summary = FALSE)

# Save posterior samples
m <- readRDS(file="stanout/experiments/exp1/MoG.rda")
(param <- names(m)[!grepl("log_|_tilde|w\\[|lp__|u\\[", names(m))])
traceplot(m, param)

(param <- names(m)[grepl("prob", names(m))])

extract(m, param) %>%
  as_tibble() %>%
  select(-prob_tilde, -prob) %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "values") %>%
  mutate(group = ifelse(grepl("_diff", param), "diff", "nondiff")) %>%
  ggplot(aes(x = values, colour = param, fill = param)) +
  geom_density(alpha = .25) +
  facet_grid(~group, scales = "free") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_minimal() +
  theme(strip.text = element_blank())
