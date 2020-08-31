library(magrittr)
library(tidyverse)
library(tidybayes)
library(ggeffects)
library(ggthemes)
library(rstan)

d <- read_csv("stanout/posterior/LMMmeta.csv")

d %>% mutate_at(vars(matches("beta")), list(~alpha + .)) %>%
  select(alpha, starts_with("beta")) %>%
  mutate_all(~exp(.)) %>%
  transmute_at(vars(matches("beta")), ~.-alpha) -> betas

data <- read_csv("dfs/PooledData.csv") %>%
  filter(Exp != "Martin et al. (2010, Exp. 4b)") %>% # remove Martin et al. (2010, 4b)
  mutate(ExpID = as.integer(factor(ExpID)),
         subj = as.integer(factor(subj)),
         structure = factor(structure, levels = c("simple", "complex"), ordered = T),
         nounphrase = as.integer(structure) - 1) 

names(betas)[-1] <- levels(factor(data$Exp))

betas %>%
  gather(Param, value) %>%
  group_by(Param) %>% summarise(M = mean(value)) %>% 
  filter(Param == "beta") %>% pull(M) -> overall

betas %>%
  gather(Param, value) %>%
  mutate(Param = recode(Param, beta = "Overall effect")) %>%
  mutate(group = ifelse(Param == "Overall effect", "a", "b")) %>%
  mutate(Param = factor(Param, levels = c("Overall effect", rev(levels(factor(data$Exp)))), ordered = T)) %>%
  ggplot(aes(x = value, y = Param, fill = group )) +
  geom_vline(xintercept = 0, linetype = "solid", colour = "grey") +
  geom_vline(xintercept = overall, linetype = "dashed", colour = "grey") +
  stat_gradientintervalh(show.legend = F) +
  scale_fill_ptol() +
  theme_minimal() +
  labs(y = "", x = bquote("Onset latency"~Delta*hat(beta)~"[in msecs] with 66% and 95% PIs"), 
       title = "NP complexity effect", 
       subtitle = "Meta analysis") +
  theme(axis.title.y = element_text(angle = 360))

ggsave("plots/NPeffectMeta.pdf", width = 8, height = 6, device = cairo_pdf)            
