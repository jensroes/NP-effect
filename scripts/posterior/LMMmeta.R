library(magrittr)
library(tidyverse)
library(tidybayes)
library(ggeffects)
library(ggthemes)
library(rstan)

# Raw data
data <- read_csv("dfs/PooledData.csv") %>%
  filter(Exp != "Martin et al. (2010, Exp. 4b)") %>% # remove Martin et al. (2010, 4b)
  mutate(ExpID = as.integer(factor(ExpID)),
         subj = as.integer(factor(subj)),
         structure = factor(structure, levels = c("simple", "complex"), ordered = T),
         nounphrase = as.integer(structure) - 1) 

# Posterior samples
d <- read_csv("stanout/posterior/LMMmeta.csv")

d %>% 
  select(starts_with("alpha")) %>%
  mutate_all(~exp(.)) %>%
  pivot_longer(everything(), names_to = "param", values_to = "values") %>%
  separate(param, into = c("param", "study"), sep = "\\[") %>%
  mutate(study = gsub("\\]", "", study),
         study = ifelse(is.na(study), "Overall", study),
         param = gsub("e|_mu","", param)) %>%
  pivot_wider(names_from = param, values_from = values) %>%
  unnest() %>% mutate(beta = alpha2 - alpha) %>% select(-starts_with("a")) -> betas

betas %<>%
  mutate(study = factor(study),
         study = plyr::mapvalues(study, from = c("Overall", 1:8), 
                                 to = c("Overall", levels(factor(data$Exp))) ))
betas %>%
  filter(study == "Overall") %>%
  summarise(M = median(beta)) %>% pull(M) -> overall

betas %>%
  mutate(study = recode(study, Overall = "Overall effect")) %>%
  mutate(group = ifelse(study == "Overall effect", "a", "b")) %>%
  mutate(study = factor(study, levels = c("Overall effect", rev(levels(factor(data$Exp)))), ordered = T)) %>%
  mutate(study = plyr::mapvalues(study, from = c("Overall effect", levels(factor(data$Exp))), to = c("Overall effect", 1:8) )) %>%
  ggplot(aes(x = beta, y = study, fill = group )) +
  geom_vline(xintercept = 0, linetype = "solid", colour = "grey") +
  geom_vline(xintercept = overall, linetype = "dashed", colour = "grey") +
  stat_gradientintervalh(show.legend = F, .width = c(.66, .95)) +
  scale_fill_ptol() +
  theme_minimal() +
  labs(y = "Study", x = bquote("Onset latency"~Delta*hat(beta)~"[in msecs] with 66% and 95% PIs"), 
       subtitle = "NP complexity effect", 
       title = "Meta analysis") +
  theme(axis.title.y = element_text(angle = 360))

ggsave("plots/NPeffectMeta.pdf", width = 7, height = 6, device = cairo_pdf)            
