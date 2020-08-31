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
d <- read_csv("stanout/posterior/MoGmeta.csv")

d %>% select(starts_with("prob")) %>%
  pivot_longer(everything(), names_to = "param", values_to = "values") %>%
  filter(!str_detect(param, pattern = "diff")) %>%
  separate(param, into = c("param", "study"), sep = "\\_e\\[") %>%
  mutate(study = gsub("\\]", "", study),
         study = ifelse(is.na(study), "Overall", study)) %>%
  mutate(prob = gsub(pattern = "prob_", "", param)) %>%
  filter(prob != "prob") %>% select(-param) -> thetas

thetas %<>%
  mutate(study = factor(study),
         study = plyr::mapvalues(study, from = c("Overall", 1:8), 
                                 to = c("Overall", levels(factor(data$Exp))) ))
thetas %>%
  filter(study == "Overall") %>%
  group_by(prob) %>%
  summarise(M = median(values)) %>% pull(M) -> overall

thetas %>%
  mutate(study = recode(study, Overall = "Overall effect")) %>%
  mutate(group = ifelse(study == "Overall effect", "a", "b")) %>%
  mutate(study = factor(study, levels = c("Overall effect", rev(levels(factor(data$Exp)))), ordered = T)) %>%
  mutate(study = plyr::mapvalues(study, from = c("Overall effect", levels(factor(data$Exp))), to = c("Overall effect", 1:8) )) %>%
  group_by(study, group, prob) %>%
  summarise(M = median(values),
            lo = quantile(values, .025),
            up = quantile(values, .975)) %>%
  ggplot(aes(y = M, x = study, linetype = prob, shape = prob)) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey") +
  geom_hline(yintercept = overall, linetype = "dashed", colour = "grey") +
  geom_pointrange(aes(ymin = lo, ymax = up), position = position_dodge(.5)) +
  scale_fill_ptol() +
  theme_minimal() +
  scale_shape("Subject NP") +
  scale_linetype("Subject NP") +
  labs(x = "Study", 
       y = bquote("Probability of long values"~hat(theta)~"with 95% PIs"), 
       subtitle = "Probability of long values", 
       title = "Meta analysis") +
  theme(axis.title.y = element_text(angle = 360),
        legend.position = "bottom",
        legend.justification = "right",
        legend.key.height = unit(.75, "cm")) +
  coord_flip()

ggsave("plots/NPMoGMeta.pdf", width = 7, height = 6, device = cairo_pdf)            
