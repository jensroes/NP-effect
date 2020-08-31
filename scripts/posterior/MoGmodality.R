library(magrittr)
library(tidyverse)
library(tidybayes)
library(ggeffects)
library(ggthemes)
library(rstan)

# Raw data
data <- read_csv("dfs/RoeserEtAl2019.csv")  

# Posterior samples
d <- read_csv("stanout/modality/MoGmodality_posterior.csv")

d %>% select(starts_with("prob")) %>%
  pivot_longer(everything(), names_to = "param", values_to = "values") %>%
  filter(!str_detect(param, pattern = "diff")) %>%
  filter(param != "prob") %>%
  separate(param, into = c("param", "study"), sep = "\\_e\\[") %>%
  separate(param, into = c("remove", "NP", "Modality")) %>%
  mutate(study = gsub("\\]", "", study),
         study = ifelse(is.na(study), "Overall", study),
         Modality = ifelse(is.na(Modality), "Overall", Modality)) %>%
  select(-remove) -> thetas

thetas %>% count(NP, Modality, study)

thetas %>%
  filter(study == "Overall") %>%
  group_by(NP, Modality) %>%
  summarise(M = median(values)) %>% pull(M) -> overall


thetas %>%
  filter(Modality != "Overall") %>%
  group_by(NP, study, Modality) %>%
  summarise(M = median(values),
            lo = quantile(values, .025),
            up = quantile(values, .975)) %>%
  ungroup() %>%
  mutate(study = recode(study, Overall = "Overall effects"),
         study = factor(study, levels = rev(c(1:3, "Overall effects")), ordered = T)) %>%
  ggplot(aes(y = M, x = study, linetype = NP, shape = NP, color = Modality)) +
#  geom_hline(yintercept = overall, linetype = "dashed", colour = "grey") +
  geom_pointrange(aes(ymin = lo, ymax = up), position = position_dodge(.5)) +
  scale_color_wsj() +
  theme_minimal() +
  scale_shape("Subject NP") +
  scale_linetype("Subject NP") +
  scale_y_continuous(limits = c(0,.3)) +
  labs(x = "Study", 
       y = bquote("Probability"~hat(theta)~"with 95% PIs"), 
       subtitle = "Phrase scope probability", 
       title = "Modality differences") +
  theme(axis.title.y = element_text(angle = 360),
        legend.position = "bottom",
        legend.justification = "right",
        legend.key.height = unit(1, "cm")) +
  coord_flip()

ggsave("plots/NPMoGModality.pdf", width = 7, height = 6, device = cairo_pdf)            
