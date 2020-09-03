library(magrittr)
library(tidyverse)
library(tidybayes)
library(ggeffects)
library(ggthemes)
library(ggstance)
library(stringi)

theme_set(theme_light(base_size = 16) + 
            theme(panel.grid = element_blank(),
                  axis.ticks = element_blank()))

# Posterior samples
exp1 <- readRDS("stanout/experiments/exp1/LMM.rda") 
(pars <- names(exp1)[str_detect(names(exp1), "alpha")])
extract(exp1, pars = pars) %>%
  as_tibble() %>%
  mutate(Experiment = "1") -> exp1_samps

exp2 <- readRDS("stanout/experiments/exp2/LMM.rda")
(pars <- names(exp2)[str_detect(names(exp2), "alpha")])
extract(exp2, pars = pars) %>%
  as_tibble() %>%
  mutate(Experiment = "2") -> exp2_samps

bind_rows(exp1_samps, exp2_samps) %>%
  select(-starts_with("alpha_")) %>%
  pivot_longer(cols = starts_with("alp"), names_to = "NP", values_to = "values") %>%
  mutate(values = exp(values),
         Model = "LMM",
         NP = recode(NP, alpha = "A moved above B and C",
                     alpha2 = "A and B moved above C"),
         NP = factor(NP, levels = unique(NP), ordered = T)) -> lmm

lmm %>% mutate(Experiment = recode(Experiment, "1" = "1\ne.g. 'A and the B moved\nabove the C'",
                             "2" = "2\ne.g. 'A, the B, the C'")) %>%
  ggplot(aes(x = values, fill = NP, colour = NP)) +
  geom_density(alpha = .5) +
  facet_wrap(~Experiment, nrow = 2, strip.position = "left", labeller = label_both)  +
  scale_colour_manual("Stimulus", values = c("black", "firebrick4")) +
  scale_fill_manual("Stimulus", values = c("grey", "firebrick4")) +
  scale_x_continuous(limits = c(800, 1400)) +
  labs(y = "", x = "Onset latency [in msecs]") +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0),
        plot.margin = margin(.1,.1,.1,-.75, "cm"),
        strip.text.y.left = element_text(angle = 360, hjust = 0, vjust = 1, colour =  "black"),
        strip.background = element_rect(fill = "transparent"),
        legend.justification = "right") 

ggsave("presentations/amlap-2020/gfx/NPexp.pdf", width = 8.5, height = 5.25, device = cairo_pdf)


# Difference
bind_rows(exp1_samps, exp2_samps) %>%
  select(-starts_with("alpha_")) %>%
  pivot_longer(cols = starts_with("alp"), names_to = "NP", values_to = "values") %>%
  mutate(values = exp(values)) %>%
  group_by(Experiment) %>%
  pivot_wider(names_from = NP, values_from = values) %>%
  unnest() %>%
  transmute(diff = alpha2 - alpha) %>%
  median_qi() %>%
  mutate_if(is.numeric, round, 0) %>%
  group_by(Experiment) %>%
  transmute(summary = paste0(diff, "ms; ", "PI [", .lower, ", ", .upper, "]")) %>%
  ungroup() -> summary

