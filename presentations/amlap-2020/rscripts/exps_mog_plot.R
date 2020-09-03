library(magrittr)
library(tidyverse)
library(tidybayes)
library(ggeffects)
library(ggthemes)
library(ggstance)

theme_set(theme_light(base_size = 16) + 
            theme(panel.grid = element_blank(),
                  axis.ticks = element_blank()))

# Posterior samples
exp1 <- readRDS("stanout/experiments/exp1/MoG.rda") 
(pars <- names(exp1)[str_detect(names(exp1), "prob")])
rstan::extract(exp1, pars = pars) %>%
  as_tibble() %>%
  mutate(Experiment = 1) -> exp1_samps

exp2 <- readRDS("stanout/experiments/exp2/MoG.rda")
(pars <- names(exp2)[str_detect(names(exp2), "prob")])
rstan::extract(exp2, pars = pars) %>%
  as_tibble() %>%
  mutate(Experiment = 2) -> exp2_samps

bind_rows(exp1_samps, exp2_samps) %>% select(-prob_diff) %>%
  pivot_longer(starts_with("prob"), names_to = "NP", values_to = "values") %>%
  separate(NP, into = c("prob", "NP")) %>% select(-prob) %>%
  mutate(NP = recode(NP, simple = "A moved above B and C",
                     complex = "A and B moved above C"),
         NP = factor(NP, levels = unique(NP), ordered = T)) -> mog

mog %>% mutate(Experiment = recode(Experiment, "1" = "1\ne.g. 'A and the B moved\nabove the C'",
                             "2" = "2\ne.g. 'A, the B, the C'")) %>%
  ggplot(aes(x = values, fill = NP, colour = NP)) +
  geom_density(alpha = .5) +
  facet_wrap(~Experiment, nrow = 2, strip.position = "left", labeller = label_both)  +
  scale_colour_manual("Stimulus", values = c("black", "firebrick4")) +
  scale_fill_manual("Stimulus", values = c("grey", "firebrick4")) +
  scale_x_continuous(limits = c(0, .4)) +
  scale_y_continuous(limits = c(0, 20)) +
  labs(y = "", x = "Probability of long onset latencies") +
  theme(legend.position = "bottom",
        axis.title = element_text(hjust = 0),
        axis.text.y = element_blank(),
        plot.margin = margin(.1,.1,.1,-.75, "cm"),
        strip.text.y.left = element_text(angle = 360, hjust = 0, vjust = 1, colour =  "black"),
        strip.background = element_rect(fill = "transparent"),
        legend.justification = "right")

ggsave("presentations/amlap-2020/gfx/NPexpmog.pdf", width = 8.5, height = 5.25, device = cairo_pdf)


bind_rows(exp1_samps, exp2_samps) %>%
  group_by(Experiment) %>%
  median_qi(prob_diff) %>%
  mutate_if(is.numeric, round, 2) %>%
  group_by(Experiment) %>%
  transmute(summary = paste0(prob_diff, " PI [", .lower, ", ", .upper, "]")) %>%
  ungroup() -> summary;summary
