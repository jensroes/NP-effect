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
# Posterior samples
exp1 <- readRDS("stanout/experiments/exp1/MoG.rda") 
pars <- c("alpha", "alpha2")
rstan::extract(exp1, pars = pars) %>% as_tibble() %>%
  mutate(Experiment = "1") -> exp1_samps

exp2 <- readRDS("stanout/experiments/exp2/MoG.rda")
rstan::extract(exp2, pars = pars) %>% as_tibble() %>%
  mutate(Experiment = "2") -> exp2_samps

bind_rows(exp1_samps, exp2_samps) %>%
  mutate_if(is.double, exp) %>%
  mutate(delta =  alpha2 - alpha) %>%
  pivot_longer(c(alpha, alpha2, delta), names_to = "Parameter", values_to = "values") %>%
  mutate(group = ifelse(Parameter == "delta", TRUE, FALSE)) -> mog

mog %>%
  group_by(Experiment, Parameter) %>%
  median_qi(values) %>%
  mutate_if(is.numeric, round, 0) %>%
  group_by(Experiment, Parameter) %>%
  transmute(summary = paste0(values, "ms; ", "PI [", .lower, ", ", .upper, "]")) %>%
  ungroup() -> summary; summary


mog %>% mutate(Experiment = recode(Experiment, `1` = "Exp. 1: e.g. 'A and the B moved above the C'",
                                         `2` = "Exp. 2: e.g. 'A, the B, the C'")) %>%
  group_by(Experiment, Parameter) %>% median_qi(values) %>%
  mutate(Parameter = recode_factor(Parameter, alpha = "alpha", delta = "delta", alpha2 = "alpha + delta", .ordered = TRUE)) %>%
  ggplot(aes(x = values, xmin = .lower, xmax = .upper, 
             y = Parameter, colour = Experiment, linetype = Experiment)) +
  geom_pointrangeh(position = position_dodgev(-.75), alpha = .75) +
  scale_y_discrete(labels = c("alpha" = bquote(atop("           Average", "onset latency"~alpha)),
                              "delta" = bquote(atop("       Slowdown for ", "long latencies"~delta)),
                              "alpha + delta" = expression(alpha + delta))) +
  scale_colour_manual("", values = c("black", "firebrick4")) +
  scale_linetype_manual("", values = c("dashed", "solid")) +
  scale_x_continuous(limits = c(250, 2100), breaks = seq(250, 2000, 250)) +
  labs(x = "Onset latency [in msecs]", y = "") +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.justification = "right",
        legend.margin = margin(-.5,.25,0,0, "cm"),
        legend.key.width = unit(1.75, "cm"),
        legend.key.height = unit(1.15, "cm"),
        plot.margin = margin(.1,.6,.1,-.5, "cm"),
        axis.title = element_text(hjust = 0),
        axis.text = element_text(vjust = .25, hjust = .1))

ggsave("presentations/amlap-2020/gfx/NPexpmoglong.pdf", width = 8, height = 5.25, device = cairo_pdf)


