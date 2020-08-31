library(magrittr)
library(tidyverse)
library(tidybayes)
library(ggeffects)
library(ggthemes)
library(ggstance)

theme_set(theme_light(base_size = 16) + 
            theme(panel.grid = element_blank(),
                  axis.ticks = element_blank()))
        
# Raw data
data <- read_csv("dfs/PooledData.csv") %>%
  filter(Exp != "Martin et al. (2010, Exp. 4b)") %>% # remove Martin et al. (2010, 4b)
  mutate(Exp = gsub("submitted", "2019", Exp)) 

# Posterior samples
read_csv("stanout/posterior/MoGmetadelta.csv") %>% 
  select(starts_with("prob")) %>% select(-starts_with("prob_diff")) %>%
  pivot_longer(everything(), names_to = "param", values_to = "values") %>%
  separate(param, into = c("param", "study"), sep = "\\[") %>%
  mutate(study = gsub("\\]", "", study),
         study = factor(ifelse(is.na(study), "Pooled coefficient", study)),
         param = gsub("_e|prob_","", param),
         param = recode_factor(param, simple = "simple", complex = "conjoined", .ordered = TRUE),
         study = plyr::mapvalues(study, from = c("Pooled coefficient", 1:8), 
                                 to = c("Pooled coefficient", levels(factor(data$Exp))[c(2,1,3:8)])),
         study = factor(study, levels = c("Pooled coefficient", levels(factor(data$Exp))[c(8:1)]), ordered = TRUE)) %>%
  rename(NP = param) -> thetas;thetas
  
thetas %>% group_by(NP, study) %>% median_qi() %>%
  ggplot(aes(x = values, y = study, xmin = .lower, xmax = .upper, colour = NP, linetype = NP)) +
  geom_pointintervalh(position = position_dodgev(height = -.75), 
                      fatten_point = 4.5, interval_size_range = c(0,.9)) +
  scale_colour_manual("NP type", values = c("black", "firebrick4")) +
  scale_linetype_manual("NP type", values = c("dotted", "solid")) +
  scale_x_continuous(limits = c(0, .45)) +
  labs(y = "", x = "Probability of long latencies", linetype = "NP type" ) +
  theme(legend.position = c(.8,.875),
        legend.key.width = unit(1.5, "cm"),
        plot.margin = margin(.1,.1,.1,-.5, "cm"),
        axis.title = element_text(hjust = 0)) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave("presentations/amlap-2020/gfx/mogmeta.pdf", width = 8, height = 5, device = cairo_pdf)


thetas %>% group_by(study) %>%
  pivot_wider(names_from = NP, values_from = values) %>%
  unnest() %>% transmute(values = conjoined - simple) %>% median_qi() %>%
  ggplot(aes(x = values, y = study, xmin = .lower, xmax = .upper)) +
  geom_pointintervalh(fatten_point = 4.5, interval_size_range = c(0,.9)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey") +
  scale_x_continuous(limits = c(-.2, .3)) +
  labs(y = "", x = "Probability of long latencies") +
  theme(axis.title = element_text(hjust = 0),
        plot.margin = margin(.1,.1,.1,-.5, "cm")) 
