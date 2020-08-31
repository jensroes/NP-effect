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
  select(starts_with("alpha"), starts_with("delta")) %>% 
  pivot_longer(everything()) %>%
  filter(!grepl("alpha2", name)) %>%
  separate(name, into = c("name", "study"), sep = "\\[") %>%
  mutate(study = gsub("\\]", "", study),
         study = factor(ifelse(is.na(study), "Pooled coefficient", study)),
         name = gsub("_mu","", name),
         study = plyr::mapvalues(study, from = c("Pooled coefficient", 1:8), 
                                 to = c("Pooled coefficient", levels(factor(data$Exp))[c(2,1,3:8)])),
         study = factor(study, levels = c("Pooled coefficient", levels(factor(data$Exp))[c(8:1)]), ordered = TRUE)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  unnest() %>%
  mutate(alpha2 = alpha + delta) %>%
  mutate_if(is.double, exp) %>%
  mutate(delta = alpha2 - alpha) %>%
  pivot_longer(-study) -> deltas;deltas

deltas %>% filter(name != "delta") %>%
  group_by(name, study) %>% median_qi() %>% ungroup() %>%
  ggplot(aes(x = value, y = study, xmin = .lower, xmax = .upper, colour = name, linetype = name)) +
  geom_pointintervalh(position = position_dodgev(height = -.75), alpha = .75,
                      fatten_point = 4.5, interval_size_range = c(0,.9)) +
  scale_colour_manual("", values = c("black", "firebrick4"),
                      breaks=c("alpha", "alpha2"),
                      labels = c(bquote(atop("Average onset", "latency"~alpha)),
                                 bquote(atop(alpha~+~"slowdown for", "long latencies"~delta)))) +
  scale_linetype_manual("", values = c("dashed", "solid"),
                        breaks=c("alpha", "alpha2"),
                        labels = c(bquote(atop("Average onset", "latency"~alpha)),
                                   bquote(atop(alpha~+~"slowdown for", "long latencies"~delta)))) +
  scale_x_continuous(limits = c(750, 2500), breaks = seq(250, 2500, 250)) +
  labs(x = "Onset latency [in msecs]", y = "") +
  theme(legend.position =c(.775,.85),
          legend.margin = margin(-1,.25,0,0, "cm"),
          legend.key.width = unit(1.75, "cm"),
          legend.key.height = unit(1.65, "cm"),
          plot.margin = margin(.1,.6,.1,-.5, "cm"),
          axis.title = element_text(hjust = 0),
          axis.text = element_text(vjust = .25, hjust = .1)) +
  guides(color = guide_legend(override.aes = list(size = 4)))
  
ggsave("presentations/amlap-2020/gfx/NPmetamoglong.pdf", width = 8, height = 5.25, device = cairo_pdf)
  
  
