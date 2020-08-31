library(magrittr)
library(tidyverse)
library(ggeffects)
library(ggthemes)
library(ggstance)

theme_set(theme_light(base_size = 16) + theme(panel.grid = element_blank(),
                  axis.ticks = element_blank()))

# Raw data
data <- read_csv("dfs/PooledData.csv") %>%
  filter(Exp != "Martin et al. (2010, Exp. 4b)") %>% # remove Martin et al. (2010, 4b)
  mutate(Exp = gsub("submitted", "2019", Exp),
         subj = as.integer(factor(subj)),
         structure = recode_factor(structure, simple = "simple", complex = "conjoined", .ordered = T),
         Exp = factor(Exp, levels = unique(Exp)[c(5, 6, 4, 3:1, 7, 8)], ordered = T)) %>%
  rename(np = structure) %>%
  select(onset, subj, Exp, np, item); data 
  
data %>% select(Exp, subj) %>% unique() %>%
  group_by(Exp) %>% count() 

ggplot(data, aes(y = onset, x = Exp, colour = np )) +
  geom_point(size = .01, alpha = .45, position = position_jitterdodge(jitter.width = .5, dodge.width = -.75)) +
  geom_boxplot(size = .25, outlier.shape = NA, position = position_dodge(-.75), notch = T, width = .5) +
  scale_colour_manual("NP type", values = c("black", "firebrick4")) +
  scale_y_continuous(limits = c(0, 7500)) +
  labs(x = "", y = "Onset latency [in msecs]" ) +
  theme(legend.position = c(.85,.875),
        axis.title = element_text(hjust = 0)) + coord_flip() 

ggsave("presentations/amlap-2020/gfx/metaraw.pdf", width = 8, height = 5, device = cairo_pdf)
