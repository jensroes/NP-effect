source("functions/get_data.R")
library(tidyverse)
library(magrittr)

theme_set(theme_light(base_size = 16) + theme(panel.grid = element_blank(),
                                              axis.ticks = element_blank(),
                                              axis.title = element_text(hjust = 0)))

# Load data
file <- "dfs/exp1.csv"
exp1 <- get_data(file, print_summary = TRUE)
exp1$Experiment <- "1\ne.g. 'A and the B moved\nabove the C'"
file <- "dfs/exp2.csv"
exp2 <- get_data(file, print_summary = TRUE)
exp2$Experiment <- "2\ne.g. 'A, the B, the C'"

data <- bind_rows(exp1, exp2) %>%
  select(onset, Experiment, structure) %>%
  mutate(structure = recode_factor(structure, `simple-complex` = "A moved above\nB and C",
                                              `complex-simple` = "A and B moved\nabove C"),
         Experiment = paste0("Experiment ", Experiment)) %>%
  rename(Stimulus = structure)

ggplot(data, aes(y = onset, x = Experiment, colour = Stimulus)) +
  geom_point(position = position_jitterdodge(jitter.width = .25, dodge.width = .5), 
             size = .1, alpha = .5) +
  geom_boxplot(outlier.shape = NA, size = .25, width = .25, position = position_dodge(.5)) +
  scale_colour_manual(values = c("black", "firebrick4")) +
  scale_y_continuous(limits = c(0, 7500)) +
  labs(y = "Onset latency [in msecs]", x = "") +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(hjust = 0, colour = "black"),
        legend.margin = margin(-1,0,0,0, "cm"),
        legend.position = "bottom",
        legend.justification = "right")
  
ggsave("presentations/amlap-2020/gfx/expsraw.pdf", width = 8, height = 5, device = cairo_pdf)
