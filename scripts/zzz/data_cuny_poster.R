library(tidyverse)
library(gridExtra)

data <- read_delim("dfs/onsetdata.txt", delim = "\t")  
alpha = .5
p <- data %>% 
  mutate(structure = ifelse(structure == "complex", "conjoined NP", "simple NP"),
         experiment = paste0("Experiment ", experiment)) %>%
  ggplot(aes(y=onset/1000, structure)) +
  geom_violin(size = .1, width = .2, alpha =.1) +
  geom_jitter(size = .05, alpha = .2, width = .2) +
  facet_grid(~experiment) +
  theme_linedraw() +
#  scale_y_continuous(trans='log',breaks = seq(100, 10000, 1000)) +
  ggtitle("Raw onset latency data") +
  scale_y_continuous(breaks = seq(0, 15, 1), limits = c(0,10)) +
  labs(x="",
       y="onset latency (in secs)") +
  theme(axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size =9, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 12,face = "bold"));p

cairo_pdf("~/Documents/poster_cuny2019/gfx/data.pdf", width = 7.75, height = 4,  bg = "transparent")
p
dev.off()


data %>% 
  mutate(x = onset,
         structure = ifelse(structure == "complex", "conjoined", structure)) %>%
  ggplot() +
  geom_histogram(aes(x, ..density.., fill = structure), binwidth = 60, color = "white",  alpha = alpha, position = "identity") +
  scale_fill_hue("NP type:", c =30, l =25  ) +
  ylab("") +
  xlab("onset latency (in ms)") +
  #  ggtitle("Sentence production onset") +
  scale_x_continuous(breaks = seq(0, 6000, 1000), limits = c(50, 6000)) + 
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        legend.position = c(.95,.95),#"top",
        legend.justification = c("right", "top"),
        #        legend.box.just = "right",
#        legend.margin=margin(t = 1, unit='cm'),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 9),
        legend.key.height =  unit(1,"line"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text( size = 10),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size = 13, face="bold"),
        text=element_text(family="ArialMT"))

ggsave("slides-ntu-bps/gfx/resultsdata1.pdf", width = 5, height = 4,  bg = "transparent")



