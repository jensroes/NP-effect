library(tidyverse)
library(gridExtra)

data <- read_delim("dfs/onsetdata.txt", delim = "\t")  

alpha = .5
data %>% 
  ggplot() + 
  geom_histogram(aes(x=onset, ..density..), 
                 binwidth = 70, color = "white", alpha = alpha, fill= "grey5", position = "identity") +
  ylab("") +
  xlab("onset latency (in ms)") +
#  ggtitle(" ",subtitle = 
#            bquote(italic("M") ==~.(round(mean(data$onset)))~italic("SD") ==~.(round((sd(data$onset)))))) +    
  scale_x_continuous(breaks = seq(0, 6000, 1000), limits = c(50, 6000)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size =12),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size=10))

ggsave("slides-ntu-bps/gfx/data.pdf", width = 5, height = 4,  bg = "transparent")



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



