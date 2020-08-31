library(tidyverse)
library(bayesplot)
library(magrittr)
library(rstanarm)
library(plyr)
library(gridExtra)

# Data
data <- read_delim("dfs/onsetdata.txt", delim = "\t")  %>% 
  mutate(NP = mapvalues(structure, from = unique(structure), to = c("simple NP", "conjoined NP"))) %>%
  mutate(NP = factor(NP, levels = c("simple NP","conjoined NP"), ordered = TRUE))

p <- data %>% ggplot() +
  stat_density(aes(x=onset, group = NP, color = NP), fill = "transparent", geom = "line", size = .25, position = "identity") +
  scale_x_continuous(limits = c(0,9000)) +
#  scale_x_continuous(limits = c(5,9)) +
  scale_color_manual("", values = c("darkred", "forestgreen")) +
  theme_minimal() +
  xlab("onset latency (in ms)") +
  ylab(" ") +
  theme(legend.position = c(.75,.85),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.ticks.y = element_blank());p

pdf("/home/jens/Documents/conferences/abstracts/cuny2019/gfx/density.pdf", width = 3, height = 3, bg = "transparent", family = "Helvetica")
print(p)
dev.off()


p_raw <- data %>% ggplot(aes(y=onset, x=NP)) +
  geom_violin(size=.3) +
  geom_point(position = position_jitter(width = 0.03), alpha = 0.1, size = .1) +
  theme_minimal() +
  ylab("onset latency (in ms)") +
  xlab(" ") +
  theme(axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank());p_raw

p_log <- data %>% ggplot(aes(y=log(onset), x=NP)) +
  geom_violin(size=.3) +
  geom_point(position = position_jitter(width = 0.03), alpha = 0.1, size = .1) +
  theme_minimal() +
  ylab("onset latency (in log ms)") +
  xlab(" ") +
  theme(axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank());p_log

pdf("/home/jens/Documents/conferences/abstracts/cuny2019/gfx/density.pdf", width = 3, height = 4, bg = "transparent", family = "Helvetica")
grid.arrange(p_raw,p_log, ncol =1)
dev.off()


# For overlaying previous plot with models
#m.ct.mog <- readRDS(file="stanout/mog-ct-allcomponents-re.rda")

y_rep <- as.data.frame(m.ct.mog, pars = "y_tilde")
dim(y_rep)
y_samp <- as.data.frame(t(y_rep[sample(1:9000,200),]))
rownames(y_samp) <- NULL

dyrep <- cbind(d, y_samp)
fits <- paste0("Fit",1:200)
colnames(dyrep)[7:206] <- fits
