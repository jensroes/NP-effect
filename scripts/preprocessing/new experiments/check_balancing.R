library(stringdist)
library(tidyverse)
library(magrittr)

(files <- dir("results"))
files <- files[which(files %in% 60:200)]

length(files)

d <- map2(files, files, ~read_delim(paste0("results/",.,"/RESULTS_FILE.txt"), delim = "\t") %>%
            mutate(subj = .y)) %>%
  bind_rows() 

d %>% count(structure, movement, p1_image) %>% as.data.frame()

d %>% count(type)

d %>% filter(type == "stimuli") %>%
  count(p1_end_y,
        p2_end_y,
        p3_end_y, structure)

d %>% filter(type == "stimuli") %>%
  count(p1_y,
        p2_y,
        p3_y, structure)

d %>% filter(type == "stimuli") %>%
  count(p1_vis,
        p2_vis,
        p3_vis, structure)

d %>% filter(type == "stimuli") %>%
  count(p1_file,
        p2_file,
        p3_file)
