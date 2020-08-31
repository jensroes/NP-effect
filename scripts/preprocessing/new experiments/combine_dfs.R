library(tidyverse)
library(magrittr)
library(stringdist)

(files <- dir("results"))
length(files)

d <- map(files, ~read_delim(paste0("results/",.,"/RESULTS_FILE.txt"), delim = "\t")) %>% 
  bind_rows() %>% 
  filter(structure %in% c("simple-complex", "complex-simple"))

d %>% select(structure, movement, text)

d %<>% mutate(structure = ifelse(subj %in% 1:59 & structure == "simple-complex" & movement == "up", "complex-simple", structure))

d %<>%
  mutate(subj = as.integer(factor(paste(Session_Name_, subj)))) %>%
  select(-Session_Name_, -trial_id, -fix_x,-fix_y, -ends_with("vis"),-list,-type)

d %<>%
  mutate(target = ifelse(structure == "simple-complex" & movement == "up",
                         paste0(p1_image, " moved ", "above the ", p2_image, " and the ", p3_image),
                         ifelse(structure == "simple-complex" & movement != "up",
                                paste0(p1_image, " moved ", "below the ", p2_image, " and the ", p3_image),
                                ifelse(structure == "complex-simple" & movement == "up", 
                                       paste0(p1_image, " and the", p2_image, " moved above the ", p3_image),
                                       ifelse(structure == "complex-simple" & movement != "up", 
                                              paste0(p1_image, " and the ", p2_image, " moved below the ", p3_image), NA)))))# %>% select(target)
d %<>% select(-starts_with("p"))

d %<>%
  mutate(
    text = tolower(text),
    target = tolower(target),
    text = trimws(text),
    ld = stringdist(text, target, method = "dl"))

d %>% select(ld, text, target) %>%
  filter(ld > 9)

d %>% 
  filter(ld < 9, onset < 20000) %>%
  group_by(structure) %>%
  summarise(M = mean(onset),
            SD = sd(onset),
            N = n())

d %>%
  filter(ld < 10) %>%
  filter(onset < 20000) %>%
  ggplot(aes(y = onset, x = structure)) +
  geom_boxplot()

d %>% filter(onset < 1)

m <- lme4::lmer(log(onset)~structure*movement + (1|subj), d, subset = ld < 9 & onset < 15000)
summary(m) %>% coef() %>% round(2)

