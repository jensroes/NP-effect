library(stringdist)
library(tidyverse)
library(magrittr)

data <- read_csv( "dfs/exp1.csv") 
total <- length(unique(data$subj)) * length(unique(data$item))

self_corr <- unique(data[data$key_id %in% c(-2,-1,0),]$trial_id); 

data %>% mutate(first_word = word(text, 1)) %>%
  select(trial_id, first_word, text) %>%
  filter(!(first_word %in% c("peter", "pater", "petter", "tanua", "oeter",
                             "tania", "tana", "tnaia", "tainia", "tani", "taia", "tanie",
                             "taina", "tanian", "talia", "tanya", "tanai", "tanina",
                             "tanis", "peterand", "taniaand"))) %>%
  unique() %>% pull(trial_id) -> incorrect_start

data %>% filter( ld >= 10) %>% 
  select(trial_id, text, ld) %>% 
  unique() %>% select(trial_id) -> large_ld

data %>% filter(key_id == 1, BACKSPACE == FALSE, onset != 0) %>%
  filter(!(trial_id %in% c(self_corr, incorrect_start, large_ld))) -> d2

d2 %>% count(subj) %>% arrange(n) %>% filter(n < 10) %>% pull(subj) -> remove_subj

d2 %>% filter(!(subj %in% remove_subj)) %>%
  mutate(subj = as.integer(factor(subj)),
         item = as.integer(factor(item))) %>%
  select(onset, subj, item, movement, structure, text) -> dfin

print(paste0("Self corrections: ", 
             length(self_corr), " trials (", round(length(self_corr)/total*100,2) ,"%)"))

print(paste0("Incorrect first word: ", 
             length(incorrect_start), " trials (", round(length(incorrect_start)/total*100,2) ,"%)"))

print(paste0("Final text Levenshtein distance >= 10: ", 
             length(large_ld), " trial(s) (", round(length(large_ld)/total*100,2) ,"%)"))

print(paste0("Participants that contributed < 10 trials: ", length(remove_subj)))

print(paste0("Remaining participants: ", length(unique(dfin$subj))))

return(dfin)

dfin %>% group_by(structure) %>% summarise(N = n())
summary(dfin$onset)
m <- lme4::lmer(log(onset) ~ structure + (1|subj), dfin)
summary(m)
