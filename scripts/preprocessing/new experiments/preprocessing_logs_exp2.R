library(stringdist)
library(tidyverse)
library(magrittr)

condnames <-  c("simple-complex", "complex-simple")

(files <- dir("results/exp2"))
length(files)

d <- map2(files, files, ~read_delim(paste0("results/exp2/",.,"/RESULTS_FILE.txt"), delim = "\t") %>%
            mutate(subj = .y)) %>%
  bind_rows() %>% 
  filter(structure %in% condnames) %>%
  select(-Session_Name_, -fix_x, -fix_y, -ends_with("vis"),-list,-type)

d[is.na(d$structure),]

logfiles <- map2(files, files, ~read_delim(paste0("results/exp2/",.x,"/LogFile.txt"), delim = "\t") %>% 
  mutate(subj = .y)) %>%
  bind_rows() %>%
  filter(cond %in% condnames) 


logfiles %<>%
  rename(structure = cond) %>%
  left_join(d, by = c("subj", "trial_id", "structure")) %>%
  select(subj, time, event, key_updown, key, keycode, modified, onset, structure, item, movement, ends_with("image"), trial_id, text) 

logfiles[is.na(logfiles$structure),]
# !!!!!!!!!!!!!!!!!!!!!!!!

logfiles %<>%
  filter(key_updown == "down") %>%
  mutate(
    trial_id = paste(subj, trial_id, sep ="_")) %>%
  group_by(trial_id) %>%
  mutate(IKI = time - lag(time)) %>%
  ungroup() %>% 
  group_by(trial_id) %>%
  mutate(event_id = 1:n()) %>%
  mutate(IKI = ifelse(event_id == 1 & is.na(IKI), onset, IKI)) %>%
  select(-onset, -time)

logfiles %<>%
  mutate(target = paste(p1_image, p2_image, p3_image, sep = " ")) %>%
  mutate(
    text = tolower(text),
    target = tolower(target),
    text = trimws(text),
    ld = stringdist(text, target, method = "dl")) %>%
  mutate(key = ifelse(is.na(key), " ", key)) %>%
  mutate(space = key == " ",
         event_id = 1:n()) %>%
  mutate(concat_text = accumulate(key, ~ paste(.x, .y, sep = ""))) %>%
  select(-space,-key_updown,-modified, -ends_with("image"))

logfiles %<>% 
  group_by(trial_id) %>%
  mutate(n_rep_keys = sequence(rle(key)$lengths)-1,
         n_keys = n(),
         word_id = str_count(concat_text, '\\w+'),
         key = gsub(pattern = " ","_", key )) %>%
  #  select(ikis_per_total, key, n_word) %>%
  ungroup() %>%
  mutate(subj = as.integer(factor(subj)))


#logfiles %>%
#  select(subj, ld, n_rep_keys, total_backspace, n_keys, total_delete, total_iki, key_per_sec)

logfiles[is.na(logfiles$subj),]
unique(logfiles$subj)

write_csv(logfiles, "data/exp2_alllogfiles.csv")

