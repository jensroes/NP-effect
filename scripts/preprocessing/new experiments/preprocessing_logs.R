library(stringdist)
library(tidyverse)
library(magrittr)

condnames <-  c("simple-complex", "complex-simple")

(files <- dir("dfs/results/exp1"))
length(files)

d <- map2(files, files, ~read_delim(paste0("dfs/results/exp1/",.,"/RESULTS_FILE.txt"), delim = "\t") %>%
            mutate(subj = .y)) %>%
  bind_rows() %>% 
  filter(structure %in% condnames) %>%
  select(-Session_Name_, -fix_x, -fix_y, -ends_with("vis"),-list,-type)

d[is.na(d$structure),]

logfiles <- map2(files, files, ~read_delim(paste0("dfs/results/exp1/",.x,"/LogFile.txt"), delim = "\t") %>% 
  mutate(subj = .y)) %>%
  bind_rows() %>%
  filter(cond %in% condnames) 


logfiles %<>%
  rename(structure = cond) %>%
  left_join(d, by = c("subj", "trial_id", "structure")) %>%
  select(subj, time, event, key_updown, key, keycode, modified, onset, structure, item, movement, ends_with("image"), trial_id, text) 

logfiles[is.na(logfiles$structure),]

logfiles$structure[logfiles$subj %in% 1:59 & logfiles$structure == "simple-complex" & logfiles$movement == "up"] <- "complex-simple"

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
  mutate(target = ifelse(structure == "simple-complex" & movement == "up",
                         paste0(p1_image, " moved ", "above the ", p2_image, " and the ", p3_image),
                         ifelse(structure == "simple-complex" & movement != "up",
                                paste0(p1_image, " moved ", "below the ", p2_image, " and the ", p3_image),
                                ifelse(structure == "complex-simple" & movement == "up", 
                                       paste0(p1_image, " and the ", p2_image, " moved above the ", p3_image),
                                       ifelse(structure == "complex-simple" & movement != "up", 
                                              paste0(p1_image, " and the ", p2_image, " moved below the ", p3_image), NA))))) %>%
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

# correction that involved going back to the sentence or deleting a space
corr <- c("Delete", "Left")
logfiles %<>% mutate(major_corr = key %in% corr,
                    n_corr = sequence(rle(major_corr)$lengths)-1,
                    n_corr = ifelse(major_corr == FALSE, 0, n_corr)) 
table(logfiles$n_corr)

write_csv(logfiles, "dfs/exp1.csv")

