library(stringdist)
library(tidyverse)
library(magrittr)

condnames <-  c("simple-complex", "complex-simple")
(files <- dir("dfs/results/exp2"))
length(files)

d <- map2(files, files, ~read_delim(paste0("dfs/results/exp2/",.,"/RESULTS_FILE.txt"), delim = "\t") %>%
            mutate(subj = .y)) %>%
  bind_rows() %>% 
  filter(structure %in% condnames) %>%
  select(-Session_Name_, -fix_x, -fix_y, -ends_with("vis"),-list,-type)

logfiles <- map2(files, files, ~read_delim(paste0("dfs/results/exp2/",.x,"/LogFile.txt"), delim = "\t") %>% 
  mutate(subj = .y)) %>%
  bind_rows() %>%
  filter(cond %in% condnames) 

logfiles %<>%
  rename(structure = cond) %>%
  left_join(d, by = c("subj", "trial_id", "structure")) %>%
  select(subj, time, onset, event, key, structure, item, movement, ends_with("image"), trial_id, text) %>%
  mutate(trial_id = as.numeric(factor(paste(subj, trial_id))))

logfiles %<>%
  mutate(event = recode(event, "KEYBOARD_firstkey" = "KEYBOARD_down"),
         key = ifelse(is.na(key), " ", key)) %>%
  group_by(trial_id) %>%
  mutate(time = time - min(time),
         time = time + onset) %>% select(-onset) %>%
  arrange(trial_id, key, time, event) %>%
#  select(trial_id, key, time, event) %>%
  group_by(trial_id, key, event) %>%
  mutate(event_id = 1:n()) %>%
  pivot_wider(names_from = event, values_from = time) %>%
  rename(down = KEYBOARD_down,
         up = KEYBOARD_up) %>%
  arrange(trial_id, down) %>%
  group_by(trial_id) %>%
  mutate(event_id = 1:n()) %>%
  ungroup()

ggplot(logfiles, aes(x = down, y = event_id, group = trial_id)) +
  geom_point() +
  geom_line(alpha = .2) 

  
# BACKSPACE
table(logfiles$key)

logfiles %<>%
  mutate(key = ifelse(is.na(key), " ", key)) %>%
  mutate(BACKSPACE = ifelse(key == "Backspace", TRUE, FALSE)) %>%
  group_by(trial_id) %>%
  mutate(key_id = 1:n()) %>%
  mutate(n_bs = sequence(rle(BACKSPACE)$lengths),
         key_id2 = key_id - n_bs
         ) %>%
  group_by(trial_id, BACKSPACE) %>%
  mutate(key_id3 = 1:n()) %>%
  ungroup() %>%
  group_by(trial_id) %>%
  mutate(key_id4 = ifelse(BACKSPACE == TRUE, (key_id2 - (key_id3)), (key_id3)),
         diff = abs(lag(key_id4)-key_id4),
         key_id5 = ifelse(!(is.na(diff) | diff == 1), key_id4 - (diff - 1) , key_id4)) %>%
  ungroup() %>%
  mutate(key_id = key_id5) %>% select(-key_id5, -key_id4,-key_id3, -key_id2, -diff, -n_bs) 


logfiles %<>%
  mutate(target = ifelse(grepl("the", text), paste0(p1_image, " the ",  p2_image, " the ", p3_image),
                         ifelse(!grepl("the", text), paste0(p1_image, " ",  p2_image, " ", p3_image), NA))) %>%
  mutate(text = tolower(text),
         target = tolower(target),
         text = trimws(text),
         ld = stringdist(text, target, method = "dl")) %>%
  mutate(concat_text = accumulate(key, ~ paste(.x, .y, sep = ""))) %>%
  select(-ends_with("image")) %>%
  group_by(trial_id) %>%
  mutate(n_rep_keys = sequence(rle(key)$lengths)-1,
         word = word(concat_text, -1),
         word_id = str_count(concat_text, '\\w+'),
         key = gsub(pattern = " ","_", key )) %>%
  ungroup() %>%
  mutate(subj = as.integer(factor(subj)))

logfiles %<>%
  group_by(trial_id) %>%
  mutate(onset = down[1],
         IKI = down - lag(down),
         IKI = ifelse(key_id == 1, onset, IKI)) 

glimpse(logfiles)

logfiles[is.na(logfiles$subj),]
length(unique(logfiles$subj))

# correction that involved going back to the sentence or deleting a space
mcorr <- c("Delete", "Left")
logfiles %<>% mutate(corr = key %in% mcorr,
                     n_corr = sequence(rle(corr)$lengths)-1,
                     n_corr = ifelse(corr == FALSE, 0, n_corr)) 

write_csv(logfiles, "dfs/exp2.csv")

logfiles %>%
  filter(key_id > 0, key_id < 15, ld < 15, BACKSPACE == FALSE) %>%
  group_by(structure, key_id) %>%
  summarise(IKI = median(IKI)) %>%
  ggplot(aes(y = IKI, x = key_id, color = structure)) +
  geom_point() +
  geom_line() +
  scale_y_log10()
