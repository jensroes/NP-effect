library(stringdist)
library(tidyverse)
library(magrittr)

condnames <-  c("simple-complex", "complex-simple")
(files <- dir("dfs/results/exp1"))
#files <- files[!(files %in% 1:59)]
length(files)

d <- map2(files, files, ~read_delim(paste0("dfs/results/exp1/",.,"/RESULTS_FILE.txt"), delim = "\t") %>%
            mutate(subj = .y)) %>%
  bind_rows() %>% 
  filter(structure %in% condnames) %>%
  select(-Session_Name_, -fix_x, -fix_y, -ends_with("vis"),-list,-type)

logfiles <- map2(files, files, ~read_delim(paste0("dfs/results/exp1/",.x,"/LogFile.txt"), delim = "\t") %>% 
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

logfiles %>%
  group_by(trial_id) %>%
  mutate(time = down - min(down) + 1) %>%
  select(event_id, time, trial_id) %>%
  ungroup() %>%
  ggplot(aes(x = time, y = event_id, group = trial_id)) +
  geom_point(size = .5) +
  geom_line(alpha = .2) +
  scale_x_log10() 


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

logfiles$structure[logfiles$subj %in% 1:59 & logfiles$structure == "simple-complex" & logfiles$movement == "up"] <- "complex-simple"

logfiles %<>%
  mutate(target = ifelse(structure == "simple-complex" & movement == "up",
                         paste0(p1_image, " moved ", "above the ", p2_image, " and the ", p3_image),
                         ifelse(structure == "simple-complex" & movement != "up",
                                paste0(p1_image, " moved ", "below the ", p2_image, " and the ", p3_image),
                                ifelse(structure == "complex-simple" & movement == "up", 
                                       paste0(p1_image, " and the ", p2_image, " moved above the ", p3_image),
                                       ifelse(structure == "complex-simple" & movement != "up", 
                                              paste0(p1_image, " and the ", p2_image, " moved below the ", p3_image), NA))))) %>%
  mutate(text = tolower(text),
         target = tolower(target),
         text = trimws(text),
         ld = stringdist(text, target, method = "dl")) %>%
  select(-ends_with("image")) %>%
  group_by(trial_id) %>%
  mutate(concat_text = accumulate(key, ~ paste(.x, .y, sep = ""))) %>%
  mutate(n_rep_keys = sequence(rle(key)$lengths)-1,
         word = word(concat_text, -1),
         word_id = str_count(concat_text, '\\w+'),
         key = gsub(pattern = " ","_", key )) %>%
  ungroup() %>%
  mutate(subj = as.integer(factor(subj)))

logfiles$concat_text

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

write_csv(logfiles, "dfs/exp1.csv")

logfiles %>%
  filter( ld < 5, key_id > 0) %>%
  group_by(structure, key_id) %>%
  summarise(M = median(IKI, na.rm = T),
            se = sd(IKI, na.rm = T)/sqrt(n())) %>%
  as.data.frame() %>%
  ggplot(aes(y = M, x = key_id, color = structure)) +
  geom_pointrange(aes(ymin = M - se, ymax = M + se), fatten = 1) +
  geom_line() +
  scale_y_log10()
