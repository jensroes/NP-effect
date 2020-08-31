library(tidyverse)



d %>%
  mutate(first_key = ifelse(lag(key)=="_" | event_id == 1, TRUE, FALSE)) %>%
  select(word_id, event_id, key, trial_id, IKI, structure, first_key, concat_text) %>%
  mutate(word = word(concat_text, -1))  %>%
  mutate(word_id = ifelse(first_key == TRUE & word_id ==1, 0, word_id)) %>%
  group_by(trial_id, word_id, first_key, structure) %>%
  summarise(M = median(IKI))  %>%
  arrange(trial_id, word_id, desc(first_key)) %>%
  group_by(trial_id) %>%
  mutate(word_id = 1:n(),
         N = n()) %>%
  ggplot(aes(y = M, x = word_id, color = structure)) +
  geom_point() +
  scale_y_log10() +
  facet_wrap(~N) +
  geom_smooth()


d %>%
  group_by(trial_id) %>% 
  mutate(max_corr = max(n_corr),
         max_pause = ifelse(event == "KEYBOARD_down", max(IKI), NA),
         max_pause = ifelse(is.na(max_pause), lead(max_pause, 1), max_pause)) %>%
  ungroup() %>%
  rename(onset = IKI) %>%
  filter(event == "KEYBOARD_firstkey",
         max_corr < 10,
         max_pause < 5000,
         ld < 11,
         onset > 50) %>%
  select(subj, structure, item, movement, onset) -> data

ggplot(data, aes(y = onset, x = structure)) +
  geom_point(size = .5, alpha = .5, position = position_jitter(width = .25)) +
  geom_boxplot(outlier.shape = NA, width = .5) +
  scale_y_log10()

