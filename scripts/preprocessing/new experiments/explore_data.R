library(dplyr)

data_extra <- read.csv("data/logfiles_extra.csv", stringsAsFactors = F)

## summarize data
sumdata <- data_extra %>%
  filter(!(event %in% c("trial_start", "DISPLAY_RECALL_PREONSET") )) %>%
  group_by(subj, item, trial_id, syn) %>%
  summarize(
    n_keystrokes = sum(grepl("KEYBOARD", event)),
    response = first(response),
    target = first(target),
    edit_distance = adist(tolower(response), tolower(target)),
    median_iki = median(IKI, na.rm = T),
    sd_iki = sd(IKI, na.rm = T),
    n_backspaces = sum(key == "Backspace", na.rm = T),
    n_comma = sum(grep(",", response)),
    iki_before_first_verb = first(IKI[before_verb & startword]),
    iki_before_second_verb = first(IKI[after_verb & startword]),
    iki_after_second_verb = last(IKI[after_verb & endword]),
    iki_before_noun = first(IKI[noun & startword]),
    iki_after_noun = first(IKI[noun & endword])
    )


## Summarize per condition
sumsyn <- sumdata %>%
  group_by(syn) %>%
  summarize(
    mean_edit_distance = mean(edit_distance, na.rm = T),
    sd_edit_distance = sd(edit_distance, na.rm = T),
    mean_backspaces = mean(n_backspaces, na.rm = T),
    sd_backspaces = sd(n_backspaces, na.rm = T),
    mean_iki = mean(median_iki, na.rm = T),
    sd_iki = mean(sd_iki, na.rm = T),
    sum_comma = sum(n_comma),
    n = n(),
    mean_iki_before_first_verb = mean(iki_before_first_verb, na.rm = T),
    sd_iki_before_first_verb = sd(iki_before_first_verb, na.rm = T),
    mean_iki_before_second_verb = mean(iki_before_second_verb, na.rm = T),
    sd_iki_before_second_verb = sd(iki_before_second_verb, na.rm = T),
    mean_iki_after_second_verb= mean(iki_after_second_verb, na.rm = T),
    sd_iki_iki_after_second_verb = sd(iki_after_second_verb, na.rm = T),
    mean_iki_before_noun = mean(iki_before_noun, na.rm = T),
    sd_iki_before_noun = sd(iki_before_noun, na.rm = T),
    mean_iki_after_noun = mean(iki_after_noun, na.rm = T),
    sd_iki_after_noun = sd(iki_after_noun, na.rm = T)
    )


#summarize revisions per student
revisions <- data_extra %>%
  group_by(subj, item, trial_id) %>%
  summarize(max_rev = max(revision_no)) %>%
  group_by(subj) %>%
  summarize(n_revision_trials = sum(max_rev > 0))

  