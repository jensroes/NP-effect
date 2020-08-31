library(stringdist)
library(tidyverse)
library(magrittr)

d <- read_csv( "data/alllogfiles.csv")
glimpse(d)

d[is.na(d$structure),]$subj

d %<>%
  group_by(trial_id) %>%
  mutate(total_backspace = length(which(key == "Backspace")),
         total_left = length(which(key == "Left")),
         total_delete = length(which(key == "Delete")),
         total_edit = total_backspace + total_delete,
         total_iki = sum(IKI[-1]),
         keys_per_sec = n_keys/total_iki*1000) %>%
  #  select(ikis_per_total, key, n_word) %>%
  ungroup() %>%
  mutate(subj = as.integer(factor(subj)))

unique(d$subj)

#d %<>% filter(!(subj %in% 1:59))

d %>% select(subj, ld, key, n_rep_keys, total_backspace,  total_iki, keys_per_sec, word_id)

d %>%
  mutate(subj = factor(subj)) %>%
  #mutate(subj_id = fct_reorder(subj, keys_per_sec, .fun='median')) %>%
  ggplot(aes(y = keys_per_sec, x=reorder(subj, keys_per_sec), fill=subj)) +
  geom_violin(show.legend = F)

d %>%
  group_by(subj) %>%
  summarise(keys_per_sec = median(keys_per_sec),
            total_edit = median(total_edit)) -> d.subj

# Calculating density: d
dens <- density(d.subj$keys_per_sec)

# Use which.max() to calculate mode
(mode <- dens$x[which.max(dens$y)])

# Finish the ggplot call
ggplot(d.subj, aes(x = keys_per_sec)) +
  geom_rug() +
  geom_density() +
  geom_vline(xintercept = mode, col = "red")


# Calculating density: d
dens_edit <- density(d.subj$total_edit)

# Use which.max() to calculate mode
(mode_edit <- dens_edit$x[which.max(dens_edit$y)])

# Finish the ggplot call
ggplot(d.subj, aes(x = total_edit)) +
  geom_rug() +
  geom_density() +
  geom_vline(xintercept = mode_edit, col = "red")


# Finish the ggplot call
ggplot(d.subj, aes(x = keys_per_sec)) +
  geom_rug() +
  geom_density() +
  geom_vline(xintercept = mode, col = "red")


d %>%
  group_by(subj) %>%
  mutate(m_keys_per_sec = median(keys_per_sec),
         m_total_edit = median(total_edit)) %>%
  mutate(slow_typing = m_keys_per_sec < mode,
         many_edits = m_total_edit > mode_edit) %>%
  mutate(struggling = many_edits == slow_typing) -> d2

d[is.na(d$structure),]

d2 %>%
  ungroup() %>%
  filter(event == "KEYBOARD_firstkey" & !is.na(structure), ld < 9) %>%
  filter(IKI > 0,
         IKI < 15000,
         ld < 10) %>%
  mutate(struggling = factor(ifelse(struggling == TRUE, "struggling writers", "not struggling writers")),
         structure = factor(structure)) -> d3

contrasts(d3$structure) <- contr.sum(2)
contrasts(d3$struggling) <- contr.sum(2)

m <- lme4::lmer(log(IKI) ~ structure + (1|subj) + (1|item), d3)
m %>% summary() %>% coef() %>%
  round(2)


m <- lme4::lmer(log(IKI) ~ struggling*structure + (1|subj) + (1|item), d3)
m %>% summary() %>% coef() %>%
  round(2)

d3$COND <- paste(d3$structure, d3$struggling, sep = "_")
m <- lme4::lmer(log(IKI) ~ 0 + COND + (1|subj) + (1|item), d3)
m %>% summary() %>% coef() %>%
  round(2) %>%
  as.data.frame() %>%  
  rownames_to_column("COND") %>%
  mutate(COND = gsub("COND", "", COND)) %>%
  separate(COND, into = c("structure", "writer"), sep = "_") %>%
  ggplot(aes(y = Estimate, x = structure, color = writer, ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`)) +
  geom_pointrange()




d2 %>% select(trial_id, struggling, event, IKI, structure, event_id) %>%
  group_by(event, struggling, subj, structure) %>%
  summarise(IKI = median(IKI)) %>%
  spread(event, IKI) %>%
  ggplot(aes(x = KEYBOARD_down, y= KEYBOARD_firstkey, color = structure, linetype = struggling)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_point() +
  geom_smooth(method = "lm", fullrange = T, se = F) 


d3 %>%
  ggplot(aes(y = IKI, x = structure)) +
  facet_grid(~struggling) +
  geom_boxplot() +
  geom_violin() +
  scale_y_log10()
  

d3 %>%
  ggplot(aes(y = IKI, x = structure)) +
  facet_grid(~struggling) +
  stat_summary(fun.data = "mean_cl_normal",
               geom = "crossbar",
               width = 0.1,
               col = "red") 


d2 %>%
#  filter(!is.na(structure), ld < 9) %>%
#  filter(IKI < 10000,
#         IKI > 0) %>%
  mutate(struggling = factor(ifelse(struggling == TRUE, "struggling writers", "not struggling writers")),
         structure = factor(structure)) %>%
#  select(trial_id, word_id)
  group_by(trial_id) %>%
  mutate(max_word_id = max(word_id)) %>%
  filter(word_id > 0, max_word_id == 8) %>%
  ggplot(aes(y=IKI, x = word_id, fill = structure)) +
  facet_grid(~struggling) +
  stat_summary(fun.data = "mean_cl_normal",
               geom = "crossbar",
               width = 0.1,
               col = "red") 

d3 %>% select(trial_id, ld, text, target) %>%
  filter(ld > 5)


d3 %>% filter(ld < 10) -> d4
m <- lme4::lmer(log(IKI) ~ structure + (1|subj) + (1|item), d4)
m %>% summary() %>% coef() %>%
  round(2)

           
