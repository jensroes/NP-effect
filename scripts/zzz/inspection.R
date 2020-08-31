library(tidyverse)

d <- read_csv("dfs/MartinEtAl2010.csv")

d$position <- factor(d$position)
d$ExpID <- factor(d$ExpID)
d$structure <- factor(d$structure)

contrasts(d$structure) <- contr.sum(2)
contrasts(d$position) <- contr.sum(2)
contrasts(d$ExpID) <- contr.sum(4)

m <- lme4::lmer(log(onset) ~ structure*ExpID*position + (1|ExpID/subj) + (1|item), d)
anova(m)
m %>% summary() %>% coef()


d %>% ggplot(aes(y = onset, x = structure, colour = position)) +
  facet_grid(~Exp) +
  geom_jitter(size =.05, 
              position = position_jitterdodge(jitter.width = .2, 
                                              dodge.width = .5)) +
  scale_y_log10() +
  geom_boxplot(outlier.shape = NA, width =.2, position = position_dodge(.5))

d %>% 
  ggplot(aes(y = onset, x = Exp, colour = position)) +
  geom_jitter(size =.05, 
              position = position_jitterdodge(jitter.width = .2, 
                                              dodge.width = .5)) +
  scale_y_log10() +
  geom_boxplot(outlier.shape = NA, width =.2, position = position_dodge(.5))


d %>% 
  ggplot(aes(y = onset, x = structure)) +
  geom_jitter(size =.05, width = .1) +
  scale_y_log10() +
  geom_boxplot(outlier.shape = NA, width =.2)


m <- lme4::lmer(log(onset) ~ structure*Exp + (1|Exp/subj) + (1|item), d)
m %>% summary() %>% coef()
anova(m)
