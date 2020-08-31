library(stringdist)
library(tidyverse)
library(magrittr)

condnames <-  c("simple-complex", "complex-simple")

(files <- dir("results/exp1"))
length(files)

d <- map2(files, files, ~read_delim(paste0("results/exp1/",.,"/RESULTS_FILE.txt"), delim = "\t") %>%
            mutate(subj = .y)) %>%
  bind_rows() %>% 
  filter(structure %in% condnames) %>%
  select(-Session_Name_, -ends_with("x"), -ends_with("y"), -ends_with("file"), -ends_with("vis"),-list,-type)


d %<>%
  mutate(structure = ifelse(subj %in% 1:59 & structure == "simple-complex" & movement == "up", "complex-simple", structure)) 


d %<>%
  mutate(target = ifelse(structure == "simple-complex" & movement == "up",
                         paste0(p1_image, " moved ", "above the ", p2_image, " and the ", p3_image),
                         ifelse(structure == "simple-complex" & movement != "up",
                                paste0(p1_image, " moved ", "below the ", p2_image, " and the ", p3_image),
                                ifelse(structure == "complex-simple" & movement == "up", 
                                       paste0(p1_image, " and the", p2_image, " moved above the ", p3_image),
                                       ifelse(structure == "complex-simple" & movement != "up", 
                                              paste0(p1_image, " and the ", p2_image, " moved below the ", p3_image), NA))))) %>%
  mutate(text = tolower(text),
    target = tolower(target),
    text = trimws(text),
    ld = stringdist(text, target, method = "dl")) %>%
  select(-ends_with("image") )

d %>%
  filter(onset > 100) %>%
  ggplot(aes(y = onset, x = structure, fill = movement)) +
  geom_boxplot() +
  scale_y_log10()

summary(d$onset)

m <- lme4::lmer(log(onset) ~ structure + (structure|subj) + (structure|item), d, subset = onset > 0)

m %>% summary() %>% coef()
 
coefs <- m %>% coef()
coefs$subj %>% rownames_to_column("subj") %>%
  rename(Int = `(Intercept)`,
         Slope = `structuresimple-complex`) %>%
  mutate(Effect = exp(Int) - exp(Int + Slope)) %>%
  ggplot(aes(x = Effect)) +
  geom_density()

coefs$item %>% rownames_to_column("item") %>%
  rename(Int = `(Intercept)`,
         Slope = `structuresimple-complex`) %>%
  mutate(Effect = exp(Int) - exp(Int + Slope)) %>%
  ggplot(aes(x = Effect)) +
  geom_density()
