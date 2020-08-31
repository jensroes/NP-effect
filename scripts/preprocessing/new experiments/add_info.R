library(dplyr)
library(stringr)
library(stringi)

alldata <- read.csv("data/alllogfiles.csv", stringsAsFactors = F) %>%
  mutate(key = ifelse(event %in% c("KEYBOARD_firstkey", "KEYBOARD_down") & 
                 is.na(key), " ", key))

#long ugly loop to get text so far
# for every trial do 
for(j in 1:max(alldata$trial_id)){
  newdf <- alldata %>%
    filter(trial_id == j,
           event %in% c("KEYBOARD_firstkey", "KEYBOARD_down"))
  newdf$trans_text <- ""
  
  # for every keystroke in a trial do
  for(i in 1:(nrow(newdf))){
    # check first keystroke is non character
    if(i==1 & newdf$key[i] %in% c("\033", "Capslock", "Down", "Up",
                                     "Insert",  "Ralt", "Rctrl", "Rshift",
                                     "Lalt", "Lctrl", "Lshift", "Backspace",
                                     "Delete", "Left", "Right")){
      newdf$trans_text[i] <- ""
    } else if (i==1) {
      newdf$trans_text[i] <- newdf$key[i]
    } else if(newdf$key[i] == "Backspace") {
      newdf$trans_text[i] <- substr(newdf$trans_text[i-1], 1, 
                                    nchar(newdf$trans_text[i-1])-1)
      # right and left are only used 151/1043 times respectively (out of 86995),
      # so we do not look at them (yet)
    } else if(newdf$key[i] %in% c("\033", "Capslock", "Down", "Up",
                                     "Insert",  "Ralt", "Rctrl", "Rshift",
                                    "Lalt", "Lctrl", "Lshift", "Backspace",
                                    "Delete", "Left", "Right")) {  
      newdf$trans_text[i] <- newdf$trans_text[i-1]
    } else {
      newdf$trans_text[i] <- paste0(newdf$trans_text[i-1], newdf$key[i])
    }
  } 
  if(j==1){
    newdf2 <- newdf
  }else{
    newdf2 <- rbind(newdf2, newdf)
  }
}

# get ambiguous info 
items <- alldata %>%
  filter(syn == "unambiguous") %>% 
  group_by(item) %>%
  select(item, target) %>%
  slice(1) %>%
  mutate(before_comma = gsub("^(.*?),.*", "\\1", target),
         before_verb = word(before_comma,-1),
         after_comma = sub("^(.*?),", "", target),
         after_np = ifelse(item == 60, word(after_comma, 1,4),
                           word(after_comma, 1,3)),
         after_verb = ifelse(item == 60, word(after_comma, 5),
                word(after_comma, 4)), 
         before_verb_position = which(grepl(before_verb,
                                            unlist(str_split(target, "\\s+")))),
         after_verb_position = which(grepl(after_verb,
                                           unlist(str_split(target, "\\s+"))))
         ) %>%
  select(-target)


# add extra features to df
data_extra <- alldata %>%
  left_join(newdf2) %>%
  left_join(items) %>%
  group_by(subj, item, trial_id) %>%
  mutate(trans_text = ifelse(is.na(trans_text),
                            lag(trans_text), trans_text),
         startword = (event == "KEYBOARD_firstkey" | 
                        (lag(event) == "KEYBOARD_down" & lag(key) == " ")),
         endword = (lead(event) == "KEYBOARD_down" & lead(key) == " "),
         n_trans_words = str_count(trans_text, "\\S+"),
         before_verb = (before_verb_position == n_trans_words),
         after_verb = (after_verb_position == n_trans_words),
         noun = (after_verb_position - 1 == n_trans_words),
         revision = (key == "Backspace" & !is.na(key)),
         revision_no = ifelse(revision == 1,
                                cumsum(revision-lag(revision) > 0 |
                                         row_number() == 1) -1,  0)
         )
           
           
write.csv(data_extra, "data/logfiles_extra.csv")       

           