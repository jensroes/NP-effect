get_data <- function(file, print_summary = TRUE){
  data <- read_csv(file) 
  total <- length(unique(data$subj)) * length(unique(data$item))
  self_corr <- unique(data[data$key_id %in% c(-2,-1,0),]$trial_id)
  
  data %>% mutate(first_word = word(text, 1)) %>%
    select(trial_id, first_word, text) %>%
    filter(!(first_word %in% c("peter", "pater", "petter", "tanua", "oeter", "tnia", "tania,", "peter,",
                               "tania", "tana", "tnaia", "tainia", "tani", "taia", "tanie", "peer", "peper", "pateter", "5ania", "Oeter", "tan8a",
                               "taina", "tanian", "talia", "tanya", "tanai", "tanina", "0eter", "tasnia", "taniafridge", "[eter",
                               "tanis", "peterand", "taniaand"))) %>% unique() %>% pull(trial_id) -> incorrect_start
  
  data %>% filter( ld >= 10) %>% pull(trial_id) %>% unique() -> large_ld
  
  data %>% filter(key_id == 1, BACKSPACE == FALSE, onset != 0) %>%
    filter(!(trial_id %in% c(self_corr, incorrect_start, large_ld))) -> d2
  
  d2 %>% count(subj) %>% arrange(n) %>% filter(n < 10) %>% pull(subj) -> remove_subj
  
  d2 %>% filter(!(subj %in% remove_subj)) %>%
    mutate(subj = as.integer(factor(subj)),
           item = as.integer(factor(item))) %>%
    select(onset, subj, item, movement, structure, text) -> dfin
  
  if(print_summary){
    print(paste0("Self corrections: ", 
                 length(self_corr), " trials (", round(length(self_corr)/total*100,2) ,"%)"), quote = F)
    print(paste0("Incorrect first word: ", 
                 length(incorrect_start), " trials (", round(length(incorrect_start)/total*100,2) ,"%)"), quote = F)
    print(paste0("Final text Levenshtein distance >= 10: ", 
                 length(large_ld), " trial(s) (", round(length(large_ld)/total*100,2) ,"%)"), quote = F)
    print(paste0("Participants that contributed < 10 trials: ", length(remove_subj)), quote = F)
    print(paste0("Remaining participants: ", length(unique(dfin$subj))), quote = F)
  }
  return(dfin)
}
