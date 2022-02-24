#SORT MICROBIAL MAPPIGNS INTO PUTATIVELY ANCIENT vs CONTAMINATION
#Locally check for microbial damage  - had trouble getting dependencies to run remotely
#Could easily set up on different super computing cluster
#I think there's an error somewhere in this script

#SOMETHING is WRONG WITH THIS SCRIPT
library(tidyverse)


#we're going to do this for file in list files
#download all of the mapdamage data for each ancient sample locally
setwd("~/Box Sync/Projects/aDNA_local/review_response/microbe_mapdam/")
all_mp <- list.files()
#fp <- "S17463.e2e.trim.filter30.Endolith_87.sorted.mapDamage/misincorporation.txt"
check_me <- function(fp){
tidy_misincorp <- read.table(fp, 
                             header = T, stringsAsFactors = F)
colnames(tidy_misincorp)[colnames(tidy_misincorp) == "T"] <- "Tt"
freqs_5p <- tidy_misincorp %>% filter(End=="5p", Pos <= 15) %>% group_by(Pos) %>%
  summarise(across(A:S, sum)) %>% mutate(ct_freq = C.T/C) %>% mutate(ag_freq = A.G/A, ga_freq = G.A/G, tc_freq = T.C/Tt, gc_freq = G.C/G, ac_freq = A.C/A, 
  at_freq = A.T/A, cg_freq = C.G/C, ca_freq= C.A/C, tg_freq = T.G/Tt, ta_freq=T.A/Tt, gt_freq = G.T/G) %>%
  mutate(firstpos = ifelse(Pos == 1, "Y", "N" ))
freqs_3p <- tidy_misincorp %>% filter(End=="3p", Pos <= 15) %>% group_by(Pos) %>%
  summarise(across(A:S, sum)) %>% mutate(ct_freq = C.T/C) %>% mutate(ag_freq = A.G/A, ga_freq = G.A/G, tc_freq = T.C/Tt, gc_freq = G.C/G, ac_freq = A.C/A, 
                                                                     at_freq = A.T/A, cg_freq = C.G/C, ca_freq= C.A/C, tg_freq = T.G/Tt, ta_freq=T.A/Tt, gt_freq = G.T/G) %>%
  mutate(firstpos = ifelse(Pos == 1, "Y", "N" ))




hard_checks_passed <- 0
soft_checks_passed <- 0

  if(is.na(freqs_5p$ct_freq[1]) == F & is.na(freqs_5p$ct_freq[2]) == F){
    if((freqs_5p$ct_freq[1] > freqs_5p$ct_freq[2]) == T){
      hard_checks_passed <- hard_checks_passed + 1
    }
    if(length(which((freqs_5p$ct_freq[1] > freqs_5p$ct_freq[3:10]) == T)) == 8){
      hard_checks_passed <- hard_checks_passed + 1
    }
  }
  if(is.na(freqs_3p$ga_freq[1]) == F & is.na(freqs_3p$ga_freq[2]) == F ){
    if((freqs_3p$ga_freq[1] > freqs_3p$ga_freq[2]) == T){
      hard_checks_passed <- hard_checks_passed + 1
    }
    if(length(which((freqs_3p$ga_freq[1] > freqs_3p$ga_freq[3:10]) == T)) == 8){
      hard_checks_passed <- hard_checks_passed + 1
    }
  }
  


  #max hard checks to pass = 4

  avg_mis5p <- freqs_5p %>% select(ag_freq:gt_freq)
  avg_mis3p <- freqs_3p %>% select(ct_freq, ag_freq, tc_freq:gt_freq)
  threshold <- function(df){
  m <- colMeans(df)
  s <- apply(df, 2, sd)
  return(m + 2*s)
  }

  high_check1 <- freqs_5p$ct_freq[1] > threshold(avg_mis5p) #first pos is higher than the rest as well 
  high_check2 <- freqs_3p$ct_freq[1] > threshold(avg_mis3p)
  together <- c(high_check1, high_check2)
  if(length(together[together == T]) > 18){
    soft_checks_passed <- soft_checks_passed + 1
  }
  
  low_check1 <- mean(freqs_5p$ct_freq[2:10]) > threshold(avg_mis5p) #first pos is higher than the rest as well 
  low_check2 <- mean(freqs_3p$ct_freq[2:10]) > threshold(avg_mis3p)
  together1 <- c(low_check1, low_check2)
  if(length(together1[together1 == T]) > 18){
    soft_checks_passed <- soft_checks_passed + 1
  }
  #max soft checks which can be passed = 2
  return(c(strsplit(fp, "[.]")[[1]][1], strsplit(fp, "[.]")[[1]][4], hard_checks_passed, soft_checks_passed))

}

results <- data.frame(name = NA, microbe = NA, hard_checks = NA, soft_checks = NA)
for(i in 1:length(all_mp)){
  x <- check_me(paste(all_mp[i], "/misincorporation.txt", sep = ''))
  results <- rbind(results, x)
}
results <- results[-1,]

#get your results, I set a relatively soft threshold of passing all hard checks
#and no soft checks. I'm not sure the soft check system is working
#I then manually examined these results to check what was putatively ancient
unique(results[results$hard_checks >= 4 & results$soft_checks >= 0,]$microbe)

#(2)  5’ C → T and 3’ G → A misincorporation rates higher at the terminal position 
#than the average misincorporation rates +/- 2 standard deviations at this position, 
#and (3) the C → T/G → A misincorporation rates must be less than the average misincorporation rates 
#+/- 2 standard deviations at the next ten base positions. R


