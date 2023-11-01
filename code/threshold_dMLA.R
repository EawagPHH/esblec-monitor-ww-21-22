#Calculate the detection limit of dMLA reaction------------------------
##Load packages--------
library (ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)

##Formatting df---------------
#Import dataset
<<<<<<< HEAD
setwd("/Users/conforsh/switchdrive/Institution/Manuscripts/ESBL_ecoli_ww_monitoring")
=======
setwd("/Users/conforsh/switchdrive/Institution/Manuscripts/ESBL_ecoli_ww_monitoring/esblec-monitor-ww-21-22/data")
>>>>>>> 8c4ee92ba793b20d4a07b1b7748e35f9e4271a52
reads <- read.csv("sequencing_targets.csv", header = TRUE)
reads

# Create a reference dataframe with all levels of X2 present in the reads df
complete_X2 <- data.frame(X2 = levels(factor(reads$X2)))

##Calculate detection limit-----------------------
#Create a dataset with number of reads in negative controls of the reaction
neg = reads[reads$X1 %in% c("ACTGTGT", "TAGAACG", "GTTTCGG"), ]

#Compute the detection limit of the reads in the negative through the sum of the mean + 3 times the sd.
summary_stat = neg %>% 
  group_by(X2) %>% 
  summarise(reads_mean = mean(n), reads_sd = sd(n))%>%
  mutate(threshold = reads_mean + 3 * reads_sd)

# Perform left join to include all levels in summary_stat and fill other columns with NA and threshold=0
summary_stat <- left_join(complete_X2, summary_stat, by = "X2") %>%
  mutate(reads_mean = ifelse(is.na(reads_mean), NA, reads_mean),
         reads_sd = ifelse(is.na(reads_sd), NA, reads_sd),
         threshold = ifelse(is.na(threshold), NA, threshold)) %>%
          replace_na(list(threshold = 0))

# Left join reads with summary_stat on X2
reads <- left_join(reads, summary_stat, by = "X2")

# Calculate real_n based on the formula n - threshold, ensuring it is not lower than 0 and remove variables reads_mean and reads_sd
reads <- reads %>%
  mutate(real_n = ifelse(is.na(threshold), n, pmax(n - threshold, 0)))%>%
  select(-reads_mean, -reads_sd)

# Export data
write.csv(reads, file = "reads_counts_pre_processing.csv", row.names = FALSE)

