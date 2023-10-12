#-----ESBL-Ec in WW------
#1) Load required packages
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(geomtextpath)
library(tidyr)
library(stats)
library(tidyverse)
library(lubridate)
library(lme4)
library(plyr)
library(ggpubr)
library(gridExtra)
library(ggpmisc)
library(mgcv)
library(data.table)
library(scales)

##Formatting df-------------------
#Import datasets 
setwd("~/switchdrive/Institution/Manuscripts/ESBLEc_Monitoring_Pictures")
df = read.table("20231011_E_coli_counts.txt", header= TRUE)
df

#Convert the date column to a date format
df$date <- as.Date(df$date, format = "%d_%m_%Y")

#Create a new column for month
df$month <- as.factor(format(df$date, "%b"))

#Extract the year from the date column
df$year <- as.numeric(format(df$date, "%Y"))

#Create a new column for season based on specific dates
df$season <- ifelse(df$date >= as.Date(paste0(df$year, "-03-23")) & df$date < as.Date(paste0(df$year, "-06-23")), "Spring",
                    ifelse(df$date >= as.Date(paste0(df$year, "-06-23")) & df$date < as.Date(paste0(df$year, "-09-23")), "Summer",
                           ifelse(df$date >= as.Date(paste0(df$year, "-09-22")) & df$date < as.Date(paste0(df$year, "-12-22")), "Fall",
                                  "Winter")))

#Create a new column for seasons specific of different years
df$season_year <- paste(df$season, df$year, sep = " ")
df$month_year = paste(df$month, df$year, sep = " ")

#Change name of WWTPs
df$wwtp[which(df$wwtp=="Altenrhein")] <- "ARA Altenrhein"
df$wwtp[which(df$wwtp=="Lugano")] <- "IDA CDA Lugano"
df$wwtp[which(df$wwtp=="Geneva")] <- "STEP d'Aïre Genève"
df$wwtp[which(df$wwtp=="Laupen")] <- "ARA Sensetal Laupen"
df$wwtp[which(df$wwtp=="Zurich")] <- "ARA Werdhölzli Zürich"
df$wwtp[which(df$wwtp=="Chur")] <- "ARA Chur"

#Convert E.coli data as numeric
df$esblEc_percentage_a=as.numeric(df$esblEc_percentage_a)
df$esblEc_percentage_b=as.numeric(df$esblEc_percentage_b)
df$esblEc_loads_a=as.numeric(df$esblEc_loads_a)
df$esblEc_loads_b=as.numeric(df$esblEc_loads_b)
df$totalEc_loads_a=as.numeric(df$totalEc_loads_a)
df$totalEc_loads_b=as.numeric(df$totalEc_loads_b)
df$totalEc_cfu_100ml_a=as.numeric(df$totalEc_cfu_100ml_a)
df$totalEc_cfu_100ml_b=as.numeric(df$totalEc_cfu_100ml_b)
df$esblEc_cfu_100ml_a=as.numeric(df$esblEc_cfu_100ml_a)
df$esblEc_cfu_100ml_b=as.numeric(df$esblEc_cfu_100ml_b)

#Create an average value of the E.coli data replicates and use the non-missing value if one of the two replicates is missing
df$average_ESBL_Ec <- as.numeric(ifelse(is.na(df$esblEc_percentage_a) | is.na(df$esblEc_percentage_b), 
                                ifelse(is.na(df$esblEc_percentage_a), df$esblEc_percentage_b, df$esblEc_percentage_a), 
                                rowMeans(df[,c("esblEc_percentage_a", "esblEc_percentage_b")], na.rm = TRUE)))

df$average_loads_tot_Ec <- as.numeric(ifelse(is.na(df$totalEc_loads_a) | is.na(df$totalEc_loads_b), 
                                   ifelse(is.na(df$totalEc_loads_a), df$totalEc_loads_b, df$totalEc_loads_a), 
                                   rowMeans(df[,c("totalEc_loads_a", "totalEc_loads_b")], na.rm = TRUE)))

df$average_loads_ESBL_Ec <- as.numeric(ifelse(is.na(df$esblEc_loads_a) | is.na(df$esblEc_loads_b), 
                                  ifelse(is.na(df$esblEc_loads_a), df$esblEc_loads_b, df$esblEc_loads_a), 
                                  rowMeans(df[,c("esblEc_loads_a", "esblEc_loads_b")], na.rm = TRUE)))

df$counts_tot_Ec_1 = (df$totalEc_cfu_100ml_a)/100 #CFU/mL
df$counts_tot_Ec_2 =(df$totalEc_cfu_100ml_b)/100 #CFU/mL
df$average_counts_tot_Ec <- as.numeric(ifelse(is.na(df$counts_tot_Ec_1) | is.na(df$counts_tot_Ec_2), 
                                              ifelse(is.na(df$counts_tot_Ec_1), df$counts_tot_Ec_2, df$counts_tot_Ec_1), 
                                              rowMeans(df[,c("counts_tot_Ec_1", "counts_tot_Ec_2")], na.rm = TRUE)))

df$counts_ESBL_Ec_1 = (df$esblEc_cfu_100ml_a)/100 #CFU/mL
df$counts_ESBL_Ec_2 =(df$esblEc_cfu_100ml_b)/100 #CFU/mL
df$average_counts_ESBL_Ec <- as.numeric(ifelse(is.na(df$counts_ESBL_Ec_1) | is.na(df$counts_ESBL_Ec_2), 
                                               ifelse(is.na(df$counts_ESBL_Ec_1), df$counts_ESBL_Ec_2, df$counts_ESBL_Ec_1), 
                                               rowMeans(df[,c("counts_ESBL_Ec_1", "counts_ESBL_Ec_2")], na.rm = TRUE)))


#Subset df for WWTPs
df_sen = df %>% filter(wwtp == "ARA Sensetal Laupen")
df_lug = df %>% filter(wwtp == "IDA CDA Lugano")
df_alt = df %>% filter(wwtp == "ARA Altenrhein")
df_chu = df %>% filter(wwtp == "ARA Chur")
df_zur = df %>% filter(wwtp == "ARA Werdhölzli Zürich")
df_gen = df %>% filter(wwtp == "STEP d'Aïre Genève")

            # Function to calculate the mean for each bootstrap sample
            data = df_sen$average_ESBL_Ec
            mean_func <- function(data, indices) {
              mean(data[indices])
            }
            
            # Set the number of bootstrap samples
            num_boot_samples <- 1000
            
            # Perform bootstrapping to calculate the 95% CI for the mean
            boot_results <- boot(data, mean_func, R = num_boot_samples)
            
            # Calculate the percentile-based confidence interval
            confidence_interval <- quantile(boot_results$t, c(0.025, 0.975))
            
            # Print the confidence interval
            print(confidence_interval)

######SHEENA YOU NEED TO ADD THE VARIABLES OF COUNTS CFU/ML IN THE OVERALL STATISTICS CODE BELOW
##Calculate overall statistics------------------------
# Define a function to calculate mean, standard deviation, median, quantiles, minimum, and maximum
calculate_stats <- function(data, log_transform = FALSE) {
  if (log_transform) {
    data <- log10(data)
    unit <- "log10"
  } else {
    unit <- ""
  }
  
  mean_value <- round(mean(data, na.rm = TRUE), 1)
  sd_value <- round(sd(data, na.rm = TRUE), 1)
  median_value <- round(median(data, na.rm = TRUE), 1)
  quantile_values <- quantile(data, na.rm = TRUE)
  
  formatted_quantiles <- sapply(quantile_values, function(x) {
    if (x == 0) {
      return(0)
    } else {
      power <- floor(log10(abs(x)))
      mantissa <- round(x / 10^power, 1)
      return(paste(mantissa, "x 10^", power, unit))
    }
  })
  
  min_value <- min(data, na.rm = TRUE)
  max_value <- max(data, na.rm = TRUE)
  
  formatted_min <- ifelse(log_transform, paste("10^", sprintf("%.2f", min_value), unit), sprintf("%.2f", min_value))
  formatted_max <- ifelse(log_transform, paste("10^", sprintf("%.2f", max_value), unit), sprintf("%.2f", max_value))
  
  return(list(mean = mean_value, sd = sd_value, median = median_value, quantiles = formatted_quantiles, min = formatted_min, max = formatted_max))
}
# List of data frames
data_frames <- list(df = df, df_sen = df_sen, df_lug = df_lug, df_alt = df_alt, df_chu = df_chu, df_zur = df_zur, df_gen = df_gen)

# Loop through data frames and variables to calculate statistics
for (data_frame_name in names(data_frames)) {
  data_frame <- data_frames[[data_frame_name]]
  cat("Data Frame:", data_frame_name, "\n")
  
  # Calculate statistics for average_ESBL_Ec without log transformation
  stats_esbl <- calculate_stats(data_frame$average_ESBL_Ec)
  cat("ESBL-Ec Percentage - Mean:", stats_esbl$mean, "- SD:", stats_esbl$sd, "- Median:", stats_esbl$median, "- Quantiles:", stats_esbl$quantiles, "- Min:", stats_esbl$min, "- Max:", stats_esbl$max, "\n")
  
  # Calculate statistics for average_loads_ESBL_Ec and average_loads_tot_Ec with log transformation
  stats_loads_esbl <- calculate_stats(data_frame$average_loads_ESBL_Ec, TRUE)
  cat("ESBL-Ec Loads (log10 transformed) - Mean:", stats_loads_esbl$mean, "- SD:", stats_loads_esbl$sd, "- Median:", stats_loads_esbl$median, "- Quantiles:", stats_loads_esbl$quantiles, "- Min:", stats_loads_esbl$min, "- Max:", stats_loads_esbl$max, "\n")
  
  stats_loads_tot <- calculate_stats(data_frame$average_loads_tot_Ec, TRUE)
  cat("Total-Ec Loads (log10 transformed) - Mean:", stats_loads_tot$mean, "- SD:", stats_loads_tot$sd, "- Median:", stats_loads_tot$median, "- Quantiles:", stats_loads_tot$quantiles, "- Min:", stats_loads_tot$min, "- Max:", stats_loads_tot$max, "\n\n")
}

            
##Differences between months-------------------------------
###Plots----------------------
#percentage_ESBL_Ec
per1=ggplot(data=df) +
  geom_boxplot(aes(y=(average_ESBL_Ec), x=reorder(format(as.Date(date), "%b %Y"), date)), outlier.colour = NA, outlier.shape = NA) +
  geom_point(aes(y=(average_ESBL_Ec), x=reorder(format(as.Date(date), "%b %Y"), date), color = ifelse(is.na(esblEc_percentage_a) | is.na(esblEc_percentage_b), "Single replicate", "Averaged on two replicates"), shape = ifelse(is.na(esblEc_percentage_a) | is.na(esblEc_percentage_b), "Single replicate", "Averaged on two replicates"))) +
  facet_wrap(wwtp~., ncol=1) +
  theme(axis.text.x=element_text(angle=90, vjust=1, hjust=1), axis.title.y=element_text(size=11),legend.position="bottom") +
  ylab(expression(paste("Percentage of ESBL- ", italic("E. coli"), " (%)"))) +
  scale_color_manual(name="", values = c("Single replicate" = "darksalmon", "Averaged on two replicates" = "black"))+
  scale_shape_manual(name="", values = c("Single replicate" = 17, "Averaged on two replicates" = 16))+
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5))+
  xlab("")

#loads_ESBL_Ec
lre1=ggplot(data=df) +
  geom_boxplot(aes(y=log10(average_loads_ESBL_Ec), x=reorder(format(as.Date(date), "%b %Y"), date)), outlier.colour = NA, outlier.shape = NA) +
  geom_point(aes(y=log10(average_loads_ESBL_Ec), x=reorder(format(as.Date(date), "%b %Y"), date), color = ifelse(is.na(esblEc_loads_a) | is.na(esblEc_loads_b), "Unique replicate", "Averaged on two replicates"), shape = ifelse(is.na(esblEc_percentage_a) | is.na(esblEc_percentage_b), "Single replicate", "Averaged on two replicates")), size=1.5) +
  facet_wrap(wwtp~., ncol=1) +
  theme(axis.text.x=element_text(angle=90, vjust=1, hjust=1), axis.title.y=element_text(size=11),legend.position="none") +
  ylab(expression(paste("loads ESBL- ", italic("E. coli"), " log10(CFUs/(person-day))"))) +
  scale_color_manual(name="", values = c("Unique replicate" = "darksalmon", "Averaged on two replicates" = "black"))+
  scale_shape_manual(name="", values = c("Single replicate" = 17, "Averaged on two replicates" = 16))+
  scale_y_continuous(breaks = c(7, 8, 9), labels=c("7" = expression(10^{7}), "8" =expression(10^{8}),
                                                   "9" = expression(10^{9})))+
  xlab("")

#loads_tot_Ec
lto1=ggplot(data=df) +
  geom_boxplot(aes(y=log10(average_loads_tot_Ec), x=reorder(format(as.Date(date), "%b %Y"), date)), outlier.colour = NA, outlier.shape = NA) +
  geom_point(aes(y=log10(average_loads_tot_Ec), x=reorder(format(as.Date(date), "%b %Y"), date), color = ifelse(is.na(totalEc_loads_a) | is.na(totalEc_loads_b), "Unique replicate", "Averaged on two replicates"), shape = ifelse(is.na(esblEc_percentage_a) | is.na(esblEc_percentage_b), "Single replicate", "Averaged on two replicates")), size=1.5) +
  facet_wrap(wwtp~., ncol=1) +
  theme(axis.text.x=element_text(angle=90, vjust=1, hjust=1), axis.title.y=element_text(size=11), legend.position="none") +
  ylab(expression(paste("loads total ", italic("E. coli"), " log10(CFUs/(person-day))"))) +
  scale_color_manual(name="", values = c("Unique replicate" = "darksalmon", "Averaged on two replicates" = "black"))+
  scale_shape_manual(name="", values = c("Single replicate" = 17, "Averaged on two replicates" = 16))+
  scale_y_continuous(breaks = c(9, 10, 11), labels=c("9" = expression(10^{9}), "10" =expression(10^{10}),
                                                     "11" = expression(10^{11})))+
  xlab("")

ggarrange(per1, lto1, lre1, ncol = 3, labels = c("A","B", "C"), common.legend = TRUE) #aggregate the three plots together

###ESBL-Ec to compare between wwtps------
#percentage_ESBL_Ec
#calculate median in each WWTP 
df_sen = df %>% filter(wwtp == "ARA Sensetal Laupen")
median(df_sen$average_ESBL_Ec, na.rm = T)

df_lug = df %>% filter(wwtp == "IDA CDA Lugano")
median(df_lug$average_ESBL_Ec, na.rm = T)

df_alt = df %>% filter(wwtp == "ARA Altenrhein")
median(df_alt$average_ESBL_Ec, na.rm = T)

df_chu = df %>% filter(wwtp == "ARA Chur")
median(df_chu$average_ESBL_Ec, na.rm = T)

df_zur = df %>% filter(wwtp == "ARA Werdhölzli Zürich")
median(df_zur$average_ESBL_Ec, na.rm = T)

df_gen = df %>% filter(wwtp == "STEP d'Aïre Genève")
median(df_gen$average_ESBL_Ec, na.rm = T)

m_o = c("STEP d'Aïre Genève", "IDA CDA Lugano","ARA Werdhölzli Zürich","ARA Altenrhein", "ARA Chur", "ARA Sensetal Laupen") #ordering WWTPs from highest to lowest percentage.
df$wwtp <- factor(df$wwtp, levels = m_o)

perpl=ggplot(data=df) +
  geom_boxplot(aes(y=(average_ESBL_Ec), x=wwtp), outlier.colour = NA, outlier.shape = NA) +
  geom_point(aes(y=(average_ESBL_Ec), x=wwtp, 
                 color = ifelse(is.na(esblEc_percentage_a) | is.na(esblEc_percentage_b), "Single replicate", "Averaged on two replicates"), 
                 shape = ifelse(is.na(esblEc_percentage_a) | is.na(esblEc_percentage_b), "Single replicate", "Averaged on two replicates")),
             position=position_jitter(width=0.1), size=1.9,alpha=0.6) +
  theme(axis.text.x=element_text(angle=20, vjust=0.8, hjust=0.9, colour = "black"), axis.title.y=element_text(size=8.5),axis.text.y=element_text(colour="black"), legend.position="bottom") +
  ylab(expression(paste("Percentage of ESBL- ", italic("E. coli"), " (%)"))) +
  scale_color_manual(name="", values = c("Single replicate" = "darksalmon", "Averaged on two replicates" = "black"))+
  scale_shape_manual(name="", values = c("Single replicate" = 17, "Averaged on two replicates" = 16))+
  scale_y_continuous(breaks = c(1, 2, 3, 4,5))+
  xlab("")

#loads_ESBL_Ec
#calculate median in each WWTP 
median(df_sen$average_loads_ESBL_Ec, na.rm = T)

median(df_lug$average_loads_ESBL_Ec, na.rm = T)

median(df_alt$average_loads_ESBL_Ec, na.rm = T)

median(df_chu$average_loads_ESBL_Ec, na.rm = T)

median(df_zur$average_loads_ESBL_Ec, na.rm = T)

median(df_gen$average_loads_ESBL_Ec, na.rm = T)

m_t = c("ARA Altenrhein","ARA Werdhölzli Zürich","STEP d'Aïre Genève", "ARA Sensetal Laupen", "ARA Chur", "IDA CDA Lugano") #ordering WWTPs from highest to lowest ESBL-Ec loads

lrepl=ggplot(data=df) +
  geom_boxplot(aes(y=log10(average_loads_ESBL_Ec), x=wwtp), outlier.colour = NA, outlier.shape = NA) +
  geom_point(aes(y=log10(average_loads_ESBL_Ec), x=wwtp,
                 color = ifelse(is.na(esblEc_loads_a) | is.na(esblEc_loads_b), "Unique replicate", "Averaged on two replicates"), 
                 shape = ifelse(is.na(esblEc_percentage_a) | is.na(esblEc_percentage_b), "Single replicate", "Averaged on two replicates")), 
            position=position_jitter(width=0.1), size=1.9,alpha=0.6) +
  theme(axis.text.x=element_text(angle=20, vjust=0.8, hjust=0.9, colour = "black"), axis.title.y=element_text(size=8.5),axis.text.y=element_text(colour="black"), legend.position="none") +
  ylab(expression(paste("loads ESBL- ", italic("E. coli"), " log10(CFUs/(person-day))"))) +
  scale_color_manual(name="", values = c("Unique replicate" = "darksalmon", "Averaged on two replicates" = "black"))+
  scale_shape_manual(name="", values = c("Single replicate" = 17, "Averaged on two replicates" = 16))+
  scale_y_continuous(breaks = c(7, 8, 9), labels=c("7" = expression(10^{7}), "8" =expression(10^{8}),
                                                   "9" = expression(10^{9})))+
  xlab("")

#loads_tot_Ec
#calculate median in each WWTP 
median(df_sen$average_loads_tot_Ec, na.rm = T)

median(df_lug$average_loads_tot_Ec, na.rm = T)

median(df_alt$average_loads_tot_Ec, na.rm = T)

median(df_chu$average_loads_tot_Ec, na.rm = T)

median(df_zur$average_loads_tot_Ec, na.rm = T)

median(df_gen$average_loads_tot_Ec, na.rm = T)

m_l = c("ARA Altenrhein", "ARA Chur","ARA Sensetal Laupen","ARA Werdhölzli Zürich", "STEP d'Aïre Genève", "IDA CDA Lugano") #ordering WWTPs from highest to lowest total Ec loads

ltopl=ggplot(data=df) +
  geom_boxplot(aes(y=log10(average_loads_tot_Ec), x=wwtp), outlier.colour = NA, outlier.shape = NA) +
  geom_point(aes(y=log10(average_loads_tot_Ec), x=wwtp, 
                 color = ifelse(is.na(totalEc_loads_a) | is.na(totalEc_loads_b), "Unique replicate", "Averaged on two replicates"), 
                 shape = ifelse(is.na(esblEc_percentage_a) | is.na(esblEc_percentage_b), "Single replicate", "Averaged on two replicates")), 
             position=position_jitter(width=0.1), size=1.9, alpha=0.6) +
  theme(axis.text.x=element_text(angle=20, vjust=0.8, hjust=0.9, colour = "black"), axis.title.y=element_text(size=8.5),axis.text.y=element_text(colour="black"), legend.position="none") +
  ylab(expression(paste("loads total ", italic("E. coli"), " log10(CFUs/(person-day))"))) +
  scale_color_manual(name="", values = c("Unique replicate" = "darksalmon", "Averaged on two replicates" = "black"))+
  scale_shape_manual(name="", values = c("Single replicate" = 17, "Averaged on two replicates" = 16))+
  scale_y_continuous(breaks = c(9, 10, 11), labels=c("9" = expression(10^{9}), "10" =expression(10^{10}),
                                                     "11" = expression(10^{11})))+
  xlab("")

ggarrange(perpl, ggarrange(ltopl, lrepl, ncol = 1, labels = c("B", "C")), 
          labels = c("A",""), ncol = 2)

### calculate CDF -------
#Plot observations of ESBL-Ec to compare between wwtps 
custom_palette <- c("#fc8d62ff", "#8da0cbff", "#e78ac3ff", "#ffd92fff", "#a6d854ff", "#66c2a5ff") #line colours for plots

#percentage_ESBL_Ec
df = filter(df, !is.na(average_ESBL_Ec)) #remove na.values from dataset

CDF_list <- list() #Create a list to store CDF data for each wwtp

for (wwtp in unique(df$wwtp)) { 
  CDF_list[[wwtp]] <- ecdf(df$average_ESBL_Ec[df$wwtp == wwtp])
} # Iterate through each wwtp and calculate its CDF using ecdf()

df_CDF <- data.frame(x = seq(0, max(df$average_ESBL_Ec), length.out = 1000)) # Combine CDF data for all wwtps into a single data frame
for (wwtp in unique(df$wwtp)) {
  df_CDF[, wwtp] <- CDF_list[[wwtp]](df_CDF$x)
}

df_CDF_long <- tidyr::gather(df_CDF, key = "wwtp", value = "CDF", -x) # Convert data from wide to long format for ggplot2

# Plot CDF curves for each wwtp using ggplot2
cdfper=ggplot(df_CDF_long, aes(x = x, y = CDF, color = wwtp)) +
  geom_line(linewidth=1.2) +
  theme_minimal()+
  theme(legend.position= c(0.99, 0.9),legend.justification = c(1, 1))+
  xlab(expression(paste("Percentage of ESBL- ", italic("E. coli"), " (%)"))) +
  ylab("Cumulative Distribution Function") +
  theme(axis.title.y = element_text(size=8.5), axis.title.x = element_text(size=9), axis.text.y = element_text(colour = "black"), axis.text.x = element_text(colour = "black"))+
  scale_color_manual(values = custom_palette)+
  labs(color = "WWTP")

#loads_ESBL_Ec
df = filter(df, !is.na(average_loads_ESBL_Ec)) #remove na.values from dataset

CDF_list <- list() #Create a list to store CDF data for each wwtp

for (wwtp in unique(df$wwtp)) {
  CDF_list[[wwtp]] <- ecdf(log10(df$average_loads_ESBL_Ec)[df$wwtp == wwtp])
} # Iterate through each wwtp and calculate its CDF using ecdf()

df_CDF <- data.frame(x = seq(6.5, max(log10(df$average_loads_ESBL_Ec)), length.out = 1000)) # Combine CDF data for all wwtps into a single data frame
for (wwtp in unique(df$wwtp)) {
  df_CDF[, wwtp] <- CDF_list[[wwtp]](df_CDF$x)
}

df_CDF_long <- tidyr::gather(df_CDF, key = "wwtp", value = "CDF", -x) # Convert data from wide to long format for ggplot2

# Plot CDF curves for each wwtp using ggplot2
cdfesbl=ggplot(df_CDF_long, aes(x = x, y = CDF, color = wwtp)) +
  geom_line(linewidth=1.2) +
  theme_minimal()+
  theme(legend.position= c(0.3, 0.9),legend.justification = c(1, 1))+
  xlab(expression(paste("loads ESBL- ", italic("E. coli"), " log10(CFUs/(person-day))"))) +
  ylab("Cumulative Distribution Function") +
  scale_color_manual(values = custom_palette)+
  labs(color = "WWTP")+
  scale_x_continuous(breaks = c(7, 8, 9), labels=c("7" = expression(10^{7}), "8" =expression(10^{8}),
                                                                       "9" = expression(10^{9}))) +
  theme(axis.title.y = element_text(size=8.5), axis.title.x = element_text(size=9), axis.text.y = element_text(colour = "black"), axis.text.x = element_text(colour = "black"))

#loads_total_Ec
df = filter(df, !is.na(log10(average_loads_tot_Ec))) #remove na.values from dataset

CDF_list <- list()# Create a list to store CDF data for each wwtp

for (wwtp in unique(df$wwtp)) { 
  CDF_list[[wwtp]] <- ecdf(log10(df$average_loads_tot_Ec)[df$wwtp == wwtp])
} # Iterate through each wwtp and calculate its CDF using ecdf()

df_CDF <- data.frame(x = seq(8.5, max(log10(df$average_loads_tot_Ec)), length.out = 1000)) # Combine CDF data for all wwtps into a single data frame
for (wwtp in unique(df$wwtp)) {
  df_CDF[, wwtp] <- CDF_list[[wwtp]](df_CDF$x)
}

df_CDF_long <- tidyr::gather(df_CDF, key = "wwtp", value = "CDF", -x) # Convert data from wide to long format for ggplot2

# Plot CDF curves for each wwtp using ggplot2
cdftot=ggplot(df_CDF_long, aes(x = x, y = CDF, color = wwtp)) +
  geom_line(linewidth=1.2) +
  theme_minimal()+
  theme(legend.position= c(0.3, 0.9),legend.justification = c(1, 1))+
  xlab(expression(paste("loads total ", italic("E. coli"), " log10(CFUs/(person-day))"))) +
  ylab("Cumulative Distribution Function") +
  scale_color_manual(values = custom_palette)+
  labs(color = "WWTP")+
  scale_x_continuous(breaks = c(9, 10, 11), labels=c("9" = expression(10^{9}), "10" =expression(10^{10}),
                                                     "11" = expression(10^{11})))+
  theme(axis.title.y = element_text(size=8.5), axis.title.x = element_text(size=9), axis.text.y = element_text(colour = "black"), axis.text.x = element_text(colour = "black"))

#####Plot CDF and Boxplots of WWTP differences----------------
ggarrange(perpl, cdfper, ltopl, cdftot, lrepl, cdfesbl,nrow=3, ncol=2, labels=c("A", "B", "C", "D", "E", "F"), common.legend = TRUE)

### Environmental variables-------
df$temperature=as.numeric(df$temperature)
df$precipitations_24h_sum_mm=as.numeric(df$precipitations_24h_sum_mm)
df$precipitations_96h_sum_mm=as.numeric(df$precipitations_96h_sum_mm)
df$wwtp <- as.factor(df$wwtp)
unique_locations <- unique(df$wwtp) # Unique locations in your dataset

####ESBL-Ec percentage---------------------
#####Temperature-------------------------------
  correlation_results <- list()
  alpha <- 0.05  # Set your desired significance level (e.g., 0.05)
  
  for (loc in unique_locations) {
    data_loc <- subset(df, df$wwtp == loc)  # Subset data for the specific location
    
    correlation_result <- cor.test(data_loc$average_ESBL_Ec, data_loc$temperature, method = "spearman")
    
    # Apply Bonferroni correction
    p_value_corrected <- correlation_result$p.value * length(unique_locations)
    significant <- p_value_corrected <= alpha
    
    correlation_results[[loc]] <- list(
      correlation_coefficient = correlation_result$estimate,
      p_value = correlation_result$p.value,
      p_value_corrected = p_value_corrected,
      significant = significant
    )
  }
  
  # Access the correlation coefficients, p-values, and corrected p-values for each location
  for (loc in unique_locations) {
    result <- correlation_results[[loc]]
    correlation_coefficient <- result$correlation_coefficient
    p_value <- result$p_value
    p_value_corrected <- result$p_value_corrected
    significant <- result$significant
    
    print(paste("Correlation coefficient for", loc, ":", correlation_coefficient))
    print(paste("Raw p-value for", loc, ":", p_value))
    print(paste("Bonferroni-corrected p-value for", loc, ":", p_value_corrected))
    print(paste("Significant for", loc, ":", significant))
  }

#####Precipitation 24h------------------------------ 
  correlation_results <- list()
  alpha <- 0.05  # Set your desired significance level (e.g., 0.05)
  
  for (loc in unique_locations) {
    data_loc <- subset(df, df$wwtp == loc)  # Subset data for the specific location
    
    correlation_result <- cor.test(data_loc$average_ESBL_Ec, data_loc$precipitations_24h_sum_mm, method = "spearman")
    
    # Apply Bonferroni correction
    p_value_corrected <- correlation_result$p.value * length(unique_locations)
    significant <- p_value_corrected <= alpha
    
    correlation_results[[loc]] <- list(
      correlation_coefficient = correlation_result$estimate,
      p_value = correlation_result$p.value,
      p_value_corrected = p_value_corrected,
      significant = significant
    )
  }
  
  # Access the correlation coefficients, p-values, and corrected p-values for each location
  for (loc in unique_locations) {
    result <- correlation_results[[loc]]
    correlation_coefficient <- result$correlation_coefficient
    p_value <- result$p_value
    p_value_corrected <- result$p_value_corrected
    significant <- result$significant
    
    print(paste("Correlation coefficient for", loc, ":", correlation_coefficient))
    print(paste("Raw p-value for", loc, ":", p_value))
    print(paste("Bonferroni-corrected p-value for", loc, ":", p_value_corrected))
    print(paste("Significant for", loc, ":", significant))
  }

#####Precipitation 96h------------------------------------ 
  correlation_results <- list()
  alpha <- 0.05  # Set your desired significance level (e.g., 0.05)
  
  for (loc in unique_locations) {
    data_loc <- subset(df, df$wwtp == loc)  # Subset data for the specific location
    
    correlation_result <- cor.test(data_loc$average_ESBL_Ec, data_loc$precipitations_96h_sum_mm, method = "spearman")
    
    # Apply Bonferroni correction
    p_value_corrected <- correlation_result$p.value * length(unique_locations)
    significant <- p_value_corrected <= alpha
    
    correlation_results[[loc]] <- list(
      correlation_coefficient = correlation_result$estimate,
      p_value = correlation_result$p.value,
      p_value_corrected = p_value_corrected,
      significant = significant
    )
  }
  
  # Access the correlation coefficients, p-values, and corrected p-values for each location
  for (loc in unique_locations) {
    result <- correlation_results[[loc]]
    correlation_coefficient <- result$correlation_coefficient
    p_value <- result$p_value
    p_value_corrected <- result$p_value_corrected
    significant <- result$significant
    
    print(paste("Correlation coefficient for", loc, ":", correlation_coefficient))
    print(paste("Raw p-value for", loc, ":", p_value))
    print(paste("Bonferroni-corrected p-value for", loc, ":", p_value_corrected))
    print(paste("Significant for", loc, ":", significant))
  }

#####Plot heatmap--------------------------------------
# Generate df variables
correlation_coefficients <- c(0.394534630931989, 0.450110865205845, 0.234993325880179, 0.320655248832457, -0.149250864559337, 0.160802567895004,
                              0.12032107193718,-0.358154604398249, -0.115771803283113, -0.281567040706516, 0.0719901192324829, -0.258667627619305,
                              0.231378437899986,-0.211805457094709, -0.203158048829787, -0.0303622030279818, 0.289718938645438, -0.105004537864799)
p_values <- c(0.0274760373394784, 0.0055169435967361, 0.581584596772695, 0.130704485988326, 1.80547628831787, 1.68147485999875,
              2.43140657888691, 0.0591849647310879,2.51107406169802, 0.271940685990817, 3.71590778298805, 0.474940838462424,
              0.635613286812865, 0.814067106478304,0.916594316023928,4.99496999327398,0.247621185391402,2.89444421427645)
locations <- c("ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", "IDA CDA Lugano", "ARA Werdhölzli Zürich",
               "ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", "IDA CDA Lugano", "ARA Werdhölzli Zürich",
               "ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", "IDA CDA Lugano", "ARA Werdhölzli Zürich")
variables <- c( "Temperature (°C)","Temperature (°C)","Temperature (°C)","Temperature (°C)","Temperature (°C)","Temperature (°C)",
                "Precipitation (24h sum)","Precipitation (24h sum)","Precipitation (24h sum)","Precipitation (24h sum)","Precipitation (24h sum)","Precipitation (24h sum)",
                "Precipitation (96h sum)","Precipitation (96h sum)","Precipitation (96h sum)","Precipitation (96h sum)","Precipitation (96h sum)","Precipitation (96h sum)")

# Generate df for the heatmap
heatmap_data_1 <- data.frame(Location = locations, Variable = variables)
heatmap_data_1$Correlation <- correlation_coefficients
heatmap_data_1$P_Value <- p_values

# Plot heatmap
heatmap_plot <- ggplot(heatmap_data_1, aes(x = Variable, y = Location)) +
  geom_tile(aes(fill = Correlation), color = "white") +
  geom_text(data = subset(heatmap_data_1, P_Value < 0.05), aes(label = "*"), color = "black", size = 10) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "darkgreen", midpoint = 0, limits = c(-1, 1), na.value = "grey") +
  labs(x = "",
       y = "",
       fill = "Spearman's correlation coefficient") +
  theme_minimal() +
  ggtitle(expression(paste("Percentage of ESBL- ", italic("E. coli"), " (%)")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10), axis.text.y=element_text(size = 10))

####ESBL-Ec loads-------------------------
#####Temperature--------------------------------
correlation_results <- list()
alpha <- 0.05  # Set your desired significance level (e.g., 0.05)

for (loc in unique_locations) {
  data_loc <- subset(df, df$wwtp == loc)  # Subset data for the specific location
  
  correlation_result <- cor.test(data_loc$average_loads_ESBL_Ec, data_loc$temperature, method = "spearman")
  
  # Apply Bonferroni correction
  p_value_corrected <- correlation_result$p.value * length(unique_locations)
  significant <- p_value_corrected <= alpha
  
  correlation_results[[loc]] <- list(
    correlation_coefficient = correlation_result$estimate,
    p_value = correlation_result$p.value,
    p_value_corrected = p_value_corrected,
    significant = significant
  )
}

# Access the correlation coefficients, p-values, and corrected p-values for each location
for (loc in unique_locations) {
  result <- correlation_results[[loc]]
  correlation_coefficient <- result$correlation_coefficient
  p_value <- result$p_value
  p_value_corrected <- result$p_value_corrected
  significant <- result$significant
  
  print(paste("Correlation coefficient for", loc, ":", correlation_coefficient))
  print(paste("Raw p-value for", loc, ":", p_value))
  print(paste("Bonferroni-corrected p-value for", loc, ":", p_value_corrected))
  print(paste("Significant for", loc, ":", significant))
}

#####Precipitation 24h-----------------------
correlation_results <- list()
alpha <- 0.05  # Set your desired significance level (e.g., 0.05)

for (loc in unique_locations) {
  data_loc <- subset(df, df$wwtp == loc)  # Subset data for the specific location
  
  correlation_result <- cor.test(data_loc$average_loads_ESBL_Ec, data_loc$precipitations_24h_sum_mm, method = "spearman")
  
  # Apply Bonferroni correction
  p_value_corrected <- correlation_result$p.value * length(unique_locations)
  significant <- p_value_corrected <= alpha
  
  correlation_results[[loc]] <- list(
    correlation_coefficient = correlation_result$estimate,
    p_value = correlation_result$p.value,
    p_value_corrected = p_value_corrected,
    significant = significant
  )
}

# Access the correlation coefficients, p-values, and corrected p-values for each location
for (loc in unique_locations) {
  result <- correlation_results[[loc]]
  correlation_coefficient <- result$correlation_coefficient
  p_value <- result$p_value
  p_value_corrected <- result$p_value_corrected
  significant <- result$significant
  
  print(paste("Correlation coefficient for", loc, ":", correlation_coefficient))
  print(paste("Raw p-value for", loc, ":", p_value))
  print(paste("Bonferroni-corrected p-value for", loc, ":", p_value_corrected))
  print(paste("Significant for", loc, ":", significant))
}
  
#####Precipitation 96h-----------------------------------
correlation_results <- list()
alpha <- 0.05  # Set your desired significance level (e.g., 0.05)

for (loc in unique_locations) {
  data_loc <- subset(df, df$wwtp == loc)  # Subset data for the specific location
  
  correlation_result <- cor.test(data_loc$average_loads_ESBL_Ec, data_loc$precipitations_96h_sum_mm, method = "spearman")
  
  # Apply Bonferroni correction
  p_value_corrected <- correlation_result$p.value * length(unique_locations)
  significant <- p_value_corrected <= alpha
  
  correlation_results[[loc]] <- list(
    correlation_coefficient = correlation_result$estimate,
    p_value = correlation_result$p.value,
    p_value_corrected = p_value_corrected,
    significant = significant
  )
}

# Access the correlation coefficients, p-values, and corrected p-values for each location
for (loc in unique_locations) {
  result <- correlation_results[[loc]]
  correlation_coefficient <- result$correlation_coefficient
  p_value <- result$p_value
  p_value_corrected <- result$p_value_corrected
  significant <- result$significant
  
  print(paste("Correlation coefficient for", loc, ":", correlation_coefficient))
  print(paste("Raw p-value for", loc, ":", p_value))
  print(paste("Bonferroni-corrected p-value for", loc, ":", p_value_corrected))
  print(paste("Significant for", loc, ":", significant))
}

#####Plot heatmap--------------------------------------
# Generate df variables
correlation_coefficients <- c(0.41235232394182, 0.167156885098059, 0.425006223120863, 0.245808539607382, -0.153956972901298, 0.212057045803126,
                              -0.0131692074340034,-0.119732523543736, 0.341990819685203, -0.0268531529562696, -0.0405144749376039, 0.0903403710100557,
                              0.196686127344386,-0.092404198550732, 0.196509906857249, 0.0384100158787722, -0.113658582064334, 0.0460723074698408)
p_values <- c(0.0175414313648677, 1.44618145807016, 0.0112728163494952, 0.492533864128993, 1.71451557704347, 0.86117627910845,
              5.56605707223795, 2.41596064862275,0.0842211721812315, 5.10973176673581, 4.67986032970795, 3.2221348578552,
              1.02595922331321, 3.11387280009258,1.0016665410599,4.73403957505083,2.59151544699789,4.51954255053951)
locations <- c("ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", "IDA CDA Lugano", "ARA Werdhölzli Zürich",
               "ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", "IDA CDA Lugano", "ARA Werdhölzli Zürich",
               "ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", "IDA CDA Lugano", "ARA Werdhölzli Zürich")
variables <- c( "Temperature (°C)","Temperature (°C)","Temperature (°C)","Temperature (°C)","Temperature (°C)","Temperature (°C)",
                "Precipitation (24h sum)","Precipitation (24h sum)","Precipitation (24h sum)","Precipitation (24h sum)","Precipitation (24h sum)","Precipitation (24h sum)",
                "Precipitation (96h sum)","Precipitation (96h sum)","Precipitation (96h sum)","Precipitation (96h sum)","Precipitation (96h sum)","Precipitation (96h sum)")

#Generate df for heatmap
heatmap_data_2 <- data.frame(Location = locations, Variable = variables)
heatmap_data_2$Correlation <- correlation_coefficients
heatmap_data_2$P_Value <- p_values

# Plot heatmap
heatmap_plot_loads_ESBL <- ggplot(heatmap_data_2, aes(x = Variable, y = Location)) +
  geom_tile(aes(fill = Correlation), color = "white") +
  geom_text(data = subset(heatmap_data_2, P_Value < 0.05), aes(label = "*"), color = "black", size = 10) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "darkgreen", midpoint = 0, limits = c(-1, 1), na.value = "grey") +
  labs(x = "",
       y = "",
       fill = "Spearman's correlation coefficient") +
  theme_minimal() +
  ggtitle(expression(paste("Loads of ESBL- ", italic("E. coli"), " (CFUs/person-day)")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10), axis.text.y=element_text(size = 10))

####Total Ec loads------------------------
#####Temperature-----------------------
correlation_results <- list()
alpha <- 0.05  # Set your desired significance level (e.g., 0.05)

for (loc in unique_locations) {
  data_loc <- subset(df, df$wwtp == loc)  # Subset data for the specific location
  
  correlation_result <- cor.test(data_loc$average_loads_tot_Ec, data_loc$temperature, method = "spearman")
  
  # Apply Bonferroni correction
  p_value_corrected <- correlation_result$p.value * length(unique_locations)
  significant <- p_value_corrected <= alpha
  
  correlation_results[[loc]] <- list(
    correlation_coefficient = correlation_result$estimate,
    p_value = correlation_result$p.value,
    p_value_corrected = p_value_corrected,
    significant = significant
  )
}

# Access the correlation coefficients, p-values, and corrected p-values for each location
for (loc in unique_locations) {
  result <- correlation_results[[loc]]
  correlation_coefficient <- result$correlation_coefficient
  p_value <- result$p_value
  p_value_corrected <- result$p_value_corrected
  significant <- result$significant
  
  print(paste("Correlation coefficient for", loc, ":", correlation_coefficient))
  print(paste("Raw p-value for", loc, ":", p_value))
  print(paste("Bonferroni-corrected p-value for", loc, ":", p_value_corrected))
  print(paste("Significant for", loc, ":", significant))
}

#####Precipitation 24h-------------------------------
correlation_results <- list()
alpha <- 0.05  # Set your desired significance level (e.g., 0.05)

for (loc in unique_locations) {
  data_loc <- subset(df, df$wwtp == loc)  # Subset data for the specific location
  
  correlation_result <- cor.test(data_loc$average_loads_tot_Ec, data_loc$precipitations_24h_sum_mm, method = "spearman")
  
  # Apply Bonferroni correction
  p_value_corrected <- correlation_result$p.value * length(unique_locations)
  significant <- p_value_corrected <= alpha
  
  correlation_results[[loc]] <- list(
    correlation_coefficient = correlation_result$estimate,
    p_value = correlation_result$p.value,
    p_value_corrected = p_value_corrected,
    significant = significant
  )
}

# Access the correlation coefficients, p-values, and corrected p-values for each location
for (loc in unique_locations) {
  result <- correlation_results[[loc]]
  correlation_coefficient <- result$correlation_coefficient
  p_value <- result$p_value
  p_value_corrected <- result$p_value_corrected
  significant <- result$significant
  
  print(paste("Correlation coefficient for", loc, ":", correlation_coefficient))
  print(paste("Raw p-value for", loc, ":", p_value))
  print(paste("Bonferroni-corrected p-value for", loc, ":", p_value_corrected))
  print(paste("Significant for", loc, ":", significant))
}

#####Precipitation 96h------------------------------------
correlation_results <- list()
alpha <- 0.05  # Set your desired significance level (e.g., 0.05)

for (loc in unique_locations) {
  data_loc <- subset(df, df$wwtp == loc)  # Subset data for the specific location
  
  correlation_result <- cor.test(data_loc$average_loads_tot_Ec, data_loc$precipitations_96h_sum_mm, method = "spearman")
  
  # Apply Bonferroni correction
  p_value_corrected <- correlation_result$p.value * length(unique_locations)
  significant <- p_value_corrected <= alpha
  
  correlation_results[[loc]] <- list(
    correlation_coefficient = correlation_result$estimate,
    p_value = correlation_result$p.value,
    p_value_corrected = p_value_corrected,
    significant = significant
  )
}

# Access the correlation coefficients, p-values, and corrected p-values for each location
for (loc in unique_locations) {
  result <- correlation_results[[loc]]
  correlation_coefficient <- result$correlation_coefficient
  p_value <- result$p_value
  p_value_corrected <- result$p_value_corrected
  significant <- result$significant
  
  print(paste("Correlation coefficient for", loc, ":", correlation_coefficient))
  print(paste("Raw p-value for", loc, ":", p_value))
  print(paste("Bonferroni-corrected p-value for", loc, ":", p_value_corrected))
  print(paste("Significant for", loc, ":", significant))
}

#####Plot heatmap--------------------------------------
# Generate df variables
correlation_coefficients <- c(0.0696022090860442, -0.23340422667455, 0.274136254223402, 0.0276489355118023, -0.0533038801997631, 0.185608142014729,
                              -0.0631609765897568,0.162859177276606, 0.294373770503055, 0.0566784024533302, -0.086029154060545, 0.151651291036682,
                              0.0676838596499421,0.133381060408051, 0.214666211830527, -0.0286246070715612, -0.270842778373269, 0.0656060449636737)
p_values <- c(3.7646715316832, 0.59565821692541, 0.309451839594531, 5.0836590863394, 4.27883383324458, 1.26979700808417,
              3.95824346240189, 1.52108087008313,0.216045965848879, 4.15684481230589, 3.31493272259063, 1.85328403512669,
              3.82191110389342, 2.10462479237324,0.782113271419657,5.051729478814,0.342685834642609,3.9677364576379)
locations <- c("ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", "IDA CDA Lugano", "ARA Werdhölzli Zürich",
               "ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", "IDA CDA Lugano", "ARA Werdhölzli Zürich",
               "ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", "IDA CDA Lugano", "ARA Werdhölzli Zürich")
variables <- c( "Temperature (°C)","Temperature (°C)","Temperature (°C)","Temperature (°C)","Temperature (°C)","Temperature (°C)",
                "Precipitation (24h sum)","Precipitation (24h sum)","Precipitation (24h sum)","Precipitation (24h sum)","Precipitation (24h sum)","Precipitation (24h sum)",
                "Precipitation (96h sum)","Precipitation (96h sum)","Precipitation (96h sum)","Precipitation (96h sum)","Precipitation (96h sum)","Precipitation (96h sum)")

# Generate df for heatmap
heatmap_data_3 <- data.frame(Location = locations, Variable = variables)
heatmap_data_3$Correlation <- correlation_coefficients
heatmap_data_3$P_Value <- p_values

# Plot heatmap
heatmap_plot_loads_tot <- ggplot(heatmap_data_3, aes(x = Variable, y = Location)) +
  geom_tile(aes(fill = Correlation), color = "white") +
  geom_text(data = subset(heatmap_data_3, P_Value < 0.05), aes(label = "*"), color = "black", size = 10) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "darkgreen", midpoint = 0, limits = c(-1, 1), na.value = "grey") +
  labs(x = "",
       y = "",
       fill = "Spearman's correlation coefficient") +
  theme_minimal() +
  ggtitle(expression(paste("Loads of total ", italic("E. coli"), " (CFUs/person-day)")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10), axis.text.y=element_text(size = 10))

ggarrange(heatmap_plot,heatmap_plot_loads_ESBL, heatmap_plot_loads_tot, ncol=3, labels = c("A", "B", "C"), common.legend = TRUE)

###Statistical tests to check for difference in ESBL-monitoring of 1year----------------
df_filtered = filter(df, !is.na(average_ESBL_Ec)) #df_filtered contains the 300 observations for which ESBL-Ec percentage could be calculated
#shapiro test indicates that the "average_ESBL_Ec" data are not normally distributed
          #shapiro.test(df_filtered$average_ESBL_Ec)

          #Shapiro-Wilk normality test

          #W = 0.92056, p-value = 1.566e-11

#non-parametric test Kruskal-Wallis, 
  #which is used to compare the medians of two or more independent groups. 
  #The test does not assume normal distribution and is therefore suitable for 
  #non-normally distributed data.

# Load the necessary library
library(coin)

# Perform the Kruskal-Wallis test
df_filtered$wwtp=as.factor(df_filtered$wwtp)
kruskal_test(average_ESBL_Ec ~ wwtp, data = df_filtered)

  #Asymptotic Kruskal-Wallis Test
  #data:  average_ESBL_Ec by
  #wwtp (ARA Altenrhein, ARA Chur, ARA Sensetal Laupen, ARA Werdhölzli Z?rich, IDA CDA Lugano, STEP d'Aïre Genève)
  #chi-squared = 48.327, df = 5, p-value = 3.045e-09

# perform pairwise comparisons using post-hoc test Dunn's test with "Bonferroni" correction
  #to identify which groups differ significantly from each other.
library(dunn.test)
dunn_test_result <- dunn.test(df_filtered$average_ESBL_Ec, df_filtered$wwtp, 
                              method = "bonferroni")
dunn_test_result

#Check if there are different between months within each WWTP
df_per_alt = df_filtered %>% filter(wwtp == "ARA Altenrhein")
df_per_chu = df_filtered %>% filter(wwtp == "ARA Chur")
df_per_gen = df_filtered %>% filter(wwtp == "STEP d'Aïre Genève")
df_per_zur = df_filtered %>% filter(wwtp == "ARA Werdhölzli Zürich")
df_per_lug = df_filtered %>% filter(wwtp == "IDA CDA Lugano")
df_per_sen = df_filtered %>% filter(wwtp == "ARA Sensetal Laupen")

#calculate mean and sd of ESBL percentage in each wwtp
# Create a list of your data frames
df_list <- list(df_per_alt, df_per_chu, df_per_gen, df_per_zur, df_per_lug, df_per_sen)

# Create empty vectors to store the means and standard deviations
means <- c()
sds <- c()

# Loop through each data frame in the list
for (i in seq_along(df_list)) {
  # Calculate the mean and standard deviation for the variable of interest
  var_mean <- mean(df_list[[i]]$average_ESBL_Ec)
  var_sd <- sd(df_list[[i]]$average_ESBL_Ec)
  
  # Append the mean and standard deviation to their respective vectors
  means <- c(means, var_mean)
  sds <- c(sds, var_sd)
}

# Print the means and standard deviations for each data frame
cat("Means:", means, "\n")
cat("Standard deviations:", sds)

#check for normality in the ESBL-Ec percentage in each WWTP
shapiro.test(df_per_alt$average_ESBL_Ec) #Altenrhein
    #W = 0.9435, p-value = 0.01858
shapiro.test(df_per_chu$average_ESBL_Ec) #Chur
    #W = 0.79095, p-value = 4.477e-07
shapiro.test(df_per_gen$average_ESBL_Ec) #Geneva
    #W = 0.96207, p-value = 0.1022
shapiro.test(df_per_zur$average_ESBL_Ec) #Zurich
    #W = 0.84651, p-value = 2.094e-05
shapiro.test(df_per_lug$average_ESBL_Ec) #Lugano
    #W = 0.94392, p-value = 0.0193
shapiro.test(df_per_sen$average_ESBL_Ec) #Sensetal-Laupen
    #W = 0.92702, p-value = 0.003833

    #For all the WWTPs except one the normality assumption is not met. 
    #I will use a Kruskal-Wallis test with a Bonferroni adjusted Dunn's test to check for 
    #differences between months. Even if in geneva data are normally distributed, it is
    #better to use the same test as in the others WWTP to be able to compare results. 

#Perform Kruskal-wallis with Dunn's test and Bonferroni adjustmen as a loop on all the WWTP df to check for significant differences on a monthly scale
library(dunn.test)
df_list <- list(df_per_alt, df_per_chu, df_per_gen, df_per_sen, df_per_lug, df_per_zur)#Create a list of data frames for each wwtp

wwtp_names <- c("ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", 
                 "IDA CDA Lugano", "ARA Werdhölzli Zürich")#Define the names of the wwtps

# Perform the Kruskal-Wallis test and Dunn's test for each wwtp
library(writexl)

# Create an empty list to store the results
result_list <- list()

for (i in seq_along(df_list)) {
  # Perform the Kruskal-Wallis test
  kruskal_test <- kruskal.test(average_ESBL_Ec ~ month_year, data = df_list[[i]])
  
  # Perform the Dunn's test
  dunn_test <- dunn.test(x = df_list[[i]]$average_ESBL_Ec, g = df_list[[i]]$month_year, 
                         method = "bonferroni")
  
  # Create a data frame to store the results for each wwtp
  result_df <- data.frame(wwtp = wwtp_names[i],
                          Kruskal_Wallis_Test = capture.output(kruskal_test),
                          Dunn_Test = capture.output(dunn_test$comparison))
  
  # Append the result to the list
  result_list[[i]] <- result_df
}

# Save the results as an Excel file
write_xlsx(result_list, "results.xlsx")

# Confirmation message
cat("Results saved to results.xlsx\n")


for (i in seq_along(df_list)) {
  # Perform the Kruskal-Wallis test
  kruskal_test <- kruskal.test(average_ESBL_Ec ~ month_year, data = df_list[[i]])
  
  # Perform the Dunn's test
  dunn_test <- dunn.test(x = df_list[[i]]$average_ESBL_Ec, g = df_list[[i]]$month_year, 
                         method = "bonferroni")
  
  # Print the results for each wwtp
  cat("\n====================\n")
  cat(paste0("Results for ", wwtp_names[i], "\n"))
  cat("====================\n\n")
  
  # Print the Kruskal-Wallis test results
  cat("Kruskal-Wallis Test:\n")
  print(kruskal_test)
  cat("\n")
  
  # Print the Dunn's test results
  cat("Dunn's Test:\n")
  print(dunn_test$comparison)
  cat("\n")
}
# Close the file connection
close(file_conn)

# Confirmation message
cat("Results saved to", file_path, "\n")

#df_filtered_t contains the 301 observations for which loads of total-Ec could be calculated
df_filtered_t = filter(df, !is.na(average_loads_tot_Ec))
df_filtered_t$average_loads_tot_Ec= log(df_filtered_t$average_loads_tot_Ec)#log transformation of loads data

#shapiro test indicates that the "average_loads_tot_Ec" data are not normally distributed
#shapiro.test(df_filtered_t$average_loads_tot_Ec)
      #Shapiro-Wilk normality test

      #data:  df_filtered_t$average_loads_tot_Ec
      #W = 0.95418, p-value = 4.252e-08

#non-parametric test Kruskal-Wallis, 
#which is used to compare the medians of two or more independent groups. 
#The test does not assume normal distribution and is therefore suitable for 
#non-normally distributed data.

# Load the necessary library
library(coin)

# Perform the Kruskal-Wallis test
df_filtered_t$wwtp=as.factor(df_filtered_t$wwtp)
kruskal_test(average_loads_tot_Ec ~ wwtp, data = df_filtered_t)

#Asymptotic Kruskal-Wallis Test
#data:  average_loads_tot_Ec by
#wwtp (ARA Altenrhein, ARA Chur, ARA Sensetal Laupen, ARA Werdhölzli Z?rich, IDA CDA Lugano, STEP d'Aïre Genève)
#chi-squared = 40.537, df = 5, p-value = 1.163e-07

# perform pairwise comparisons using post-hoc test Dunn's test with "Bonferroni" correction
#to identify which groups differ significantly from each other.
library(dunn.test)
dunn_test_result <- dunn.test(df_filtered_t$average_loads_tot_Ec, df_filtered_t$wwtp, 
                              method = "bonferroni")
dunn_test_result

#filter df_filtered_t based on wwtp to find highest and lowest values in each WWTP
df_tot_alt = df_filtered_t %>% filter(wwtp == "ARA Altenrhein")
df_tot_chu = df_filtered_t %>% filter(wwtp == "ARA Chur")
df_tot_gen = df_filtered_t %>% filter(wwtp == "STEP d'Aïre Genève")
df_tot_zur = df_filtered_t %>% filter(wwtp == "ARA Werdhölzli Z?rich")
df_tot_lug = df_filtered_t %>% filter(wwtp == "IDA CDA Lugano")
df_tot_sen = df_filtered_t %>% filter(wwtp == "ARA Sensetal Laupen")

#calculate mean and sd of tot-Ec loads in each wwtp
# Create a list of your data frames
df_list <- list(df_tot_alt, df_tot_chu, df_tot_gen, df_tot_zur, df_tot_lug, df_tot_sen)

# Create empty vectors to store the means and standard deviations
means <- c()
sds <- c()

# Loop through each data frame in the list
for (i in seq_along(df_list)) {
  # Calculate the mean and standard deviation for the variable of interest
  var_mean <- mean(df_list[[i]]$average_loads_tot_Ec)
  var_sd <- sd(df_list[[i]]$average_loads_tot_Ec)
  
  # Append the mean and standard deviation to their respective vectors
  means <- c(means, var_mean)
  sds <- c(sds, var_sd)
}

# Print the means and standard deviations for each data frame
cat("Means:", means, "\n")
cat("Standard deviations:", sds)

#Perform Kruskal-wallis with Dunn's test and Bonferroni adjustmen as a loop on all the WWTP df to check for significant differences on a monthly scale
library(dunn.test)
df_list <- list(df_tot_alt, df_tot_chu, df_tot_gen, df_tot_sen, df_tot_lug, df_tot_zur)#Create a list of data frames for each wwtp

wwtp_names <- c("ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", 
                 "IDA CDA Lugano", "ARA Werdhölzli Zürich")#Define the names of the wwtps

# Perform the Kruskal-Wallis test and Dunn's test for each wwtp
for (i in seq_along(df_list)) {
  # Perform the Kruskal-Wallis test
  kruskal_test <- kruskal.test(average_loads_tot_Ec ~ month_year, data = df_list[[i]])
  
  # Perform the Dunn's test
  dunn_test <- dunn.test(x = df_list[[i]]$average_loads_tot_Ec, g = df_list[[i]]$month_year, 
                         method = "bonferroni")
  
  # Print the results for each wwtp
  cat("\n====================\n")
  cat(paste0("Results for ", wwtp_names[i], "\n"))
  cat("====================\n\n")
  
  # Print the Kruskal-Wallis test results
  cat("Kruskal-Wallis Test:\n")
  print(kruskal_test)
  cat("\n")
  
  # Print the Dunn's test results
  cat("Dunn's Test:\n")
  print(dunn_test$comparison)
  cat("\n")
}

#df_filtered_r contains the 302 observations for which loads of ESBL-Ec could be calculated
df_filtered_r = filter(df, !is.na(average_loads_ESBL_Ec))
df_filtered_r$average_loads_ESBL_Ec= log(df_filtered_r$average_loads_ESBL_Ec)#log transformation of loads data

#shapiro test indicates that the "average_loads_ESBL_Ec" data are not normally distributed
#shapiro.test(df_filtered_r$average_loads_ESBL_Ec)
#Shapiro-Wilk normality test

#data:  df_filtered_r$average_loads_ESBL_Ec
#W = 0.95078, p-value = 1.584e-08

#non-parametric test Kruskal-Wallis, 
#which is used to compare the medians of two or more independent groups. 
#The test does not assume normal distribution and is therefore suitable for 
#non-normally distributed data.

# Load the necessary library
library(coin)

# Perform the Kruskal-Wallis test
df_filtered_r$wwtp=as.factor(df_filtered_r$wwtp)
kruskal_test(average_loads_ESBL_Ec ~ wwtp, data = df_filtered_r)

#Asymptotic Kruskal-Wallis Test
#data:  average_loads_ESBL_Ec by
#wwtp (ARA Altenrhein, ARA Chur, ARA Sensetal Laupen, ARA Werdhölzli Z?rich, IDA CDA Lugano, STEP d'Aïre Genève)
#chi-squared = 36.239, df = 5, p-value = 8.507e-07


# perform pairwise comparisons using post-hoc test Dunn's test with "Bonferroni" correction
#to identify which groups differ significantly from each other.
library(dunn.test)
dunn_test_result <- dunn.test(df_filtered_r$average_loads_ESBL_Ec, df_filtered_r$wwtp, 
                              method = "bonferroni")
dunn_test_result

#filter df_filtered_t based on wwtp to find highest and lowest values in each WWTP
df_ESBL_alt = df_filtered_r %>% filter(wwtp == "ARA Altenrhein")
df_ESBL_chu = df_filtered_r %>% filter(wwtp == "ARA Chur")
df_ESBL_gen = df_filtered_r %>% filter(wwtp == "STEP d'Aïre Genève")
df_ESBL_zur = df_filtered_r %>% filter(wwtp == "ARA Werdhölzli Zürich")
df_ESBL_lug = df_filtered_r %>% filter(wwtp == "IDA CDA Lugano")
df_ESBL_sen = df_filtered_r %>% filter(wwtp == "ARA Sensetal Laupen")

#calculate mean and sd of tot-Ec loads in each wwtp
# Create a list of your data frames
df_list_r <- list(df_ESBL_alt, df_ESBL_chu, df_ESBL_gen, df_ESBL_zur, df_ESBL_lug, df_ESBL_sen)

# Create empty vectors to store the means and standard deviations
means <- c()
sds <- c()

# Loop through each data frame in the list
for (i in seq_along(df_list_r)) {
  # Calculate the mean and standard deviation for the variable of interest
  var_mean <- mean(df_list[[i]]$average_loads_ESBL_Ec)
  var_sd <- sd(df_list[[i]]$average_loads_ESBL_Ec)
  
  # Append the mean and standard deviation to their respective vectors
  means <- c(means, var_mean)
  sds <- c(sds, var_sd)
}

# Print the means and standard deviations for each data frame
cat("Means:", means, "\n")
cat("Standard deviations:", sds)

#Calculate min and max values of ESBL-Ec loads in each WWTP
min_values <- sapply(df_list_r, function(df) min(df$average_loads_ESBL_Ec))

# Apply the max() function to the 'x' variable in each data frame
max_values <- sapply(df_list_r, function(df) max(df$average_loads_ESBL_Ec))

# Print the results
min_values
max_values

#Perform Kruskal-wallis with Dunn's test and Bonferroni adjustmen as a loop on all the WWTP df to check for significant differences on a monthly scale
library(dunn.test)
df_list <- list(df_ESBL_alt, df_ESBL_chu, df_ESBL_gen, df_ESBL_sen, df_ESBL_lug, df_ESBL_zur)#Create a list of data frames for each wwtp

wwtp_names <- c("ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", 
                 "IDA CDA Lugano", "ARA Werdhölzli Zürich")#Define the names of the wwtps

# Perform the Kruskal-Wallis test and Dunn's test for each wwtp
for (i in seq_along(df_list)) {
  # Perform the Kruskal-Wallis test
  kruskal_test <- kruskal.test(average_loads_ESBL_Ec ~ month_year, data = df_list[[i]])
  
  # Perform the Dunn's test
  dunn_test <- dunn.test(x = df_list[[i]]$average_loads_ESBL_Ec, g = df_list[[i]]$month_year, 
                         method = "bonferroni")
  
  # Print the results for each wwtp
  cat("\n====================\n")
  cat(paste0("Results for ", wwtp_names[i], "\n"))
  cat("====================\n\n")
  
  # Print the Kruskal-Wallis test results
  cat("Kruskal-Wallis Test:\n")
  print(kruskal_test)
  cat("\n")
  
  # Print the Dunn's test results
  cat("Dunn's Test:\n")
  print(dunn_test$comparison)
  cat("\n")
}

###Carriage estimation------
####Whole Switzerland, with population adjusted ESBL-Ec percentage in WW-----
#Median and 25th-75th interquantile range of ESBL-Ec percentage in WW
  # ARA Alt = 1.57 (1.27-2.10)  64'000
  # ARA Chu = 1.29 (0.85-1.65)  55'000
  # ARA Gen = 2.07 (1.61-2.63)  454'000
  # ARA Zur = 1.69 (1.40-2.15)  471'000
  # ARA Lug = 1.89 (1.46-2.37)  124'000
  # ARA Sen = 1.26 (0.99-1.66)  62'000
                                  #tot population = 1'230'000

#Estimate median and interquantile ranges of ESBL-Ec adjusted by population in Switzerland
median = (1.57*64000+1.29*55000+2.07*454000+1.69*471000+1.89*124000+1.26*62000)/(64000+55000+454000+471000+124000+62000)
  #median = 1.81
min = (1.27*64000+0.85*55000+1.61*454000+1.40*471000+1.46*124000+0.99*62000)/(64000+55000+454000+471000+124000+62000)
  #min = 1.43
max = (2.1*64000+1.65*55000+2.63*454000+2.15*471000+2.37*124000+1.66*62000)/(64000+55000+454000+471000+124000+62000)
  #max = 2.3

# Create data frame for Switzerland
data_ch <- data.frame(x = rep(seq(0.01, 1, by = 0.0001)),  # 100 values for "x" repeated 6 times
                      a = rep(c(0.0181)))

data_ch_min <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                          a = c(0.0141))

data_ch_max <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                          a = c(0.0230))

# Calculate Y based on the function Y = a/x
data_ch$Y <- data_ch$a / data_ch$x
data_ch$Ymin <- data_ch_min$a / data_ch_min$x
data_ch$Ymax <- data_ch_max$a /data_ch_max$x

# Set Y values greater than 1 to 1
data_ch$Y[data_ch$Y > 1] <- 1
data_ch$Ymin[data_ch$Ymin > 1] <- 1
data_ch$Ymax[data_ch$Ymax > 1] <- 1

#Plot
ch <- ggplot(data_ch[data_ch$x >= 0.0181,], aes(x = x, y = Y)) +
  geom_line(color="black", linewidth=1) +
  geom_ribbon(data=data_ch, 
  aes(ymin = Ymin, ymax = Ymax), fill = "black", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c(0.005, 0.01,0.1, 1),  
    labels = c(0.005, 0.01,0.1, 1)) +  
  labs(x = "Proportion of ESBL-Ec out of total E. coli in the gut",
       y = "Prevalence of ESBL-Ec carriage within the community") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.018", vjust = 1.3, hjust=-0.0,data = data_ch[data_ch$Y == 1, ], 
                 aes(xintercept = a), linetype = "dashed", color="red", alpha = 0.8)+
  ggtitle("Switzerland")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ch

####All together-----
# Create a data frame with different values of "a", different one for each WWTP
data <- data.frame(x = rep(seq(0.01, 1, by = 0.01), times = 6),  # 100 values for "x" repeated 6 times
                   a = rep(c(0.0157, 0.0129, 0.0207, 0.0169, 0.0189, 0.0126), each = 100))  # Different values of "a" repeated for each "x"
#Alt, Chur, Genf, Zurich, Lugano, Laupen
# Calculate Y based on the function Y = a/x
data$Y <- data$a / data$x

# Set Y values greater than 1 to 1
data$Y[data$Y > 1] <- 1

# Adjust 'x' values based on 'a'
data$x <- ifelse(data$x < data$a, data$a, data$x)

# Create the plot
p=ggplot(data, aes(x = x, y = Y, color = factor(a))) +
  geom_line() +
  geom_point() +
  labs(x = "Proportion of ESBL-Ec out of total E. coli in the gut",
       y = "Prevalence of ESBL-Ec carriage within the community",
       color = "WWTP") +
  scale_color_manual(
    values = c("#fc8d62ff", "#8da0cbff", "#66c2a5ff","#e78ac3ff", "#ffd92fff", "#a6d854ff"),
    breaks = c(0.0157, 0.0129, 0.0207, 0.0169, 0.0189, 0.0126),  # Values of 'a'
    labels = c("ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Werdhölzli Zürich", "IDA CDA Lugano", "ARA Sensetal Laupen")  # Names for legend
  ) +
  scale_y_continuous(
    limits = c(0, 1),  # Set y-axis limits
    breaks = seq(0, 1, by = 0.2)  # Set y-axis tick positions
  ) +
  scale_x_log10(
    limits = c(0.01, 1),
    breaks = c(0.01,0.1, 1),  # Set x-axis tick positions
    labels = c(0.01, 0.1,1) # Set x-axis tick labels
  ) +
  theme_minimal()+
  geom_vline(data = data[data$Y == 1, ], aes(xintercept = a, color = factor(a)), linetype = "dashed", alpha=0.5) +
  theme(legend.position= c(0.9, 0.9),legend.justification = c(1, 1))

####ARA Altenrhein with CI-------
# Create a data frame with different values of "a", different one for each WWTP
data_alt <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                   a = c(0.0157))  

data_alt_min <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                       a = c(0.0127))

data_alt_max <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(0.0210))


# Calculate Ys based on the function Y = a/x
data_alt$Y <- data_alt$a / data_alt$x
data_alt$Ymin <- data_alt_min$a / data_alt_min$x
data_alt$Ymax <- data_alt_max$a /data_alt_max$x

# Set Y values greater than 1 to 1
data_alt$Y[data_alt$Y > 1] <- 1
data_alt$Ymin[data_alt$Ymin > 1] <- 1
data_alt$Ymax[data_alt$Ymax > 1] <- 1

# Create the plot
a <- ggplot(data_alt[data_alt$x >= 0.0157,], aes(x = x, y = Y)) +
  geom_line(color="#fc8d62ff") +
  geom_ribbon(data=data_alt, 
              aes(ymin = Ymin, ymax = Ymax), fill = "#fc8d62ff", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c(0.005, 0.01,0.1, 1),  
    labels = c(0.005, 0.01,0.1, 1)) +  
  labs(x = "Proportion of ESBL-Ec out of total E. coli in the gut",
       y = "Prevalence of ESBL-Ec carriage within the community") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.0157", vjust = 1.3, hjust=-0.0,data = data_alt[data_alt$Y == 1, ], 
                 aes(xintercept = a), linetype = "dashed", color="#fc8d62ff", alpha = 0.8)+
  ggtitle("ARA Altenrhein")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  
a

####ARA Chur with CI-------
# Create a data frame with different values of "a", different one for each WWTP
data_chu <- data.frame(x = seq(0.005, 1, by = 0.0001),  
                       a = c(0.0129))  

data_chu_min <- data.frame(x = seq(0.005, 1, by = 0.0001),  
                           a = c(0.0085))

data_chu_max <- data.frame(x = seq(0.005, 1, by = 0.0001),  
                           a = c(0.0165))


# Calculate Ys based on the function Y = a/x
data_chu$Y <- data_chu$a / data_chu$x
data_chu$Ymin <- data_chu_min$a / data_chu_min$x
data_chu$Ymax <- data_chu_max$a /data_chu_max$x

# Set Y values greater than 1 to 1
data_chu$Y[data_chu$Y > 1] <- 1
data_chu$Ymin[data_chu$Ymin > 1] <- 1
data_chu$Ymax[data_chu$Ymax > 1] <- 1

# Create the plot
b <- ggplot(data_chu[data_chu$x >= 0.0129,], aes(x = x, y = Y)) +
  geom_line(color="#8da0cbff") +
  geom_ribbon(data = data_chu,aes(ymin = Ymin, ymax = Ymax), fill = "#8da0cbff", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c( 0.005, 0.01, 0.1, 1),  
    labels = c( 0.005, 0.01, 0.1, 1)) +  
  labs(x = "Proportion of ESBL-Ec out of total E. coli in the gut",
       y = "Prevalence of ESBL-Ec carriage within the community") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.0129", vjust = 1.3, hjust=-0.0,data = data_chu[data_chu$Y == 1, ], 
                 aes(xintercept = a), linetype = "dashed", color="#8da0cbff", alpha = 0.8)+
  ggtitle("ARA Chur")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
b


####STEP d'Aïre Genève with CI-------
# Create a data frame with different values of "a", different one for each WWTP
data_gen <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                       a = c(0.0207))  

data_gen_min <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(0.0161))

data_gen_max <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(0.0263))


# Calculate Ys based on the function Y = a/x
data_gen$Y <- data_gen$a / data_gen$x
data_gen$Ymin <- data_gen_min$a / data_gen_min$x
data_gen$Ymax <- data_gen_max$a /data_gen_max$x

# Set Y values greater than 1 to 1
data_gen$Y[data_gen$Y > 1] <- 1
data_gen$Ymin[data_gen$Ymin > 1] <- 1
data_gen$Ymax[data_gen$Ymax > 1] <- 1

# Create the plot
c <- ggplot(data_gen[data_gen$x >= 0.0207,], aes(x = x, y = Y)) +
  geom_line(color="#66c2a5ff") +
  geom_ribbon(data = data_gen,aes(ymin = Ymin, ymax = Ymax), fill = "#66c2a5ff", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c( 0.005, 0.01, 0.1, 1),  
    labels = c( 0.005, 0.01, 0.1, 1)) +  
  labs(x = "Proportion of ESBL-Ec out of total E. coli in the gut",
       y = "Prevalence of ESBL-Ec carriage within the community") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.0207", vjust = 1.3, hjust=-0.0,data = data_gen[data_gen$Y == 1, ], 
                 aes(xintercept = a), linetype = "dashed", color="#66c2a5ff", alpha = 0.8)+
  ggtitle("STEP d'Aïre Genève")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
c
####ARA Werdhölzli Zürich with CI-------
# Create a data frame with different values of "a", different one for each WWTP
data_zur <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                       a = c(0.0169))  

data_zur_min <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(0.0140))

data_zur_max <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(0.0215))


# Calculate Ys based on the function Y = a/x
data_zur$Y <- data_zur$a / data_zur$x
data_zur$Ymin <- data_zur_min$a / data_zur_min$x
data_zur$Ymax <- data_zur_max$a /data_zur_max$x

# Set Y values greater than 1 to 1
data_zur$Y[data_zur$Y > 1] <- 1
data_zur$Ymin[data_zur$Ymin > 1] <- 1
data_zur$Ymax[data_zur$Ymax > 1] <- 1

# Create the plot
d <- ggplot(data_zur[data_zur$x >= 0.0169,], aes(x = x, y = Y)) +
  geom_line(color="#e78ac3ff") +
  geom_ribbon(data = data_zur,aes(ymin = Ymin, ymax = Ymax), fill = "#e78ac3ff", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c( 0.005, 0.01, 0.1, 1),  
    labels = c( 0.005, 0.01, 0.1, 1)) +  
  labs(x = "Proportion of ESBL-Ec out of total E. coli in the gut",
       y = "Prevalence of ESBL-Ec carriage within the community") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.0169", vjust = 1.3, hjust=-0.0,data = data_zur[data_zur$Y == 1, ], 
                 aes(xintercept = a), linetype = "dashed", color="#e78ac3ff", alpha = 0.8)+
  ggtitle("ARA Werdhölzli Zürich")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
d

####IDA CDA Lugano with CI-------
# Create a data frame with different values of "a", different one for each WWTP
data_lug <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                       a = c(0.0189))  

data_lug_min <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(0.0146))

data_lug_max <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(0.0237))


# Calculate Ys based on the function Y = a/x
data_lug$Y <- data_lug$a / data_lug$x
data_lug$Ymin <- data_lug_min$a / data_lug_min$x
data_lug$Ymax <- data_lug_max$a /data_lug_max$x

# Set Y values greater than 1 to 1
data_lug$Y[data_lug$Y > 1] <- 1
data_lug$Ymin[data_lug$Ymin > 1] <- 1
data_lug$Ymax[data_lug$Ymax > 1] <- 1

# Create the plot
e <- ggplot(data_lug[data_lug$x >= 0.0189,], aes(x = x, y = Y)) +
  geom_line(color="#ffd92fff") +
  geom_ribbon(data = data_lug,aes(ymin = Ymin, ymax = Ymax), fill = "#ffd92fff", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c( 0.005, 0.01, 0.1, 1),  
    labels = c( 0.005, 0.01, 0.1, 1)) +  
  labs(x = "Proportion of ESBL-Ec out of total E. coli in the gut",
       y = "Prevalence of ESBL-Ec carriage within the community") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.0189", vjust = 1.3, hjust=-0.0,data = data_lug[data_lug$Y == 1, ], 
                 aes(xintercept = a), linetype = "dashed", color="#ffd92fff", alpha = 0.8)+
  ggtitle("IDA CDA Lugano")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
e
####ARA Sensetal Laupen with CI-------
# Create a data frame with different values of "a", different one for each WWTP
data_sen <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                       a = c(0.0126))  

data_sen_min <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(0.0099))

data_sen_max <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(0.0166))


# Calculate Ys based on the function Y = a/x
data_sen$Y <- data_sen$a / data_sen$x
data_sen$Ymin <- data_sen_min$a / data_sen_min$x
data_sen$Ymax <- data_sen_max$a /data_sen_max$x

# Set Y values greater than 1 to 1
data_sen$Y[data_sen$Y > 1] <- 1
data_sen$Ymin[data_sen$Ymin > 1] <- 1
data_sen$Ymax[data_sen$Ymax > 1] <- 1

# Create the plot
f <- ggplot(data_sen[data_sen$x >= 0.0126,], aes(x = x, y = Y)) +
  geom_line(color="#a6d854ff") +
  geom_ribbon(data = data_sen,aes(ymin = Ymin, ymax = Ymax), fill = "#a6d854ff", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c( 0.005, 0.01, 0.1, 1),  
    labels = c( 0.005, 0.01, 0.1, 1)) +  
  labs(x = "Proportion of ESBL-Ec out of total E. coli in the gut",
       y = "Prevalence of ESBL-Ec carriage within the community") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.0126", vjust = 1.3, hjust=-0.0,data = data_sen[data_sen$Y == 1, ], 
             aes(xintercept = a), linetype = "dashed", color="#a6d854ff", alpha = 0.7)+
  ggtitle("ARA Sensetal Laupen")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

f

ggarrange(a,b, c, d, e, f,nrow = 2, ncol = 3, labels = c("A", "B", "C", "D","E", "F"))