#ESBL-Ec in WW------
##Load required packages------
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
library(MASS)
library(coin)
library(dunn.test)

##Formatting df-------------------
#Import datasets 
setwd("/Users/conforsh/switchdrive/Institution/Manuscripts/esblec-monitor-ww-21-22/data")
df <- read.csv("ecoli_counts.csv", header = TRUE)
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

#Create a df without NAs for the variable ESBL-E. coli percentage
df_filtered = filter(df, !is.na(average_ESBL_Ec)) 
#Subset the new filtered df
df_per_alt = df_filtered %>% filter(wwtp == "ARA Altenrhein")
df_per_chu = df_filtered %>% filter(wwtp == "ARA Chur")
df_per_gen = df_filtered %>% filter(wwtp == "STEP d'Aïre Genève")
df_per_zur = df_filtered %>% filter(wwtp == "ARA Werdhölzli Zürich")
df_per_lug = df_filtered %>% filter(wwtp == "IDA CDA Lugano")
df_per_sen = df_filtered %>% filter(wwtp == "ARA Sensetal Laupen")

#Create a df without NAs for the variable total-E. coli loads
df_filtered_t = filter(df, !is.na(average_loads_tot_Ec)) 
#Subset the new filtered df
df_tot_alt = df_filtered_t %>% filter(wwtp == "ARA Altenrhein")
df_tot_chu = df_filtered_t %>% filter(wwtp == "ARA Chur")
df_tot_gen = df_filtered_t %>% filter(wwtp == "STEP d'Aïre Genève")
df_tot_zur = df_filtered_t %>% filter(wwtp == "ARA Werdhölzli Z?rich")
df_tot_lug = df_filtered_t %>% filter(wwtp == "IDA CDA Lugano")
df_tot_sen = df_filtered_t %>% filter(wwtp == "ARA Sensetal Laupen")

#Create a df without NAs for the variable ESBL-E. coli loads
df_filtered_r = filter(df, !is.na(average_loads_ESBL_Ec)) 
#Subset the new filtered df
df_ESBL_alt = df_filtered_r %>% filter(wwtp == "ARA Altenrhein")
df_ESBL_chu = df_filtered_r %>% filter(wwtp == "ARA Chur")
df_ESBL_gen = df_filtered_r %>% filter(wwtp == "STEP d'Aïre Genève")
df_ESBL_zur = df_filtered_r %>% filter(wwtp == "ARA Werdhölzli Zürich")
df_ESBL_lug = df_filtered_r %>% filter(wwtp == "IDA CDA Lugano")
df_ESBL_sen = df_filtered_r %>% filter(wwtp == "ARA Sensetal Laupen")

##Calculate overall statistics------------------------
# Define a function to calculate summary statistics
calculate_stats <- function(data, log_transform = FALSE, exclude_mean_sd = FALSE) {
  if (log_transform) {
    data <- log10(data)
  }
  
  mean_value <- ifelse(!exclude_mean_sd, mean(data, na.rm = TRUE), NA)
  sd_value <- ifelse(!exclude_mean_sd, sd(data, na.rm = TRUE), NA)
  median_value <- quantile(data, probs = 0.5, na.rm = TRUE)
  quantile_values <- quantile(data, probs = c(0.25, 0.75), na.rm = TRUE)
  min_value <- min(data, na.rm = TRUE)
  max_value <- max(data, na.rm = TRUE)
  
  return(list(mean = mean_value, sd = sd_value, median = median_value, quantiles = quantile_values, min = min_value, max = max_value))
}

# Define a function to format numbers in scientific notation
format_scientific <- function(x) {
  formatted_values <- sapply(x, function(val) {
    power <- floor(log10(abs(val)))
    mantissa <- round(val / 10^power, 2)
    if (power == 0) {
      return(paste(mantissa, "x 10^0"))
    } else {
      return(paste(mantissa, "x 10^", power))
    }
  })
  return(formatted_values)
}

# List of data frames
data_frames <- list("Switzerland" = df, "ARA Altenrhein" = df_alt, "ARA Chur" = df_chu,"STEP d'Aire Genève" = df_gen,"ARA Werdhölzli Zürich" = df_zur,
                    "IDA CDA Lugano" = df_lug,"ARA Sensetal Laupen" = df_sen)

# Define the variables for which you want to calculate summary statistics
variables <- c("average_ESBL_Ec", "average_loads_tot_Ec", "average_loads_ESBL_Ec", "average_counts_tot_Ec", "average_counts_ESBL_Ec")

# Define custom names for the variables in the output
variable_names <- c("ESBL-Ec Percentage", "Total-Ec Loads", "ESBL-Ec Loads", "Total-Ec Counts", "ESBL-Ec Counts")

# Variables for which to exclude mean and sd
exclude_mean_sd <- c("average_ESBL_Ec", "average_counts_tot_Ec", "average_counts_ESBL_Ec")

# Loop through data frames and variables to calculate and display summary statistics
for (data_frame_name in names(data_frames)) {
  data_frame <- data_frames[[data_frame_name]]
  cat("Data Frame:", data_frame_name, "\n")
  
  for (i in 1:length(variables)) {
    variable_name <- variables[i]
    output_name <- variable_names[i]
    
    log_transform <- variable_name %in% c("average_loads_tot_Ec", "average_loads_ESBL_Ec")
    exclude_mean_sd_flag <- variable_name %in% exclude_mean_sd
    
    stats <- calculate_stats(data_frame[[variable_name]], log_transform, exclude_mean_sd_flag)
    formatted_min <- format_scientific(stats$min)
    formatted_max <- format_scientific(stats$max)
    formatted_median <- format_scientific(stats$median)
    formatted_quantiles <- format_scientific(stats$quantiles)
    
    cat(output_name, "- Mean:", ifelse(!exclude_mean_sd_flag, stats$mean, "NA"), 
        "- SD:", ifelse(!exclude_mean_sd_flag, stats$sd, "NA"), "- Median:", formatted_median, 
        "- 25th Quantile:", formatted_quantiles[1], "- 75th Quantile:", formatted_quantiles[2], 
        "- Min:", formatted_min, "- Max:", formatted_max, "\n")
  }
  cat("\n")
}

##Differences between WWTPs------
###Plots--------------------
####Percentage ESBL-E. coli-----------
#Order WWTPs in descending median percentage of ESBL-E. coli
m_o = c("STEP d'Aïre Genève", "IDA CDA Lugano","ARA Werdhölzli Zürich","ARA Altenrhein", "ARA Chur", "ARA Sensetal Laupen") #ordering WWTPs from highest to lowest percentage.
df$wwtp <- factor(df$wwtp, levels = m_o)

#Generate boxplot
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

# calculate CDF
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

# Generate plot with CDF curves 
custom_palette <- c("#fc8d62ff", "#8da0cbff", "#e78ac3ff", "#ffd92fff", "#a6d854ff", "#66c2a5ff") #line colours CDF plots

cdfper=ggplot(df_CDF_long, aes(x = x, y = CDF, color = wwtp)) +
  geom_line(linewidth=1.2) +
  theme_minimal()+
  theme(legend.position= c(0.99, 0.9),legend.justification = c(1, 1))+
  xlab(expression(paste("Percentage of ESBL- ", italic("E. coli"), " (%)"))) +
  ylab("Cumulative Distribution Function") +
  theme(axis.title.y = element_text(size=8.5), axis.title.x = element_text(size=9), axis.text.y = element_text(colour = "black"), axis.text.x = element_text(colour = "black"))+
  scale_color_manual(values = custom_palette)+
  labs(color = "WWTP")

####Loads ESBL-E. coli-----------
#Order WWTPs in descending median ESBL-E. coli loads
m_t = c("ARA Altenrhein","ARA Werdhölzli Zürich","STEP d'Aïre Genève", "ARA Sensetal Laupen", "ARA Chur", "IDA CDA Lugano") #ordering WWTPs from highest to lowest ESBL-Ec loads

#Generate boxplot
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

# calculate CDF
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

# Generate plot with CDF curves
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

####Loads total E. coli-----------
#Order WWTPs in descending median total E. coli loads
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

# calculate CDF
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

# Generate plot with CDF curves
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

#Arrange all plots together (boxplots + CDF curves)
ggarrange(perpl, cdfper, ltopl, cdftot, lrepl, cdfesbl,nrow=3, ncol=2, labels=c("A", "B", "C", "D", "E", "F"), common.legend = TRUE)

###Statistical tests--------------------------------
####Percentage ESBL-E. coli-----------
  #shapiro.test(df_filtered$average_ESBL_Ec)
  #Shapiro-Wilk normality test
  #W = 0.92056, p-value = 1.566e-11 
      #shapiro test indicates that the "average_ESBL_Ec" data are not normally distributed

# Perform the Kruskal-Wallis test
df_filtered$wwtp=as.factor(df_filtered$wwtp)
kruskal_test(average_ESBL_Ec ~ wwtp, data = df_filtered)

# perform pairwise comparisons using post-hoc test Dunn's test with "Bonferroni" correction to identify which groups differ significantly from each other.
dunn_test_result_perc_ESBLEc <- dunn.test(df_filtered$average_ESBL_Ec, df_filtered$wwtp, 
                              method = "bonferroni")
dunn_test_result_perc_ESBLEc

####Loads ESBL-E. coli-----------
  #shapiro.test(df_filtered_t$average_loads_tot_Ec)
  #Shapiro-Wilk normality test
  #data:  df_filtered_t$average_loads_tot_Ec
  #W = 0.95418, p-value = 4.252e-08
    #shapiro test indicates that the "average_loads_tot_Ec" data are not normally distributed

# Perform the Kruskal-Wallis test
df_filtered_t$wwtp=as.factor(df_filtered_t$wwtp)
kruskal_test(average_loads_tot_Ec ~ wwtp, data = df_filtered_t)

# perform pairwise comparisons using post-hoc test Dunn's test with "Bonferroni" correction to identify which groups differ significantly from each other.
dunn_test_result_loads_totEc <- dunn.test(df_filtered_t$average_loads_tot_Ec, df_filtered_t$wwtp, 
                              method = "bonferroni")
dunn_test_result_loads_totEc

####Loads ESBL E. coli-----------
    #shapiro.test(df_filtered_r$average_loads_ESBL_Ec)
    #Shapiro-Wilk normality test
    #data:  df_filtered_r$average_loads_ESBL_Ec
    #W = 0.95078, p-value = 1.584e-08
      #shapiro test indicates that the "average_loads_ESBL_Ec" data are not normally distributed

# Perform the Kruskal-Wallis test
df_filtered_r$wwtp=as.factor(df_filtered_r$wwtp)
kruskal_test(average_loads_ESBL_Ec ~ wwtp, data = df_filtered_r)

# perform pairwise comparisons using post-hoc test Dunn's test with "Bonferroni" correction to identify which groups differ significantly from each other.
dunn_test_result_loads_ESBLEc <- dunn.test(df_filtered_r$average_loads_ESBL_Ec, df_filtered_r$wwtp, 
                              method = "bonferroni")
dunn_test_result_loads_ESBLEc

##Differences between months-------------------------------
###Plots----------------------
####Percentage ESBL-E.coli----------------
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

####Loads ESBL-E.coli----------------
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

####Loads total E.coli----------------
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

###Statistical tests--------------------
####Percentage ESBL-E.coli----------------
#Perform Kruskal-wallis with Dunn's test and Bonferroni adjustmen as a loop on all the WWTP df to check for significant differences on a monthly scale
df_list <- list(df_per_alt, df_per_chu, df_per_gen, df_per_sen, df_per_lug, df_per_zur) #Create a list of data frames for each wwtp

wwtp_names <- c("ARA Altenrhein", "ARA Chur", "STEP d'Aïre Genève", "ARA Sensetal Laupen", 
                "IDA CDA Lugano", "ARA Werdhölzli Zürich") #Define the names of the wwtps

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

####Loads total E.coli----------------
#Perform Kruskal-wallis with Dunn's test and Bonferroni adjustmen as a loop on all the WWTP df to check for significant differences on a monthly scale
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

####Loads ESBL-E.coli----------------
#Perform Kruskal-wallis with Dunn's test and Bonferroni adjustmen as a loop on all the WWTP df to check for significant differences on a monthly scale
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

##Environmental variables-------
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

####Loads ESBL-E. coli--------------
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

####Loads total E. coli------------------------
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


##Estimate sample size---------------------------
#Bayesian code to fit log normal distribution to my data
graphics.off()
setwd("/Users/conforsh/switchdrive/Institution/Manuscripts/esblec-monitor-ww-21-22/code")
source("DBDA2E-utilities.R")
fileNameRoot="A" 
library(rjags)
library(runjags)

###Altenrhein----------------
df_per_alt = df %>% filter(wwtp == "ARA Altenrhein")
C = as.numeric(df_per_alt[, "average_ESBL_Ec"])
N = length(C)

#### Log Normal - Model ---------------------------
# Package the data for shipping to JAGS:
dataList = list(
  C = C ,
  N = N
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

#Likelihood LN 

C[i] ~ dlnorm(muOfLog,1/sigmaOfLog^2) 
}
# Prior 
muOfLog ~ dunif( -100 , 100)
sigmaOfLog ~ dexp(0.1)

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# RUN THE CHAINS
require(rjags)
parameters = c("muOfLog","sigmaOfLog") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
perc_ESBL_alt_ln = jags.model( "model.txt" , data=dataList ,
                               n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( perc_ESBL_alt_ln , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples( perc_ESBL_alt_ln , variable.names=parameters ,
                         n.iter=nPerChain , thin=thinSteps )

##### Examine results-DIC------------------------------------------------------------------------------
# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

median(mcmcChain[,"muOfLog"])
median(mcmcChain[,"sigmaOfLog"])

dic.samples(perc_ESBL_alt_ln, 10000)
##### Plot CCDF ------------------------------------------------------------------------------
tiff(file = paste("Alt_ln.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1,1.0E+0, 1.0E+1)
par(mar=c(4.1, 4.1, 1.1, 1.1))
par(cex.lab = 0.8)

plot(point_x, point_y, col="white", pch = 19, main = expression(paste("Altenrhein LogNormal")), xlab = expression(paste("ESBL-", italic("E. coli"), " (%)")), ylab =  expression(paste("Fraction of samples above the percentage")), xaxt="n", xlim = range(myTicks),ylim=c(0.001, 1),
     font.lab = 1, log='xy')
axis(side = 1, at = myTicks)

#Plot uncertainty log-normal distribution
#Simulating uncertainty interval:
HDI_output = matrix(1:160 , nrow = 80, ncol = 2)
pb = txtProgressBar(min = 0, max = 80, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(0.1, 80, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = 1-plnorm(10 ^ (j / 4), meanlog = mcmcChain[i,"muOfLog"], sdlog = mcmcChain[i,"sigmaOfLog"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 20, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot log-normal distribution
lnorm_data=rlnorm(n=500000, meanlog=median(mcmcChain[,"muOfLog"]), sdlog=median(mcmcChain[,"sigmaOfLog"]))
lnorm_x=sort(lnorm_data)
lnorm_y = 1 - ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="blue", lwd=1,lty=1)

point_x = sort(C)
point_y = 1-ecdf(C)(sort(C))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

dev.off()

#### Gamma - Model ---------------------------
# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {
   
    #Likelihood Gamma distribution
    
    C[i] ~ dgamma(shape,rate)
    
  }
  # Prior
  shape ~ dgamma(0.01, 0.01)
  rate ~ dgamma(0.01, 0.01)
}  
  
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# RUN THE CHAINS
require(rjags)
parameters = c("shape","rate") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
perc_ESBL_alt_gamma = jags.model( "model.txt" , data=dataList ,
                                  n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update(perc_ESBL_alt_gamma, n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples(perc_ESBL_alt_gamma, variable.names=parameters ,
                        n.iter=nPerChain , thin=thinSteps )

##### Examine results-DIC-----------------------------------------------------------------------------
# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

median(mcmcChain[,"shape"])
median(mcmcChain[,"rate"])

dic.samples(perc_ESBL_alt_gamma, 10000)
##### Plot CCDF ------------------------------------------------------------------------------
tiff(file = paste("Alt_gamma.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1,1.0E+0, 1.0E+1)
par(mar=c(4.1, 4.1, 1.1, 1.1))
par(cex.lab = 0.8)

plot(point_x, point_y, col="white", pch = 19, main = expression(paste("Altenrhein Gamma")), xlab = expression(paste("ESBL-", italic("E. coli"), " (%)")), ylab =  expression(paste("Fraction of samples above the percentage")), xaxt="n", xlim = range(myTicks),ylim=c(0.001, 1),
     font.lab = 1, log='xy')
axis(side = 1, at = myTicks)


#Plot uncertainty log-normal distribution 
#Simulating uncertainty interval:
HDI_output = matrix(1:160, nrow = 80, ncol = 2)
pb = txtProgressBar(min = 0, max = 80, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 80, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = 1-pgamma(10 ^ (j / 4), shape = mcmcChain[i,"shape"], rate = mcmcChain[i,"rate"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 20, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot gamma distribution
gamma_data=rgamma(n=500000, shape=median(mcmcChain[,"shape"]), rate=median(mcmcChain[,"rate"]))
gamma_x=sort(gamma_data)
gamma_y = 1 - ecdf(gamma_data)(sort(gamma_data) )
lines(gamma_x,gamma_y, col="blue", lwd=1,lty=1)

point_x = sort(C)
point_y = 1-ecdf(C)(sort(C))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

dev.off()

###Chur----------------
df_per_chu = df %>% filter(wwtp == "ARA Chur")
C = as.numeric(df_per_chu[, "average_ESBL_Ec"])
N = length(C)

####Log Normal - Model ---------------------------
# Package the data for shipping to JAGS:
dataList = list(
  C = C ,
  N = N
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

#Likelihood LN 

C[i] ~ dlnorm(muOfLog,1/sigmaOfLog^2) 
}
# Prior 
muOfLog ~ dunif( -100 , 100)
sigmaOfLog ~ dexp(0.1)

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# RUN THE CHAINS
require(rjags)
parameters = c("muOfLog","sigmaOfLog") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
perc_ESBL_chu_ln = jags.model( "model.txt" , data=dataList ,
                               n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( perc_ESBL_chu_ln , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples( perc_ESBL_chu_ln , variable.names=parameters ,
                         n.iter=nPerChain , thin=thinSteps )

##### Examine results-DIC------------------------------------------------------------------------------
# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

median(mcmcChain[,"muOfLog"])
median(mcmcChain[,"sigmaOfLog"])

dic.samples(perc_ESBL_chu_ln, 10000)
##### Plot CCDF ------------------------------------------------------------------------------
tiff(file = paste("Chu_ln.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1,1.0E+0, 1.0E+1)
par(mar=c(4.1, 4.1, 1.1, 1.1))
par(cex.lab = 0.8)

plot(point_x, point_y, col="white", pch = 19, main = expression(paste("Chur LogNormal")), xlab = expression(paste("ESBL-", italic("E. coli"), " (%)")), ylab =  expression(paste("Fraction of samples above the percentage")), xaxt="n", xlim = range(myTicks),ylim=c(0.001, 1),
     font.lab = 1, log='xy')
axis(side = 1, at = myTicks)

#Plot uncertainty log-normal distribution
#Simulating uncertainty interval:
HDI_output = matrix(1:160 , nrow = 80, ncol = 2)
pb = txtProgressBar(min = 0, max = 80, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 80, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = 1-plnorm(10 ^ (j / 4), meanlog = mcmcChain[i,"muOfLog"], sdlog = mcmcChain[i,"sigmaOfLog"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 20, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot log-normal distribution
lnorm_data=rlnorm(n=500000, meanlog=median(mcmcChain[,"muOfLog"]), sdlog=median(mcmcChain[,"sigmaOfLog"]))
lnorm_x=sort(lnorm_data)
lnorm_y = 1 - ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="blue", lwd=1,lty=1)

point_x = sort(C)
point_y = 1-ecdf(C)(sort(C))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

dev.off()

#### Gamma - Model ---------------------------
# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {
   
    #Likelihood Gamma distribution
    
    C[i] ~ dgamma(shape,rate)
    
  }
  # Prior
  shape ~ dgamma(0.01, 0.01)
  rate ~ dgamma(0.01, 0.01)
}  
  
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# RUN THE CHAINS
require(rjags)
parameters = c("shape","rate") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
perc_ESBL_chu_gamma = jags.model( "model.txt" , data=dataList ,
                                  n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update(perc_ESBL_chu_gamma, n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples(perc_ESBL_chu_gamma, variable.names=parameters ,
                        n.iter=nPerChain , thin=thinSteps )

##### Examine results-DIC-----------------------------------------------------------------------------
# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

median(mcmcChain[,"shape"])
median(mcmcChain[,"rate"])

dic.samples(perc_ESBL_chu_gamma, 10000)
##### Plot CCDF ------------------------------------------------------------------------------
tiff(file = paste("Chu_gamma.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1,1.0E+0, 1.0E+1)
par(mar=c(4.1, 4.1, 1.1, 1.1))
par(cex.lab = 0.8)

plot(point_x, point_y, col="white", pch = 19, main = expression(paste("Chur Gamma")), xlab = expression(paste("ESBL-", italic("E. coli"), " (%)")), ylab =  expression(paste("Fraction of samples above the percentage")), xaxt="n", xlim = range(myTicks),ylim=c(0.001, 1),
     font.lab = 1, log='xy')
axis(side = 1, at = myTicks)

#Plot uncertainty log-normal distribution 
#Simulating uncertainty interval:
HDI_output = matrix(1:160, nrow = 80, ncol = 2)
pb = txtProgressBar(min = 0, max = 80, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 80, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = 1-pgamma(10 ^ (j / 4), shape = mcmcChain[i,"shape"], rate = mcmcChain[i,"rate"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 20, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot gamma distribution
gamma_data=rgamma(n=500000, shape=median(mcmcChain[,"shape"]), rate=median(mcmcChain[,"rate"]))
gamma_x=sort(gamma_data)
gamma_y = 1 - ecdf(gamma_data)(sort(gamma_data) )
lines(gamma_x,gamma_y, col="blue", lwd=1,lty=1)

point_x = sort(C)
point_y = 1-ecdf(C)(sort(C))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

dev.off()

###Zürich----------------
df_per_zur = df %>% filter(wwtp == "ARA Werdhölzli Zürich")
C = as.numeric(df_per_zur[, "average_ESBL_Ec"])
N = length(C)

#### Log Normal - Model ---------------------------
# Package the data for shipping to JAGS:
dataList = list(
  C = C ,
  N = N
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

#Likelihood LN 

C[i] ~ dlnorm(muOfLog,1/sigmaOfLog^2) 
}
# Prior 
muOfLog ~ dunif( -100 , 100)
sigmaOfLog ~ dexp(0.1)

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# RUN THE CHAINS
require(rjags)
parameters = c("muOfLog","sigmaOfLog") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
perc_ESBL_zur_ln = jags.model( "model.txt" , data=dataList ,
                               n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( perc_ESBL_zur_ln , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples( perc_ESBL_zur_ln , variable.names=parameters ,
                         n.iter=nPerChain , thin=thinSteps )

##### Examine results-DIC------------------------------------------------------------------------------
# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

median(mcmcChain[,"muOfLog"])
median(mcmcChain[,"sigmaOfLog"])

dic.samples(perc_ESBL_zur_ln, 10000)
##### Plot CCDF ------------------------------------------------------------------------------
tiff(file = paste("Zur_ln.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1,1.0E+0, 1.0E+1)
par(mar=c(4.1, 4.1, 1.1, 1.1))
par(cex.lab = 0.8)

plot(point_x, point_y, col="white", pch = 19, main = expression(paste("Zürich LogNormal")), xlab = expression(paste("ESBL-", italic("E. coli"), " (%)")), ylab =  expression(paste("Fraction of samples above the percentage")), xaxt="n", xlim = range(myTicks),ylim=c(0.001, 1),
     font.lab = 1, log='xy')
axis(side = 1, at = myTicks)

#Plot uncertainty log-normal distribution
#Simulating uncertainty interval:
HDI_output = matrix(1:160 , nrow = 80, ncol = 2)
pb = txtProgressBar(min = 0, max = 80, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 80, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = 1-plnorm(10 ^ (j / 4), meanlog = mcmcChain[i,"muOfLog"], sdlog = mcmcChain[i,"sigmaOfLog"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 20, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot log-normal distribution
lnorm_data=rlnorm(n=500000, meanlog=median(mcmcChain[,"muOfLog"]), sdlog=median(mcmcChain[,"sigmaOfLog"]))
lnorm_x=sort(lnorm_data)
lnorm_y = 1 - ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="blue", lwd=1,lty=1)

point_x = sort(C)
point_y = 1-ecdf(C)(sort(C))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

dev.off()

#### Gamma - Model ---------------------------
# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {
   
    #Likelihood Gamma distribution
    
    C[i] ~ dgamma(shape,rate)
    
  }
  # Prior
  shape ~ dgamma(0.01, 0.01)
  rate ~ dgamma(0.01, 0.01)
}  
  
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# RUN THE CHAINS
require(rjags)
parameters = c("shape","rate") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
perc_ESBL_zur_gamma = jags.model( "model.txt" , data=dataList ,
                                  n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update(perc_ESBL_zur_gamma, n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples(perc_ESBL_zur_gamma, variable.names=parameters ,
                        n.iter=nPerChain , thin=thinSteps )

##### Examine results-DIC-----------------------------------------------------------------------------
# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

median(mcmcChain[,"shape"])
median(mcmcChain[,"rate"])

dic.samples(perc_ESBL_zur_gamma, 10000)
##### Plot CCDF ------------------------------------------------------------------------------
tiff(file = paste("Zur_gamma.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1,1.0E+0, 1.0E+1)
par(mar=c(4.1, 4.1, 1.1, 1.1))
par(cex.lab = 0.8)

plot(point_x, point_y, col="white", pch = 19, main = expression(paste("Zürich Gamma")), xlab = expression(paste("ESBL-", italic("E. coli"), " (%)")), ylab =  expression(paste("Fraction of samples above the percentage")), xaxt="n", xlim = range(myTicks),ylim=c(0.001, 1),
     font.lab = 1, log='xy')
axis(side = 1, at = myTicks)

#Plot uncertainty log-normal distribution 
#Simulating uncertainty interval:
HDI_output = matrix(1:160, nrow = 80, ncol = 2)
pb = txtProgressBar(min = 0, max = 80, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 80, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = 1-pgamma(10 ^ (j / 4), shape = mcmcChain[i,"shape"], rate = mcmcChain[i,"rate"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 20, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot gamma distribution
gamma_data=rgamma(n=500000, shape=median(mcmcChain[,"shape"]), rate=median(mcmcChain[,"rate"]))
gamma_x=sort(gamma_data)
gamma_y = 1 - ecdf(gamma_data)(sort(gamma_data) )
lines(gamma_x,gamma_y, col="blue", lwd=1,lty=1)

point_x = sort(C)
point_y = 1-ecdf(C)(sort(C))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

dev.off()

###Geneva----------------
df_per_gen = df %>% filter(wwtp == "STEP d'Aïre Genève")
C = as.numeric(df_per_gen[, "average_ESBL_Ec"])
N = length(C)

#### Log Normal - Model ---------------------------
# Package the data for shipping to JAGS:
dataList = list(
  C = C ,
  N = N
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

#Likelihood LN 

C[i] ~ dlnorm(muOfLog,1/sigmaOfLog^2) 
}
# Prior 
muOfLog ~ dunif( -100 , 100)
sigmaOfLog ~ dexp(0.1)

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# RUN THE CHAINS
require(rjags)
parameters = c("muOfLog","sigmaOfLog") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
perc_ESBL_gen_ln = jags.model( "model.txt" , data=dataList ,
                               n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( perc_ESBL_gen_ln , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples( perc_ESBL_gen_ln , variable.names=parameters ,
                         n.iter=nPerChain , thin=thinSteps )

##### Examine results-DIC------------------------------------------------------------------------------
# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

median(mcmcChain[,"muOfLog"])
median(mcmcChain[,"sigmaOfLog"])

dic.samples(perc_ESBL_gen_ln, 10000)
##### Plot CCDF ------------------------------------------------------------------------------
tiff(file = paste("Gen_ln.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1,1.0E+0, 1.0E+1)
par(mar=c(4.1, 4.1, 1.1, 1.1))
par(cex.lab = 0.8)

plot(point_x, point_y, col="white", pch = 19, main = expression(paste("Geneva LogNormal")), xlab = expression(paste("ESBL-", italic("E. coli"), " (%)")), ylab =  expression(paste("Fraction of samples above the percentage")), xaxt="n", xlim = range(myTicks),ylim=c(0.001, 1),
     font.lab = 1, log='xy')
axis(side = 1, at = myTicks)

#Plot uncertainty log-normal distribution
#Simulating uncertainty interval:
HDI_output = matrix(1:160 , nrow = 80, ncol = 2)
pb = txtProgressBar(min = 0, max = 80, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 80, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = 1-plnorm(10 ^ (j / 4), meanlog = mcmcChain[i,"muOfLog"], sdlog = mcmcChain[i,"sigmaOfLog"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 20, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot log-normal distribution
lnorm_data=rlnorm(n=500000, meanlog=median(mcmcChain[,"muOfLog"]), sdlog=median(mcmcChain[,"sigmaOfLog"]))
lnorm_x=sort(lnorm_data)
lnorm_y = 1 - ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="blue", lwd=1,lty=1)

point_x = sort(C)
point_y = 1-ecdf(C)(sort(C))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

dev.off()

#### Gamma - Model ---------------------------
# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {
   
    #Likelihood Gamma distribution
    
    C[i] ~ dgamma(shape,rate)
    
  }
  # Prior
  shape ~ dgamma(0.01, 0.01)
  rate ~ dgamma(0.01, 0.01)
}  
  
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# RUN THE CHAINS
require(rjags)
parameters = c("shape","rate") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
perc_ESBL_gen_gamma = jags.model( "model.txt" , data=dataList ,
                                  n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update(perc_ESBL_gen_gamma, n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples(perc_ESBL_gen_gamma, variable.names=parameters ,
                        n.iter=nPerChain , thin=thinSteps )

##### Examine results-DIC-----------------------------------------------------------------------------
# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

median(mcmcChain[,"shape"])
median(mcmcChain[,"rate"])

dic.samples(perc_ESBL_gen_gamma, 10000)

##### Plot CCDF ------------------------------------------------------------------------------
tiff(file = paste("Gen_gamma.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1,1.0E+0, 1.0E+1)
par(mar=c(4.1, 4.1, 1.1, 1.1))
par(cex.lab = 0.8)

plot(point_x, point_y, col="white", pch = 19, main = expression(paste("Geneva Gamma")), xlab = expression(paste("ESBL-", italic("E. coli"), " (%)")), ylab =  expression(paste("Fraction of samples above the percentage")), xaxt="n", xlim = range(myTicks),ylim=c(0.001, 1),
     font.lab = 1, log='xy')
axis(side = 1, at = myTicks)

#Plot uncertainty log-normal distribution 
#Simulating uncertainty interval:
HDI_output = matrix(1:160, nrow = 80, ncol = 2)
pb = txtProgressBar(min = 0, max = 80, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 80, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = 1-pgamma(10 ^ (j / 4), shape = mcmcChain[i,"shape"], rate = mcmcChain[i,"rate"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 20, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot gamma distribution
gamma_data=rgamma(n=500000, shape=median(mcmcChain[,"shape"]), rate=median(mcmcChain[,"rate"]))
gamma_x=sort(gamma_data)
gamma_y = 1 - ecdf(gamma_data)(sort(gamma_data) )
lines(gamma_x,gamma_y, col="blue", lwd=1,lty=1)

point_x = sort(C)
point_y = 1-ecdf(C)(sort(C))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

dev.off()

###Lugano----------------
df_per_lug = df %>% filter(wwtp == "IDA CDA Lugano")
C = as.numeric(df_per_lug[, "average_ESBL_Ec"])
N = length(C)

#### Log Normal - Model ---------------------------
# Package the data for shipping to JAGS:
dataList = list(
  C = C ,
  N = N
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

#Likelihood LN 

C[i] ~ dlnorm(muOfLog,1/sigmaOfLog^2) 
}
# Prior 
muOfLog ~ dunif( -100 , 100)
sigmaOfLog ~ dexp(0.1)

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# RUN THE CHAINS
require(rjags)
parameters = c("muOfLog","sigmaOfLog") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
perc_ESBL_lug_ln = jags.model( "model.txt" , data=dataList ,
                               n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( perc_ESBL_lug_ln , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples( perc_ESBL_lug_ln , variable.names=parameters ,
                         n.iter=nPerChain , thin=thinSteps )

##### Examine results-DIC------------------------------------------------------------------------------
# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

median(mcmcChain[,"muOfLog"])
median(mcmcChain[,"sigmaOfLog"])

dic.samples(perc_ESBL_lug_ln, 10000)
##### Plot CCDF ------------------------------------------------------------------------------
tiff(file = paste("Lug_ln.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1,1.0E+0, 1.0E+1)
par(mar=c(4.1, 4.1, 1.1, 1.1))
par(cex.lab = 0.8)

plot(point_x, point_y, col="white", pch = 19, main = expression(paste("Lugano LogNormal")), xlab = expression(paste("ESBL-", italic("E. coli"), " (%)")), ylab =  expression(paste("Fraction of samples above the percentage")), xaxt="n", xlim = range(myTicks),ylim=c(0.001, 1),
     font.lab = 1, log='xy')
axis(side = 1, at = myTicks)

#Plot uncertainty log-normal distribution
#Simulating uncertainty interval:
HDI_output = matrix(1:160 , nrow = 80, ncol = 2)
pb = txtProgressBar(min = 0, max = 80, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 80, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = 1-plnorm(10 ^ (j / 4), meanlog = mcmcChain[i,"muOfLog"], sdlog = mcmcChain[i,"sigmaOfLog"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 20, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot log-normal distribution
lnorm_data=rlnorm(n=500000, meanlog=median(mcmcChain[,"muOfLog"]), sdlog=median(mcmcChain[,"sigmaOfLog"]))
lnorm_x=sort(lnorm_data)
lnorm_y = 1 - ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="blue", lwd=1,lty=1)

point_x = sort(C)
point_y = 1-ecdf(C)(sort(C))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

dev.off()

#### Gamma - Model ---------------------------
# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {
   
    #Likelihood Gamma distribution
    
    C[i] ~ dgamma(shape,rate)
    
  }
  # Prior
  shape ~ dgamma(0.01, 0.01)
  rate ~ dgamma(0.01, 0.01)
}  
  
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# RUN THE CHAINS
require(rjags)
parameters = c("shape","rate") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
perc_ESBL_lug_gamma = jags.model( "model.txt" , data=dataList ,
                                  n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update(perc_ESBL_lug_gamma, n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples(perc_ESBL_lug_gamma, variable.names=parameters ,
                        n.iter=nPerChain , thin=thinSteps )

##### Examine results-DIC-----------------------------------------------------------------------------
# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

median(mcmcChain[,"shape"])
median(mcmcChain[,"rate"])

dic.samples(perc_ESBL_lug_gamma, 10000)
##### Plot CCDF ------------------------------------------------------------------------------
tiff(file = paste("Lug_gamma.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1,1.0E+0, 1.0E+1)
par(mar=c(4.1, 4.1, 1.1, 1.1))
par(cex.lab = 0.8)

plot(point_x, point_y, col="white", pch = 19, main = expression(paste("Lugano Gamma")), xlab = expression(paste("ESBL-", italic("E. coli"), " (%)")), ylab =  expression(paste("Fraction of samples above the percentage")), xaxt="n", xlim = range(myTicks),ylim=c(0.001, 1),
     font.lab = 1, log='xy')
axis(side = 1, at = myTicks)

#Plot uncertainty log-normal distribution 
#Simulating uncertainty interval:
HDI_output = matrix(1:160, nrow = 80, ncol = 2)
pb = txtProgressBar(min = 0, max = 80, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 80, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = 1-pgamma(10 ^ (j / 4), shape = mcmcChain[i,"shape"], rate = mcmcChain[i,"rate"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 20, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot gamma distribution
gamma_data=rgamma(n=500000, shape=median(mcmcChain[,"shape"]), rate=median(mcmcChain[,"rate"]))
gamma_x=sort(gamma_data)
gamma_y = 1 - ecdf(gamma_data)(sort(gamma_data) )
lines(gamma_x,gamma_y, col="blue", lwd=1,lty=1)

point_x = sort(C)
point_y = 1-ecdf(C)(sort(C))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

dev.off()

###Sensetal Laupen----------------
df_per_sen = df %>% filter(wwtp == "ARA Sensetal Laupen")
C = as.numeric(df_per_sen[, "average_ESBL_Ec"])
N = length(C)

#### Log Normal - Model ---------------------------
# Package the data for shipping to JAGS:
dataList = list(
  C = C ,
  N = N
)

# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {

#Likelihood LN 

C[i] ~ dlnorm(muOfLog,1/sigmaOfLog^2) 
}
# Prior 
muOfLog ~ dunif( -100 , 100)
sigmaOfLog ~ dexp(0.1)

}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# RUN THE CHAINS
require(rjags)
parameters = c("muOfLog","sigmaOfLog") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
perc_ESBL_sen_ln = jags.model( "model.txt" , data=dataList ,
                               n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( perc_ESBL_sen_ln , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples( perc_ESBL_sen_ln , variable.names=parameters ,
                         n.iter=nPerChain , thin=thinSteps )

##### Examine results-DIC------------------------------------------------------------------------------
# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

median(mcmcChain[,"muOfLog"])
median(mcmcChain[,"sigmaOfLog"])

dic.samples(perc_ESBL_sen_ln, 10000)
##### Plot CCDF ------------------------------------------------------------------------------
tiff(file = paste("Sen_ln.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1,1.0E+0, 1.0E+1)
par(mar=c(4.1, 4.1, 1.1, 1.1))
par(cex.lab = 0.8)

plot(point_x, point_y, col="white", pch = 19, main = expression(paste("Sensetal LogNormal")), xlab = expression(paste("ESBL-", italic("E. coli"), " (%)")), ylab =  expression(paste("Fraction of samples above the percentage")), xaxt="n", xlim = range(myTicks),ylim=c(0.001, 1),
     font.lab = 1, log='xy')
axis(side = 1, at = myTicks)

#Plot uncertainty log-normal distribution
#Simulating uncertainty interval:
HDI_output = matrix(1:160 , nrow = 80, ncol = 2)
pb = txtProgressBar(min = 0, max = 80, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 80, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = 1-plnorm(10 ^ (j / 4), meanlog = mcmcChain[i,"muOfLog"], sdlog = mcmcChain[i,"sigmaOfLog"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 20, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot log-normal distribution
lnorm_data=rlnorm(n=500000, meanlog=median(mcmcChain[,"muOfLog"]), sdlog=median(mcmcChain[,"sigmaOfLog"]))
lnorm_x=sort(lnorm_data)
lnorm_y = 1 - ecdf(lnorm_data)(sort(lnorm_data) )
lines(lnorm_x,lnorm_y, col="blue", lwd=1,lty=1)

point_x = sort(C)
point_y = 1-ecdf(C)(sort(C))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

dev.off()

#### Gamma - Model ---------------------------
# THE MODEL.
modelstring ="
model {
for( i in 1 : N ) {
   
    #Likelihood Gamma distribution
    
    C[i] ~ dgamma(shape,rate)
    
  }
  # Prior
  shape ~ dgamma(0.01, 0.01)
  rate ~ dgamma(0.01, 0.01)
}  
  
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

# RUN THE CHAINS
require(rjags)
parameters = c("shape","rate") 
adaptSteps = 5000         # Number of steps to "tune" the samplers.
burnInSteps = 5000        # Number of steps to "burn-in" the samplers.
nChains = 3               # Number of chains to run.
numSavedSteps=10000       # Total number of steps in chains to save.
thinSteps=1               # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
perc_ESBL_sen_gamma = jags.model( "model.txt" , data=dataList ,
                                  n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update(perc_ESBL_sen_gamma, n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
mcmcCoda = coda.samples(perc_ESBL_sen_gamma, variable.names=parameters ,
                        n.iter=nPerChain , thin=thinSteps )

##### Examine results-DIC-----------------------------------------------------------------------------
# Convert coda-object codaSamples to matrix object for easier handling.
mcmcChain = as.matrix( mcmcCoda )
chainLength = NROW(mcmcChain)

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName)
}

median(mcmcChain[,"shape"])
median(mcmcChain[,"rate"])

dic.samples(perc_ESBL_sen_gamma, 10000)
##### Plot CCDF ------------------------------------------------------------------------------
tiff(file = paste("Sen_gamma.tiff", sep = ""), width = 2000, height = 2000, 
     units = "px",res = 1000, pointsize = 7)

options(scipen = 0)
point_x = 1
point_y = 1
myTicks<-c(1.0E-1,1.0E+0, 1.0E+1)
par(mar=c(4.1, 4.1, 1.1, 1.1))
par(cex.lab = 0.8)

plot(point_x, point_y, col="white", pch = 19, main = expression(paste("Sensetal Gamma")), xlab = expression(paste("ESBL-", italic("E. coli"), " (%)")), ylab =  expression(paste("Fraction of samples above the percentage")), xaxt="n", xlim = range(myTicks),ylim=c(0.001, 1),
     font.lab = 1, log='xy')
axis(side = 1, at = myTicks)

#Plot uncertainty log-normal distribution 
#Simulating uncertainty interval:
HDI_output = matrix(1:160, nrow = 80, ncol = 2)
pb = txtProgressBar(min = 0, max = 80, initial = 0, style = 3)
cat("Simulating uncertainty interval...\n")

for (j in seq(1, 80, 1)) {
  output <- numeric(length = nrow(mcmcChain))
  for (i in 1:nrow(mcmcChain)) {            
    output[[i]] = 1-pgamma(10 ^ (j / 4), shape = mcmcChain[i,"shape"], rate = mcmcChain[i,"rate"])
  }
  HDI_output[j, ] = HDIofMCMC(output, 0.95)
  setTxtProgressBar(pb,j)
}

# Representation of the interval as a surface:
MatrixA = cbind(c(10 ^ (seq(0.25, 20, 0.25))), HDI_output)
MatrixB = MatrixA[(MatrixA[,2] > 0),]
MatrixC = MatrixA[(MatrixA[,3] > 0),]
polygon(c(MatrixC[,1], rev(MatrixB[,1])), c(MatrixC[,3], rev(MatrixB[,2])), col = rgb(0, 0, 1, 0.2),border = NA) 

#Plot gamma distribution
gamma_data=rgamma(n=500000, shape=median(mcmcChain[,"shape"]), rate=median(mcmcChain[,"rate"]))
gamma_x=sort(gamma_data)
gamma_y = 1 - ecdf(gamma_data)(sort(gamma_data) )
lines(gamma_x,gamma_y, col="blue", lwd=1,lty=1)

point_x = sort(C)
point_y = 1-ecdf(C)(sort(C))
points(point_x, point_y ,col = "black", bg = "red", pch = 21, lwd = 0.5, cex = 0.75)

dev.off()

###95% CI plotting-------------------------
#df formatting
df = filter(df, !is.na(average_ESBL_Ec))

#Create dataframes
sample_sizes = seq(1, 156, by = 1) #only for graphs up to 3x/week
true_mean = 0

df_alt = df %>% filter(wwtp == "ARA Altenrhein")
alpha_alt = mean(df_alt$average_ESBL_Ec)^2/sd(df_alt$average_ESBL_Ec)^2
beta_alt = mean(df_alt$average_ESBL_Ec)/sd(df_alt$average_ESBL_Ec)^2
gmean_alt=alpha_alt/beta_alt
se_alt <- sqrt(alpha_alt / beta_alt^2 / sample_sizes)
CI_alt = se_alt*1.96

data_alt <- data.frame(sample_sizes=sample_sizes, wwtp="ARA Altenrhein",
                       CI=CI_alt)

df_chu = df %>% filter(wwtp == "ARA Chur")
alpha_chu = mean(df_chu$average_ESBL_Ec)^2/sd(df_chu$average_ESBL_Ec)^2
beta_chu = mean(df_chu$average_ESBL_Ec)/sd(df_chu$average_ESBL_Ec)^2
gmean_chu=alpha_chu/beta_chu
se_chu <- sqrt(alpha_chu / beta_chu^2 / sample_sizes)
CI_chu = se_chu*1.96

data_chu <- data.frame(sample_sizes=sample_sizes, wwtp="ARA Chur",
                       CI=CI_chu)

df_zur = df %>% filter(wwtp == "ARA Werdhölzli Zürich")
alpha_zur = mean(df_zur$average_ESBL_Ec)^2/sd(df_zur$average_ESBL_Ec)^2
beta_zur = mean(df_zur$average_ESBL_Ec)/sd(df_zur$average_ESBL_Ec)^2
gmean_zur=alpha_zur/beta_zur
se_zur <- sqrt(alpha_zur / beta_zur^2 / sample_sizes)
CI_zur = se_zur*1.96

data_zur <- data.frame(sample_sizes=sample_sizes, wwtp="ARA Werdhölzli Zürich",
                       CI=CI_zur)

df_gen = df %>% filter(wwtp == "STEP d'Aïre Genève")
alpha_gen = mean(df_gen$average_ESBL_Ec)^2/sd(df_gen$average_ESBL_Ec)^2
beta_gen = mean(df_gen$average_ESBL_Ec)/sd(df_gen$average_ESBL_Ec)^2
gmean_gen=alpha_gen/beta_gen
se_gen <- sqrt(alpha_gen / beta_gen^2 / sample_sizes)
CI_gen = se_gen*1.96

data_gen <- data.frame(sample_sizes=sample_sizes, wwtp="STEP d'Aïre Genève",
                       CI=CI_gen)

df_sen = df %>% filter(wwtp == "ARA Sensetal Laupen")
alpha_sen = mean(df_sen$average_ESBL_Ec)^2/sd(df_sen$average_ESBL_Ec)^2
beta_sen = mean(df_sen$average_ESBL_Ec)/sd(df_sen$average_ESBL_Ec)^2
gmean_sen=alpha_sen/beta_sen
se_sen <- sqrt(alpha_sen / beta_sen^2 / sample_sizes)
CI_sen = se_sen*1.96

data_sen <- data.frame(sample_sizes=sample_sizes, wwtp="ARA Sensetal Laupen",
                       CI=CI_sen)

df_lug = df %>% filter(wwtp == "IDA CDA Lugano")
alpha_lug = mean(df_lug$average_ESBL_Ec)^2/sd(df_lug$average_ESBL_Ec)^2
beta_lug = mean(df_lug$average_ESBL_Ec)/sd(df_lug$average_ESBL_Ec)^2
gmean_lug=alpha_lug/beta_lug
se_lug <- sqrt(alpha_lug / beta_lug^2 / sample_sizes)
CI_lug = se_lug*1.96

data_lug <- data.frame(sample_sizes=sample_sizes, wwtp="IDA CDA Lugano",
                       CI=CI_lug)

#Bind different df
data <- rbind(data_alt, data_chu, data_zur, data_gen, data_sen, data_lug)

# Generate plot
tot=ggplot(data, aes(x = sample_sizes)) +
  geom_line(aes(y = CI, color=wwtp), position = position_dodge(width = 0), linewidth=1, alpha=0.8) +
  theme_minimal()+
  scale_color_manual(values = c("IDA CDA Lugano" = "#a6d854ff", "ARA Altenrhein" = "#fc8d62ff", "ARA Sensetal Laupen" = "#e78ac3ff",
                                "ARA Chur" = "#8da0cbff", "STEP d'Aïre Genève" = "#66c2a5ff", "ARA Werdhölzli Zürich" = "#ffd92fff")) +
  theme(axis.text.x=element_text(angle=90), legend.position= "right")+
  scale_x_continuous(breaks = c(1,12, 24, 52, 104, 156), labels = c("1x Year","1x Month", "2x Month", "1x Week", "2x Week", "3x Week"))+
  xlab("Sampling Frequency") +
  ylab("Width of 95% CI") +
  labs(color="WWTP")+
  geom_vline(xintercept = c(1, 12, 24, 52, 104, 156), linetype = "dotted")+
  ggtitle(expression(paste("Percentage of ESBL- ", italic("E. coli"), " (%)")))

print(tot)

##Carriage estimation------
###Define estimates-----------------
#Gamma_mean of ESBL-Ec percentage in WW (%) (Xi)
galt = 1.722433   
gchu = 1.397264   
ggen = 2.019172   
gzur = 1.913655   
glug = 1.970193   
gsen = 1.373625   

      #gamma mean calculated as: gamma_mean=alpha/beta , where 
      #alpha = mean(ESBL_Ec_percentage)^2/sd(ESBL_Ec_percentage)^2
      #beta = mean(ESBL_Ec_percentage)/sd(ESBL_Ec_percentage)^2
      #These lines are not part of the code, they only inform about the formula used, variables do not exist.

#95% CI of ESBL-E.coli percentage in each WWTP 
CIalt=0.1955347
CIchu=0.2167696
CIgen=0.1807549
CIzur=0.2522777
CIlug=0.2267069
CIsen=0.1521924
  
      #95 CI of gamma mean calculated as: 1.96*se, where 
      #se = sqrt(alpha/beta^2 /sample_sizes)
      #These lines are not part of the code, they only inform about the formula used, variables do not exist.
  
#Population connected to each WWTP (n° of inhabitants) (Pi)
palt=64000
pchu=55000
pgen=454000
pzur=471000
plug=124000
psen=62000

###Whole Switzerland, with population adjusted ESBL-Ec percentage in WW----------------
#S.E of population weighted mean
#(SEMw)^2 = [n/(n-1)(∑Pi)^2]*[∑(Pi*Xi - P*Xw)^2-2*Xw*∑(Pi-P)(Pi*Xi - P*Xw)+Xw^2*∑(Pi-p)^2]

calculate_SEMw <- function(Pi, Xi, Xw, P) {
  n <- length(Pi)
  
  term1 <- n / ((n - 1) * sum(Pi)^2)
  
  term2 <- (sum((Pi * Xi - P * Xw)^2) - 
                 2 * Xw * sum((Pi - P) * (Pi * Xi - P * Xw)) + 
                 Xw^2 * sum((Pi - P)^2))
  
  SEMw <- sqrt(term1 * term2)
  
  return(SEMw)
}

# Replace Pi, Xi, Xw, and P with actual gamma mean and population size
Pi <- c(palt, pchu, pgen, pzur, plug, psen)  # Population in each site
Xi <- c(galt, gchu, ggen, gzur, glug, gsen)  # Gamma mean ofESBL-Ec percentage in each site
Xw <- weighted.mean(Xi, Pi)  # Population-weighted mean over sites
P <- mean(Pi)  # Mean of the population across sites

SEMw <- calculate_SEMw(Pi, Xi, Xw, P)
print(SEMw)

CIch = SEMw*1.96
print(CIch)

# Create data frame for Switzerland
data_ch <- data.frame(x = rep(seq(0.01, 1, by = 0.0001)),  
                      a = rep(c(Xw/100))) #a = population weighted mean of gamma_mean of each WWTP.

data_ch_min <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                          a = c(Xw/100-CIch/100)) #a = population weighted mean - CI of weighted mean

data_ch_max <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                          a = c(Xw/100+CIch/100))

# Calculate Y based on the function Y = a/x
data_ch$Y <- data_ch$a / data_ch$x
data_ch$Ymin <- data_ch_min$a / data_ch_min$x
data_ch$Ymax <- data_ch_max$a /data_ch_max$x

# Set Y values greater than 1 to 1
data_ch$Y[data_ch$Y > 1] <- 1
data_ch$Ymin[data_ch$Ymin > 1] <- 1
data_ch$Ymax[data_ch$Ymax > 1] <- 1

#Plot
ch <- ggplot(data_ch[data_ch$x >= Xw/100,], aes(x = x, y = Y)) +
  geom_line(color="black", linewidth=1) +
  geom_ribbon(data=data_ch, 
  aes(ymin = Ymin, ymax = Ymax), fill = "black", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c(0.005, 0.01,0.1, 1),  
    labels = c(0.005, 0.01,0.1, 1)) +  
  labs(x = expression(paste("Percentage of ESBL-", italic("E. coli"), " out of total ", italic("E. coli"), " in the gut")),
       y = expression(paste("Prevalence of ESBL-", italic("E. coli"), " carriage within the community"))) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.019", vjust = 1.3, hjust=-0.0,data = data_ch[data_ch$Y == 1, ], 
                 aes(xintercept = a), linetype = "dashed", color="red", alpha = 0.8)+
  ggtitle("Switzerland")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

print(ch)

###ARA Altenrhein with CI-------
# Create a data frame with different values of "a": a=gamma_mean, a=gamma_mean - CI, a=gamma_mean + CI
data_alt <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                   a = c(galt/100))  #a = gamma_mean

data_alt_min <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                       a = c(galt/100-CIalt/100)) #a = gamma_mean - CI

data_alt_max <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(galt/100+CIalt/100)) #a = gamma_mean + CI


# Calculate Ys based on the function Y = a/x
data_alt$Y <- data_alt$a / data_alt$x
data_alt$Ymin <- data_alt_min$a / data_alt_min$x
data_alt$Ymax <- data_alt_max$a /data_alt_max$x

# Set Y values greater than 1 to 1
data_alt$Y[data_alt$Y > 1] <- 1
data_alt$Ymin[data_alt$Ymin > 1] <- 1
data_alt$Ymax[data_alt$Ymax > 1] <- 1

# Create the plot
a <- ggplot(data_alt[data_alt$x >= galt/100,], aes(x = x, y = Y)) +
  geom_line(color="#fc8d62ff") +
  geom_ribbon(data=data_alt, 
              aes(ymin = Ymin, ymax = Ymax), fill = "#fc8d62ff", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c(0.005, 0.01,0.1, 1),  
    labels = c(0.005, 0.01,0.1, 1)) +  
  labs(x = expression(paste("Percentage of ESBL-", italic("E. coli"), " out of total ", italic("E. coli"), " in the gut")),
       y = expression(paste("Prevalence of ESBL-", italic("E. coli"), " carriage within the community"))) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.017", vjust = 1.3, hjust=-0.0,data = data_alt[data_alt$Y == 1, ], 
                 aes(xintercept = a), linetype = "dashed", color="#fc8d62ff", alpha = 0.8)+
  ggtitle("ARA Altenrhein")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  
a

###ARA Chur with CI-------
# Create a data frame with different values of "a": a=gamma_mean, a=gamma_mean - CI, a=gamma_mean + CI
data_chu <- data.frame(x = seq(0.005, 1, by = 0.0001),  
                       a = c(gchu/100))  #a = gamma_mean

data_chu_min <- data.frame(x = seq(0.005, 1, by = 0.0001),  
                           a = c(gchu/100-CIchu/100)) #a = gamma_mean - CI

data_chu_max <- data.frame(x = seq(0.005, 1, by = 0.0001),  
                           a = c(gchu/100+CIchu/100)) #a = gamma_mean + CI


# Calculate Ys based on the function Y = a/x
data_chu$Y <- data_chu$a / data_chu$x
data_chu$Ymin <- data_chu_min$a / data_chu_min$x
data_chu$Ymax <- data_chu_max$a /data_chu_max$x

# Set Y values greater than 1 to 1
data_chu$Y[data_chu$Y > 1] <- 1
data_chu$Ymin[data_chu$Ymin > 1] <- 1
data_chu$Ymax[data_chu$Ymax > 1] <- 1

# Create the plot
b <- ggplot(data_chu[data_chu$x >= gchu/100,], aes(x = x, y = Y)) +
  geom_line(color="#8da0cbff") +
  geom_ribbon(data = data_chu,aes(ymin = Ymin, ymax = Ymax), fill = "#8da0cbff", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c( 0.005, 0.01, 0.1, 1),  
    labels = c( 0.005, 0.01, 0.1, 1)) +  
  labs(x = expression(paste("Percentage of ESBL-", italic("E. coli"), " out of total ", italic("E. coli"), " in the gut")),
       y = expression(paste("Prevalence of ESBL-", italic("E. coli"), " carriage within the community"))) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.014", vjust = 1.3, hjust=-0.0,data = data_chu[data_chu$Y == 1, ], 
                 aes(xintercept = a), linetype = "dashed", color="#8da0cbff", alpha = 0.8)+
  ggtitle("ARA Chur")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
b


###STEP d'Aïre Genève with CI-------
# Create a data frame with different values of "a": a=gamma_mean, a=gamma_mean - CI, a=gamma_mean + CI
data_gen <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                       a = c(ggen/100))  #a = gamma_mean

data_gen_min <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(ggen/100-CIgen/100)) #a = gamma_mean - CI

data_gen_max <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(ggen/100+CIgen/100)) #a = gamma_mean + CI


# Calculate Ys based on the function Y = a/x
data_gen$Y <- data_gen$a / data_gen$x
data_gen$Ymin <- data_gen_min$a / data_gen_min$x
data_gen$Ymax <- data_gen_max$a /data_gen_max$x

# Set Y values greater than 1 to 1
data_gen$Y[data_gen$Y > 1] <- 1
data_gen$Ymin[data_gen$Ymin > 1] <- 1
data_gen$Ymax[data_gen$Ymax > 1] <- 1

# Create the plot
c <- ggplot(data_gen[data_gen$x >= ggen/100,], aes(x = x, y = Y)) +
  geom_line(color="#66c2a5ff") +
  geom_ribbon(data = data_gen,aes(ymin = Ymin, ymax = Ymax), fill = "#66c2a5ff", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c( 0.005, 0.01, 0.1, 1),  
    labels = c( 0.005, 0.01, 0.1, 1)) +  
  labs(x = expression(paste("Percentage of ESBL-", italic("E. coli"), " out of total ", italic("E. coli"), " in the gut")),
       y = expression(paste("Prevalence of ESBL-", italic("E. coli"), " carriage within the community"))) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.020", vjust = 1.3, hjust=-0.0,data = data_gen[data_gen$Y == 1, ], 
                 aes(xintercept = a), linetype = "dashed", color="#66c2a5ff", alpha = 0.8)+
  ggtitle("STEP d'Aïre Genève")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
c
###ARA Werdhölzli Zürich with CI-------
# Create a data frame with different values of "a": a=gamma_mean, a=gamma_mean - CI, a=gamma_mean + CI
data_zur <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                       a = c(gzur/100))  #a = gamma_mean

data_zur_min <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(gzur/100-CIzur/100)) #a = gamma_mean - CI

data_zur_max <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(gzur/100+CIzur/100)) #a = gamma_mean + CI


# Calculate Ys based on the function Y = a/x
data_zur$Y <- data_zur$a / data_zur$x
data_zur$Ymin <- data_zur_min$a / data_zur_min$x
data_zur$Ymax <- data_zur_max$a /data_zur_max$x

# Set Y values greater than 1 to 1
data_zur$Y[data_zur$Y > 1] <- 1
data_zur$Ymin[data_zur$Ymin > 1] <- 1
data_zur$Ymax[data_zur$Ymax > 1] <- 1

# Create the plot
d <- ggplot(data_zur[data_zur$x >= gzur/100,], aes(x = x, y = Y)) +
  geom_line(color="#e78ac3ff") +
  geom_ribbon(data = data_zur,aes(ymin = Ymin, ymax = Ymax), fill = "#e78ac3ff", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c( 0.005, 0.01, 0.1, 1),  
    labels = c( 0.005, 0.01, 0.1, 1)) +  
  labs(x = expression(paste("Percentage of ESBL-", italic("E. coli"), " out of total ", italic("E. coli"), " in the gut")),
       y = expression(paste("Prevalence of ESBL-", italic("E. coli"), " carriage within the community"))) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.019", vjust = 1.3, hjust=-0.0,data = data_zur[data_zur$Y == 1, ], 
                 aes(xintercept = a), linetype = "dashed", color="#e78ac3ff", alpha = 0.8)+
  ggtitle("ARA Werdhölzli Zürich")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
d

###IDA CDA Lugano with CI-------
# Create a data frame with different values of "a": a=gamma_mean, a=gamma_mean - CI, a=gamma_mean + CI
data_lug <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                       a = c(glug/100))  #a = gamma_mean

data_lug_min <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(glug/100-CIlug/100)) #a = gamma_mean - CI

data_lug_max <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(glug/100+CIlug/100)) #a = gamma_mean + CI


# Calculate Ys based on the function Y = a/x
data_lug$Y <- data_lug$a / data_lug$x
data_lug$Ymin <- data_lug_min$a / data_lug_min$x
data_lug$Ymax <- data_lug_max$a /data_lug_max$x

# Set Y values greater than 1 to 1
data_lug$Y[data_lug$Y > 1] <- 1
data_lug$Ymin[data_lug$Ymin > 1] <- 1
data_lug$Ymax[data_lug$Ymax > 1] <- 1

# Create the plot
e <- ggplot(data_lug[data_lug$x >= glug/100,], aes(x = x, y = Y)) +
  geom_line(color="#ffd92fff") +
  geom_ribbon(data = data_lug,aes(ymin = Ymin, ymax = Ymax), fill = "#ffd92fff", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c( 0.005, 0.01, 0.1, 1),  
    labels = c( 0.005, 0.01, 0.1, 1)) +  
  labs(x = expression(paste("Percentage of ESBL-", italic("E. coli"), " out of total ", italic("E. coli"), " in the gut")),
       y = expression(paste("Prevalence of ESBL-", italic("E. coli"), " carriage within the community"))) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.020", vjust = 1.3, hjust=-0.0,data = data_lug[data_lug$Y == 1, ], 
                 aes(xintercept = a), linetype = "dashed", color="#ffd92fff", alpha = 0.8)+
  ggtitle("IDA CDA Lugano")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
e

###ARA Sensetal Laupen with CI-------
# Create a data frame with different values of "a": a=gamma_mean, a=gamma_mean - CI, a=gamma_mean + CI
data_sen <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                       a = c(gsen/100))  #a = gamma_mean

data_sen_min <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(gsen/100-CIsen/100)) #a = gamma_mean - CI

data_sen_max <- data.frame(x = seq(0.01, 1, by = 0.0001),  
                           a = c(gsen/100+CIsen/100)) #a = gamma_mean + CI


# Calculate Ys based on the function Y = a/x
data_sen$Y <- data_sen$a / data_sen$x
data_sen$Ymin <- data_sen_min$a / data_sen_min$x
data_sen$Ymax <- data_sen_max$a /data_sen_max$x

# Set Y values greater than 1 to 1
data_sen$Y[data_sen$Y > 1] <- 1
data_sen$Ymin[data_sen$Ymin > 1] <- 1
data_sen$Ymax[data_sen$Ymax > 1] <- 1

# Create the plot
f <- ggplot(data_sen[data_sen$x >= gsen/100,], aes(x = x, y = Y)) +
  geom_line(color="#a6d854ff") +
  geom_ribbon(data = data_sen,aes(ymin = Ymin, ymax = Ymax), fill = "#a6d854ff", alpha = 0.3, colour=NA) +
  scale_x_log10(
    limits = c(0.005, 1),
    breaks = c( 0.005, 0.01, 0.1, 1),  
    labels = c( 0.005, 0.01, 0.1, 1)) +  
  labs(x = expression(paste("Percentage of ESBL-", italic("E. coli"), " out of total ", italic("E. coli"), " in the gut")),
       y = expression(paste("Prevalence of ESBL-", italic("E. coli"), " carriage within the community"))) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)) +
  geom_textvline(label = "0.014", vjust = 1.3, hjust=-0.0,data = data_sen[data_sen$Y == 1, ], 
             aes(xintercept = a), linetype = "dashed", color="#a6d854ff", alpha = 0.7)+
  ggtitle("ARA Sensetal Laupen")+
  theme_minimal()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

f

ggarrange(a,b, c, d, e, f,nrow = 2, ncol = 3, labels = c("A", "B", "C", "D","E", "F"))

##CDF, mean and median of Bangladesh ESBL-E.coli percentage in the gut----------------
#Percentage of ESBL-E. coli out of total E. coli in the gut of children
setwd("/Users/conforsh/switchdrive/Institution/Manuscripts/esblec-monitor-ww-21-22/data")
df_bang <- read.csv("intestinal_carriage_bangladesh.csv", header = TRUE)

# Create empirical CDF data
df_bang$ecdf <- ecdf(df_bang$percentage_esblEc)(df_bang$percentage_esblEc) 

# Fit lognormal distribution to data
fit <- fitdistr(df_bang$percentage_esblEc, "lognormal") 

# Generate fitted CDF values based on the lognormal parameters
df_bang$lognorm_cdf <- plnorm(df_bang$percentage_esblEc, meanlog = fit$estimate[1], sdlog = fit$estimate[2])

# Calculate empirical arithmetic mean and median
mean_val <- mean(df_bang$percentage_esblEc)
median_val <- median(df_bang$percentage_esblEc)

# Generate Plot
cdf=ggplot(df_bang, aes(percentage_esblEc)) + 
  geom_line(aes(y=ecdf, color="Empirical CDF"), size=1) +
  geom_line(aes(y=lognorm_cdf, color="Lognormal CDF"), linetype="solid", linewidth=0.5) +
  scale_x_log10(labels = scales::comma) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  labs(
    x = expression(paste("Percentage of ESBL-", italic("E. coli"), " out of total ", italic("E. coli"), " in the gut of children in Bangladesh (%)")),
    y = "CDF",
    title = paste("Lognormal Parameters: Meanlog =", round(fit$estimate[1], 3), 
                  "Sdlog =", round(fit$estimate[2], 3))) +
  geom_textvline(label = "Median = 2%", vjust = 1.3, hjust=0.01,
                 aes(xintercept = median_val), linetype = "dashed", color="black", alpha =1)+
  geom_textvline(label = "Mean = 19.3%", vjust = 1.3, hjust=0.01,
                 aes(xintercept = mean_val), linetype = "dotted", color="black", alpha =1)+
  scale_color_manual("Legend", 
                     values = c("Empirical CDF" = "blue", "Lognormal CDF" = "red"))

print(cdf)
