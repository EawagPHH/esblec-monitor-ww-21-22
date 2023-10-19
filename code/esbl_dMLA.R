# Analyses of ESBL-gene reads----------------------
## Load packages--------------------------
library (ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(ggsci)
library(ggnewscale)
library(dplyr)

## Formatting df---------------------------
#Import dataset
setwd("/Users/conforsh/switchdrive/Institution/Manuscripts/ESBL_ecoli_ww_monitoring")
reads_coli <- read.csv("reads_counts.csv", header = TRUE)

#Rename ESBL-gene families and list them from most to least abundant
reads_coli$X2=as.factor(reads_coli$X2)
levels(reads_coli$X2)=list("PDC.1"="C4aT","PDC.2"= "C4bT","VIM.1"="C10aT", "VIM.2"= "C10bT",
                           "GES.1"="C13aT", "GES.2"="C13bT", "ADC.1"="C6aT", "ADC.2"= "C6bT", 
                           "OXA-51.1"="C1aT","OXA-51.2"="C1bT","OXA-213.1"="C12aT","OXA-213.2"="C12bT",
                           "OXA-48.1"="C20aT", "OXA-48.2"="C20bT","CARB.1"="C15aT", "CARB.2"="C15bT", 
                           "CMY.1"="C3aT","CMY.2"="C3bT","CTX-M2.1"="C18aT", "CTX-M2.2"="C18bT",
                           "IMP.1"="C9aT", "IMP.2"="C9bT", "SHV.1"="C2aT", "SHV.2"="C2bT",
                           "ACT/MIR.1"="C7aT","ACT/MIR.2"="C7bT",
                           "CTX-M8 and 25.1"="C27aT","CTX-M8 and 25.2"="C27bT",
                           "CMY,FOX,MOX.1"="C14aT","CMY,FOX,MOX.2"="C14bT",
                           "TEM.1"="C0aT","TEM.2"= "C0bT", "CTX-M9.1"= "C8aT", "CTX-M9.2"="C8bT",
                            "CTX-M1.1"="C5aT", "CTX-M1.2"="C5bT")

#Create a new variable in which the two probe-pairs of each gene family have the same name
reads_coli$X3 <- factor(reads_coli$X2, levels = levels(reads_coli$X2))
levels(reads_coli$X3) = list("PDC"="PDC.1","PDC"="PDC.2","VIM"="VIM.1","VIM"= "VIM.2","GES"="GES.1","GES"="GES.2",
                             "ADC"="ADC.1", "ADC"="ADC.2","OXA-51"="OXA-51.1", "OXA-51"="OXA-51.2",
                             "OXA-213"= "OXA-213.1","OXA-213"="OXA-213.2","OXA-48"="OXA-48.1","OXA-48"="OXA-48.2",
                             "CARB"="CARB.1", "CARB"="CARB.2","CMY"="CMY.1","CMY"="CMY.2",
                             "CTX-M2"="CTX-M2.1","CTX-M2"="CTX-M2.2","IMP"="IMP.1", "IMP"="IMP.2",
                             "SHV"="SHV.1", "SHV"="SHV.2","ACT/MIR"="ACT/MIR.1","ACT/MIR"="ACT/MIR.2",
                             "CTX-M8 and 25"="CTX-M8 and 25.1","CTX-M8 and 25"="CTX-M8 and 25.2",
                             "CMY,FOX,MOX"="CMY,FOX,MOX.1","CMY,FOX,MOX"="CMY,FOX,MOX.2","TEM" ="TEM.1","TEM" ="TEM.2", 
                             "CTX-M9"="CTX-M9.1", "CTX-M9"="CTX-M9.2","CTX-M1"="CTX-M1.1", "CTX-M1"="CTX-M1.2")
#Rename the wwtp
reads_coli$wwtp[which(reads_coli$wwtp=="Altenrhein")] <- "ARA Altenrhein"
reads_coli$wwtp[which(reads_coli$wwtp=="Lugano")] <- "IDA CDA Lugano"
reads_coli$wwtp[which(reads_coli$wwtp=="Geneva")] <- "STEP d'Aïre Genève"
reads_coli$wwtp[which(reads_coli$wwtp=="Laupen")] <- "ARA Sensetal Laupen"
reads_coli$wwtp[which(reads_coli$wwtp=="Zurich")] <- "ARA Werdhölzli Zürich"
reads_coli$wwtp[which(reads_coli$wwtp=="Chur")] <- "ARA Chur"

#Rename and reorder the months in the dataframe
reads_coli$month=as.factor(reads_coli$month)
levels(reads_coli$month)=list("November 2021" = "Nov", "December 2021"="Dec", "January 2022" = "Jan",
                              "February 2022" = "Feb", "March 2022"="Mar", "April 2022"="Apr","May 2022"="May", 
                              "June 2022"="Jun", "July 2022"="Jul", "2nd August 2022"="Aug2nd", "30th August 2022"="Aug30th", 
                              "September 2022"="Sep", "October 2022"="Oct")

#Assign NA values when the number of reads is 0.
reads_coli$real_n[reads_coli$real_n == 0] <- NA

# Reorder levels of sample based on month
reads_coli$sample <- factor(reads_coli$sample, 
                                     levels = unique(reads_coli[order(reads_coli$month), ]$sample))

## Summary statistics and plots ------------------------
### Positive to at least one probe-pair-------------
#Count how many isolates are positive to at least one probe-pair for at least on gene family
print(reads_coli %>% 
        summarise(n_distinct(sample)))

#Compute the number of isolates positive to each probe-pair in each wwtp
result_1 = reads_coli %>%
  filter(real_n > 0) %>%
  group_by(wwtp, X2) %>%
  summarize(unique_samples = n_distinct(sample))

#Compute the number of isolates positive to each probe-pair in each month
result_2 = reads_coli %>%
  filter(real_n > 0) %>%
  group_by(month, X2) %>%
  summarize(unique_samples = n_distinct(sample))

# Remove rows where the column real_n has missing values (NA) or where real_n is equal to 0
reads_coli_filtered=reads_coli %>%
  filter(!is.na(real_n) & real_n != 0)

# Plot with presence absence of genes positive for at least one probe-pair
ggplot(data=reads_coli_filtered, aes(sample, X2)) + 
  geom_tile(colour = "black") +
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 0.5, vjust = 0.5), axis.text.y = element_text(size = 4)) +
  ylab("Genes clusters") +
  xlab(expression(paste(" ", italic("E. coli"), "isolates"))) +
  facet_wrap(~wwtp, ncol = 3, scales = "free_x")

### Positive only to both probe-pairs of same gene family-----------
#Subset for isolates positive to both probe pairs 
result <- reads_coli %>%  
  filter(real_n > 0)%>%
  group_by(X3, sample) %>%
  filter(n() >= 2) %>%
  ungroup() 

#Count how many isolates are positive to both probe pairs for at least one gene family
print(result %>% 
  summarise(n_distinct(sample)))

# Compute the number of isolates positive to both probe pairs of each gene family in each wwtp
result_wwtp = result %>% 
  filter(real_n > 0) %>%
  group_by(wwtp, X3)  %>%
  summarize(unique_sample = n_distinct(sample))%>%
  mutate(percentage = unique_sample / 39 * 100)

# Compute the number of isolates positive to both probe pairs of each gene family combining wwtp
result_X3 = result %>%   
  filter(real_n > 0) %>%  
  group_by(X3)  %>%
  summarize(unique_sample = n_distinct(sample))%>%
  mutate(percentage = unique_sample / 234 * 100)

# Compute the number of isolates positive to both probe pairs of each gene family in each month combining wwtp
result_month = result%>%
  filter(real_n > 0) %>%
  group_by(month, X3)  %>%
  summarize(unique_sample = n_distinct(sample))%>%
  mutate(percentage = unique_sample / 18 * 100)

# Remove rows where the column real_n has missing values (NA) or where real_n is equal to 0
reads_coli_filtered=result %>%
  filter(!is.na(real_n) & real_n != 0)

#Plot with presence absence of genes positive for both probe-pairs
ggplot(data=reads_coli_filtered, aes(sample, X2)) + 
  geom_tile(colour = "black") +
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 0.5, vjust = 0.5), axis.text.y = element_text(size = 4)) +
  ylab("Genes clusters") +
  xlab(expression(paste(" ", italic("E. coli"), "isolates"))) +
  facet_wrap(~wwtp, ncol = 3, scales = "free_x")

#Format result_wwtp df to have numeric number of samples and wwtp as factor
result_wwtp$unique_sample=as.numeric(result_wwtp$unique_sample)
result_wwtp$wwtp=as.factor(result_wwtp$wwtp)

#Plot the number of isolates positive to both probe-pair of each gene-family and facet by wwtp
ggplot(result_wwtp, aes(x = X3, y = unique_sample)) +
  geom_bar(stat = "identity", fill = "black") +
  xlab("ESBL-gene family") +
  ylab("n° isolates") +
  facet_wrap(~wwtp, ncol = 3)+
  coord_flip()

#Plot the number of isolates positive to both probe-pair of each gene-family and facet by wwtp, colourcoded
ggplot(complete(result_wwtp,X3), aes(x = X3, y = unique_sample, fill = X3)) + 
  geom_bar(stat = "identity", position = "dodge") +
  xlab("ESBL-gene family") +
  ylab(expression(paste("Number of ",italic("E. coli "), "isolates"))) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                                        "gray", "red3", "purple", "lightgreen", "gray25", "coral4",
                                        "pink", "lightblue", "darkgreen", "blue4", "black", "orange4")) +
  theme(axis.text.x =  element_text(size = 7, angle=90, hjust = 1),
        legend.position = "none") +
  facet_wrap(~wwtp) 

#Plotting MONTHS on the x axis and Facet-wrap of the same (Number of E. coli isolates)
ggplot(result_month, aes(x = X3, y = percentage, fill = X3)) + 
  geom_bar(stat = "identity", position = "dodge") +
  xlab("ESBL-gene family") +
  ylab(expression(paste("Percentage of ESBL-",italic("E. coli "), "isolates (%)"))) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                                        "gray", "red3", "purple", "lightgreen", "gray25", "coral4",
                                        "pink", "lightblue", "darkgreen", "blue4", "black", "orange4")) +
  theme(axis.text.x =  element_text(size = 10, angle=90, hjust=1, colour = "black"),
        axis.text.y =  element_text(size = 10, angle=0, hjust=1, colour = "black"),
        strip.text = element_text(size = 14),
        legend.position = "none") +
  facet_wrap(~month) 

#Calculate percentage of how many samples positive to 2 PROBE PAIRS in each MONTH for each X3
result_MONTH_per = result %>%
  filter(real_n > 0) %>%
  group_by(Month, X3) %>%
  summarize(unique_samples = n_distinct(sample))%>%
  mutate(percentage = unique_samples / 18 * 100)

#Plotting MONTHS on the x axis and Facet-wrap of the same (% of E. coli isolates)
ggplot(result_MONTH_per, aes(x = Month, y = percentage, fill = X3)) + 
  geom_bar(stat = "identity", position = "dodge") +
  xlab("") +
  ylab(expression(paste("Percentage of ",italic("E. coli "), "isolates (%)"))) +
  #ggtitle("Counts of each target gene by number of isolates and Month") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                                        "gray", "red3", "purple", "lightgreen", "gray25", "coral4",
                                        "pink", "lightblue", "darkgreen", "blue4", "black", "orange4")) +
                                          guides(fill = guide_legend(title = "ESBL-genes clusters", nrow = 1)) +
  theme(axis.text.x =  element_blank(), axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.key.size = unit(0.5, "cm"),
        legend.text.align = 0) +
  guides(fill = guide_legend(title = "ESBL-genes clusters", ncol=2, byrow = TRUE))+
  facet_wrap(~Month, scales = "free_x") 

#same as above but no legend
ggplot(result_MONTH_per, aes(x = Month, y = percentage, fill = X3)) + 
  geom_bar(stat = "identity", position = "dodge") +
  xlab("") +
  ylab(expression(paste("Percentage of ",italic("E. coli "), "isolates (%)"))) +
  #ggtitle("Counts of each target gene by number of isolates and Month") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                                        "gray", "red3", "purple", "lightgreen", "gray25", "coral4",
                                        "pink", "lightblue", "darkgreen", "blue4", "black", "orange4")) +
                                          theme(axis.text.x =  element_text(size = 10, angle=0)) +
  facet_wrap(~Month, scales = "free_x")+
  theme(legend.position = "none", axis.text.x =  element_blank(), axis.ticks.x = element_blank())

#----------------a-----------------#
#Calculate how many isolate are positive to each target gene, showing comparisons between Seasons
df_subset <- reads_coli %>% select(Season, sample, real_n, X2)

df_subset <- df_subset %>% 
  filter(real_n > 0) %>% 
  group_by(X2, Season) %>% 
  summarize(count = n_distinct(sample)) %>% 
  ungroup()

ggplot(df_subset, aes(x = X2, y = count, fill = Season)) + 
  geom_bar(stat = "identity", position = "stack") +
  xlab("Target Gene") +
  ylab(expression(paste("Number of ",italic("E. coli "), "isolates"))) +
  ggtitle("Counts of each target gene by number of isolates and season") +
  scale_fill_manual(values = c("brown", "gray", "pink", "red"), labels = c("Fall", "Winter", "Spring", "Summer")) +
  guides(fill = guide_legend(title = "Month", nrow = 1)) +
  theme(axis.text.x = element_text(angle=90))+
  theme(legend.position = "bottom")

##Subset based on wwtp
alt_sub <- subset(reads_coli, wwtp == "Altenrhein")
alt=ggplot(complete(alt_sub,X2, sample), aes(X2,sample, fill= real_n))+ 
  geom_tile(colour="white")+
  scale_fill_gradient(low = "lightblue", high = "blue", na.value = "gray95",  name ='Molecular counts')+
  theme(axis.text.x = element_text(size=5, angle=90))

chur_sub <- subset(reads_coli, wwtp == "Chur")
chur=ggplot(complete(chur_sub,X2, sample), aes(X2,sample, fill= real_n))+ 
  geom_tile(colour="white")+
  scale_fill_gradient(low = "lightblue", high = "blue", na.value = "gray95",  name ='Molecular counts')+
  theme(axis.text.x = element_text(size=5, angle=90))

lug_sub <- subset(reads_coli, wwtp == "Lugano")
lug=ggplot(complete(lug_sub,X2, sample), aes(X2,sample, fill= real_n))+ 
  geom_tile(colour="white")+
  scale_fill_gradient(low = "lightblue", high = "blue", na.value = "gray95",  name ='Molecular counts')+
  theme(axis.text.x = element_text(size=5, angle=90))

lau_sub <- subset(reads_coli, wwtp == "Laupen")
lau=ggplot(complete(lau_sub,X2, sample), aes(X2,sample, fill= real_n))+ 
  geom_tile(colour="white")+
  scale_fill_gradient(low = "lightblue", high = "blue", na.value = "gray95",  name ='Molecular counts')+
  theme(axis.text.x = element_text(size=5, angle=90))

gen_sub <- subset(reads_coli, wwtp == "Geneva")
gen=ggplot(complete(gen_sub,X2, sample), aes(X2,sample, fill= real_n))+ 
  geom_tile(colour="white")+
  scale_fill_gradient(low = "lightblue", high = "blue", na.value = "gray95",  name ='Molecular counts')+
  theme(axis.text.x = element_text(size=5, angle=90))

zur_sub <- subset(reads_coli, wwtp == "Zurich")
zur=ggplot(complete(zur_sub,X2, sample), aes(X2,sample, fill= real_n))+ 
  geom_tile(colour="white")+
  scale_fill_gradient(low = "lightblue", high = "blue", na.value = "gray95",  name ='Molecular counts')+
  theme(axis.text.x = element_text(size=5, angle=90))

grid.arrange(alt, chur, lug, lau, gen, zur, ncol=3, nrow = 2)

##Depending on place
coco=reads_coli %>%
  group_by(X2, wwtp, Isolate_n) %>%
  summarise(Isolate_n = sum(Isolate_n)) %>%
  complete( fill = list(Isolate_n = 0))%>%
  summarise(Isolate_n = sum(Isolate_n))

coco =as.data.frame(coco)
coco[coco == 0] <- NA
ggplot(complete(coco, X2, wwtp), aes(X2,wwtp, fill= Isolate_n))+ 
  geom_tile(colour="white")+
  scale_fill_gradient(low = "lightblue", high = "blue", na.value = "gray95",  
                      name ='Isolates positive to the gene', limits=c(1,13), n.breaks=12)+
  theme(axis.text.x = element_text(size=5, angle=90))

ll=ggplot(coco, aes(x = X2, y = Isolate_n)) +
  geom_bar(stat = "identity", width=.5, position = "dodge", colour="lightblue", fill="lightblue")+ 
  theme(axis.text.x = element_text(size=5, angle=90))

ll+facet_wrap(.~wwtp, ncol=3)

##Depending on month
gaga=reads_coli %>%
  group_by(X2, Month, Isolate_n) %>%
  summarise(Isolate_n = sum(Isolate_n)) %>%
  complete( fill = list(Isolate_n = 0))%>%
  summarise(Isolate_n = sum(Isolate_n))

gaga =as.data.frame(gaga)
gaga[gaga == 0] <- NA
ggplot(complete(gaga, X2, Month), aes(X2,Month, fill= Isolate_n))+ 
  geom_tile(colour="white")+
  scale_fill_gradient(low = "lightblue", high = "blue", na.value = "gray95",  
                      name ='Isolates positive to the gene', limits=c(1,15), n.breaks=15)+
  theme(axis.text.x = element_text(size=5, angle=90))

ee=ggplot(gaga, aes(x = X2, y = Isolate_n, fill=Month)) +
  geom_bar(stat = "identity", width=.5, position = "dodge")+ 
  theme(axis.text.x = element_text(size=5, angle=90))

ee+facet_wrap(.~Month, ncol=5)

ff=ggplot(gaga, aes(x = X2, y = Isolate_n, fill=X2)) +
  geom_bar(stat = "identity", width=.5, position = "dodge")+ 
  theme(axis.text.x = element_text(size=5, angle=90))

ff+facet_wrap(.~Month, ncol=5)


###Assign colours to values in ggplot heatmap: https://statisticsglobe.com/change-colors-of-ranges-in-ggplot2-heatmap-r
reads$real_n_groups=as.factor(reads$real_n_groups)
levels(reads$real_n_groups) <- list(">10'000" = "(10000,100000]","1'000-10'000" = "(1000,10000]", 
                                    "100-1'000" = "(100,1000]","10-100" = "(10,100]",
                                    "1-10" = "(1,10]", "0-1" = "(0,1]")

ggplot(data = na.omit(reads), aes(X2,wwtp, fill= real_n_groups))+ 
  geom_tile(colour="black")+
  theme_grey()+
  scale_fill_manual(breaks = levels(reads$real_n_groups), name ='n° of reads',
                    values = c("red4","red3","red","lightpink","grey70","grey90"))+
  theme(axis.text.x = element_text(size=9, angle=90, colour="black"), 
        axis.text.y = element_text(size = 5, colour="black"))+
  xlab("Target genes")

