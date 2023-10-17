#------Sample Size estimation--------
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

#Create an average value of the E.coli data replicates and use the non-missing value if one of the two replicates is missing
df$average_ESBL_Ec <- as.numeric(ifelse(is.na(df$esblEc_percentage_a) | is.na(df$esblEc_percentage_b), 
                                        ifelse(is.na(df$esblEc_percentage_a), df$esblEc_percentage_b, df$esblEc_percentage_a), 
                                        rowMeans(df[,c("esblEc_percentage_a", "esblEc_percentage_b")], na.rm = TRUE)))

#11) Filter df to keep only non NAs values
df = filter(df, !is.na(average_ESBL_Ec)) #ESBL-Ec percentage (%)

#12) Subset df based on WWTP
df_alt = df %>% filter(wwtp == "ARA Altenrhein")
df_chu = df %>% filter(wwtp == "ARA Chur")
df_gen = df %>% filter(wwtp == "STEP d'Aïre Genève")
df_zur = df %>% filter(wwtp == "ARA Werdhölzli Zürich")
df_lug = df %>% filter(wwtp == "IDA CDA Lugano")
df_sen = df %>% filter(wwtp == "ARA Sensetal Laupen")

#13) Bayesian code to fit log normal distribution to my data
graphics.off()
setwd("C:/Users/Sheena/Desktop/Sample Size")
source("DBDA2E-utilities.R")
fileNameRoot="A" 
library(rjags)
library(runjags)

##Altenrhein----------------
df_alt = df %>% filter(wwtp == "ARA Altenrhein")
C = as.numeric(df_alt[, "average_ESBL_Ec"])
N = length(C)

### Log Normal - Model ---------------------------
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

#### Examine results-DIC------------------------------------------------------------------------------
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
#### Plot CCDF ------------------------------------------------------------------------------
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

### Gamma - Model ---------------------------
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

#### Examine results-DIC-----------------------------------------------------------------------------
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
#### Plot CCDF ------------------------------------------------------------------------------
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

#Chur----------------
df_chu = df %>% filter(wwtp == "ARA Chur")
C = as.numeric(df_chu[, "average_ESBL_Ec"])
N = length(C)

### Log Normal - Model ---------------------------
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

#### Examine results-DIC------------------------------------------------------------------------------
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
#### Plot CCDF ------------------------------------------------------------------------------
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

### Gamma - Model ---------------------------
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

#### Examine results-DIC-----------------------------------------------------------------------------
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
#### Plot CCDF ------------------------------------------------------------------------------
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

#Zürich----------------
df_zur = df %>% filter(wwtp == "ARA Werdhölzli Zürich")
C = as.numeric(df_zur[, "average_ESBL_Ec"])
N = length(C)

### Log Normal - Model ---------------------------
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

#### Examine results-DIC------------------------------------------------------------------------------
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
#### Plot CCDF ------------------------------------------------------------------------------
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

### Gamma - Model ---------------------------
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

#### Examine results-DIC-----------------------------------------------------------------------------
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
#### Plot CCDF ------------------------------------------------------------------------------
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

#Geneva----------------
df_gen = df %>% filter(wwtp == "STEP d'Aïre Genève")
C = as.numeric(df_gen[, "average_ESBL_Ec"])
N = length(C)

### Log Normal - Model ---------------------------
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

#### Examine results-DIC------------------------------------------------------------------------------
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
#### Plot CCDF ------------------------------------------------------------------------------
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

### Gamma - Model ---------------------------
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

#### Examine results-DIC-----------------------------------------------------------------------------
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
#### Plot CCDF ------------------------------------------------------------------------------
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

#Lugano----------------
df_lug = df %>% filter(wwtp == "IDA CDA Lugano")
C = as.numeric(df_lug[, "average_ESBL_Ec"])
N = length(C)

### Log Normal - Model ---------------------------
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

#### Examine results-DIC------------------------------------------------------------------------------
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
#### Plot CCDF ------------------------------------------------------------------------------
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

### Gamma - Model ---------------------------
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

#### Examine results-DIC-----------------------------------------------------------------------------
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
#### Plot CCDF ------------------------------------------------------------------------------
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

#Sensetal Laupen----------------
df_sen = df %>% filter(wwtp == "ARA Sensetal Laupen")
C = as.numeric(df_sen[, "average_ESBL_Ec"])
N = length(C)

### Log Normal - Model ---------------------------
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

#### Examine results-DIC------------------------------------------------------------------------------
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
#### Plot CCDF ------------------------------------------------------------------------------
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

### Gamma - Model ---------------------------
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

#### Examine results-DIC-----------------------------------------------------------------------------
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
#### Plot CCDF ------------------------------------------------------------------------------
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

#S.E. plotting-------------------------
sample_sizes = seq(1, 156, by = 1) #only for graphs up to 3x/week
true_mean = 0

df_alt = df %>% filter(wwtp == "ARA Altenrhein")
alpha_alt = mean(df_alt$average_ESBL_Ec)^2/sd(df_alt$average_ESBL_Ec)^2
beta_alt = mean(df_alt$average_ESBL_Ec)/sd(df_alt$average_ESBL_Ec)^2
gmean_alt=alpha_alt/beta_alt
CI_alt_data <- (sqrt(alpha_alt / beta_alt^2 /50))*1.96
se_alt <- sqrt(alpha_alt / beta_alt^2 / sample_sizes)
CI_alt = se_alt*1.96

data_alt <- data.frame(sample_sizes=sample_sizes, wwtp="ARA Altenrhein",
                       CI=CI_alt)

df_chu = df %>% filter(wwtp == "ARA Chur")
alpha_chu = mean(df_chu$average_ESBL_Ec)^2/sd(df_chu$average_ESBL_Ec)^2
beta_chu = mean(df_chu$average_ESBL_Ec)/sd(df_chu$average_ESBL_Ec)^2
gmean_chu=alpha_chu/beta_chu
CI_chu_data <- (sqrt(alpha_chu / beta_chu^2 /51))*1.96
se_chu <- sqrt(alpha_chu / beta_chu^2 / sample_sizes)
CI_chu = se_chu*1.96

data_chu <- data.frame(sample_sizes=sample_sizes, wwtp="ARA Chur",
                       CI=CI_chu)

df_zur = df %>% filter(wwtp == "ARA Werdhölzli Zürich")
alpha_zur = mean(df_zur$average_ESBL_Ec)^2/sd(df_zur$average_ESBL_Ec)^2
beta_zur = mean(df_zur$average_ESBL_Ec)/sd(df_zur$average_ESBL_Ec)^2
gmean_zur=alpha_zur/beta_zur
CI_zur_data <- (sqrt(alpha_zur / beta_zur^2 / 51))*1.96
se_zur <- sqrt(alpha_zur / beta_zur^2 / sample_sizes)
CI_zur = se_zur*1.96

data_zur <- data.frame(sample_sizes=sample_sizes, wwtp="ARA Werdhölzli Zürich",
                       CI=CI_zur)

df_gen = df %>% filter(wwtp == "STEP d'Aïre Genève")
alpha_gen = mean(df_gen$average_ESBL_Ec)^2/sd(df_gen$average_ESBL_Ec)^2
beta_gen = mean(df_gen$average_ESBL_Ec)/sd(df_gen$average_ESBL_Ec)^2
gmean_gen=alpha_gen/beta_gen
CI_gen_data <- (sqrt(alpha_gen / beta_gen^2 /51))*1.96
se_gen <- sqrt(alpha_gen / beta_gen^2 / sample_sizes)
CI_gen = se_gen*1.96

data_gen <- data.frame(sample_sizes=sample_sizes, wwtp="STEP d'Aïre Genève",
                       CI=CI_gen)

df_sen = df %>% filter(wwtp == "ARA Sensetal Laupen")
alpha_sen = mean(df_sen$average_ESBL_Ec)^2/sd(df_sen$average_ESBL_Ec)^2
beta_sen = mean(df_sen$average_ESBL_Ec)/sd(df_sen$average_ESBL_Ec)^2
gmean_sen=alpha_sen/beta_sen
CI_sen_data <- (sqrt(alpha_sen / beta_sen^2 / 51))*1.96
se_sen <- sqrt(alpha_sen / beta_sen^2 / sample_sizes)
CI_sen = se_sen*1.96

data_sen <- data.frame(sample_sizes=sample_sizes, wwtp="ARA Sensetal Laupen",
                       CI=CI_sen)

df_lug = df %>% filter(wwtp == "IDA CDA Lugano")
alpha_lug = mean(df_lug$average_ESBL_Ec)^2/sd(df_lug$average_ESBL_Ec)^2
beta_lug = mean(df_lug$average_ESBL_Ec)/sd(df_lug$average_ESBL_Ec)^2
gmean_lug=alpha_lug/beta_lug
CI_lug_data <- (sqrt(alpha_lug / beta_lug^2 / 50))*1.96
se_lug <- sqrt(alpha_lug / beta_lug^2 / sample_sizes)
CI_lug = se_lug*1.96

data_lug <- data.frame(sample_sizes=sample_sizes, wwtp="IDA CDA Lugano",
                       CI=CI_lug)

#Bind different df
data <- rbind(data_alt, data_chu, data_zur, data_gen, data_sen, data_lug)
# Plot
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
