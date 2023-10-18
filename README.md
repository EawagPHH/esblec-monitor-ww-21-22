# ESBL-_E. coli_ monitoring in WW
## code 
### analyses.R 
**analyses.R** file that allows to reproduce all the figures, tables and statistical tests within the main manuscript and within the supplementary material. 

The code is subdivided in 7 main chapters:
1. Calculate overall statistics of _E. coli_ ecoli_counts
2. Investigate differences between wastewater treatment plants
    - Plots
    - Statistical tests
3. Investigate differences between months
    - Plots
    - Statistical tests
4. Investigate correlations with environmental variables
5. Estimate the sample size
6. Estimate carriage of ESBL-_E. coli_
7. Analyze distribution of intestinal carriage data in Bangladesh

In order to run the whole code, it is essential to:
- install and load the required packages running the lines 2:22
- format the df running the lines 23:132

The working directory is set three times and needs to be changed depending on the user:
- line 25: setwd where **ecoli_counts.csv** is stored
- line 963: setwd where **DBDA2E-utilities.R** file is stored
- line 2708: setwd where **intestinal_carriage_bangladesh.csv** is stored

The DBDA2E-utilities.R file contains the packages required to estimate the sample size using Bayesian inference. The statistics use the package rjags, which work properly on windows systems, but not necessarily on MacOS. 

## Load required packages
Lines 2:22

**All need to be run to make the code work.**                                                                         
##Formatting df  
Lines 23:132

**All need to be run to make the code work.**                                                                  
##Calculate overall statistics   
Lines 132:201                                              
##Differences between WWTPs                                                                   
###Plots        
Lines 202:349                                                               
####Percentage ESBL-E. coli                                                              
####Loads ESBL-E. coli                                                                
####Loads total E. coli                                                                     
###Statistical tests  
Lines 350:397                                                 
####Percentage ESBL-E. coli                                                                
####Loads ESBL-E. coli                                                                      
####Loads ESBL E. coli                                                                      
##Differences between months                                           
###Plots
Lines                                                                          
####Percentage ESBL-E.coli                                                           
####Loads ESBL-E.coli                                                                  
####Loads total E.coli                                                                 
###Statistical tests                                                               
####Percentage ESBL-E.coli                                                            
####Loads total E.coli                                                                 
####Loads ESBL-E.coli                                                                  
##Environmental variables                                                                       
####ESBL-E. coli percentage                                                            
#####Temperature                                                        
#####Precipitation 24h                                                  
#####Precipitation 96h                                           
#####Plot heatmap                                                
####Loads ESBL-E. coli                                                                   
#####Temperature                                                      
#####Precipitation 24h                                                         
#####Precipitation 96h                                             
#####Plot heatmap                                                
####Loads total E. coli                                                        
#####Temperature                                                                
#####Precipitation 24h                                                  
#####Precipitation 96h                                            
#####Plot heatmap                                                
##Estimate sample size                                                      
###Altenrhein                                                                          
#### Log Normal - Model                                                     
##### Examine results-DIC
##### Plot CCDF         
#### Gamma - Model                                                      
##### Examine results-DIC
##### Plot CCDF        
###Chur                                                                              
####Log Normal - Model                                                      
##### Examine results-DIC
##### Plot CCDF        
#### Gamma - Model                                                          
##### Examine results-DIC
##### Plot CCDF       
###Zürich                                                                        
#### Log Normal - Model                                                    
##### Examine results-DIC
##### Plot CCDF          
#### Gamma - Model                                                          
##### Examine results-DIC 
##### Plot CCDF         
###Geneva                                                                              
#### Log Normal - Model                                                     
##### Examine results-DIC
##### Plot CCDF          
#### Gamma - Model                                                         
##### Examine results-DIC
##### Plot CCDF          
###Lugano                                                                              
#### Log Normal - Model                                                     
##### Examine results-DIC
##### Plot CCDF         
#### Gamma - Model                                                        
##### Examine results-DIC 
##### Plot CCDF          
###Sensetal Laupen                                                                     
#### Log Normal - Model 
##### Examine results-DIC
##### Plot CCDF 
#### Gamma - Model                                                          
##### Examine results-DIC 
##### Plot CCDF        
###95% CI plotting                                                             
##Carriage estimation                                                                           
###Define estimates                                                                   
###Whole Switzerland, with population adjusted ESBL-E. coli percentage in WW               
###ARA Altenrhein with CI                                                                       
###ARA Chur with CI                                                                             
###STEP d'Aïre Genève with CI                                                                  
###ARA Werdhölzli Zürich with CI                                                               
###IDA CDA Lugano with CI                                                                       
###ARA Sensetal Laupen with CI                                                                 
##CDF, mean and median of Bangladesh ESBL-E.coli percentage in the gut
