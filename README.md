# ESBL-_E. coli_ monitoring in WW
## codes 
### analyses.R 
The analyses.R file that allows to reproduce all the figures, tables and statistical tests within the main manuscript and within the supplementary material. 

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

### DBDA2E-utilities.R 
The DBDA2E-utilities.R file contains the packages required to estimate the sample size using Bayesian inference. The statistics use the package rjags, which work properly on windows systems, but not necessarily on MacOS. 

The file is called by the main code file analyses.R if stored in the correct working directory. 

### process.nf
The process.nf file is a script written in Nextflow. It allows to process the raw sequencing data and to obtain a final .csv file containing the number of reads for each target (ESBL-gene family) in each sample sequenced (barcode). 

The script can be run through the local terminal by typing: **nextflow process.nf**
Nextflow needs to be installed before, as well as the dependencies required. 

The directory in which the script should be run is the one corresponding to "baseDir/" in the first line of the script. Data are stored in a folder named raw, that is nested in a folder named data ($baseDir/data/raw). 

### esbl_dMLA.R
The esbl_dMLA.R file contains the code to compute the overall statistics and the plots related to the dMLA results proposed in the manuscript. 

The code is subdivided in 4 main chapters:
1. Load packages
2. Formatting df
3. Summary statistics
4. Plots 

In order to run the whole code, it is essential to:
- install and load the required packages running the lines 3:9
- format the df running the lines 11:62

The working directory is set at line 13 and needs to be changed depending on the user. It contains the file **reads_counts.csv**.

## data
### ecoli_counts.csv
The file ecoli_counts.csv contains all the data necessary to run the first 6 chapters of the analyses.R code (all except the analyses of the intestinal carriage in Bangladesh). 

Variables in the data file (_na_ indicates that data are not available for that observation and/or that variable):
- **date**: dd_mm_yyyy
- **wwtp**: wastewater treatment plant (Altenrhein, Chur, Geneva, Laupen, Lugano)
- **capacity**: number of inhabitants connected to the wastewater treatment plant
- **temperature**: air temperature data at a height of 2 meters above ground, with a daily average from 6 UTC to 18 UTC (°C)
- **precipitation_96h_sum_mm**: cumulative sum of precipitation over 96 hours leading up to and including the wastewater collection day
- **precipitations_24h_sum_mm**: daily total precipitations from 0 UTC to 0 UTC on the wastewater collection day
- **totalEc_cfu_100ml_a**: number of colony forming units of total _E. coli_ counted in 100mL of wastewater, non-adjusted by flowrate (Replicate 1)
- **totalEc_cfu_100ml_b**: number of colony forming units of total _E. coli_ counted in 100mL of wastewater, non-adjusted by flowrate (Replicate 2)
- **esblEc_cfu_100ml_a**: number of colony forming units of ESBL-_E. coli_ counted in 100mL of wastewater, non-adjusted by flowrate (Replicate 1)
- **esblEc_cfu_100ml_2**: number of colony forming units of ESBL-_E. coli_ counted in 100mL of wastewater, non-adjusted by flowrate (Replicate 2)
- **esblEc_percentage_a**: esblEc_cfu_100ml_a/totalEc_cfu_100ml_a*100 (Replicate 1)
- **esblEc_percentage_b**: esblEc_cfu_100ml_b/totalEc_cfu_100ml_b*100 (Replicate 2)
- **flowrate_m3_day**: flow rate of the wwtp on the sampling day (as cubic meters per day)
- **totalEc_loads_a**: loads of total _E. coli_ expressed as colony forming units per day and per person, the totalEc_cfu_100ml_a variable was normalized by the flowrate and the number of inhabitants in the wwtp (Replicate 1)
- **totalEc_loads_b**: loads of total _E. coli_ expressed as colony forming units per day and per person, the totalEc_cfu_100ml_b variable was normalized by the flowrate and the number of inhabitants in the wwtp (Replicate 2)
- **esblEc_loads_a**: loads of ESBL-_E. coli_ expressed as colony forming units per day and per person, the esblEc_cfu_100ml_a variable was normalized by the flowrate and the number of inhabitants in the wwtp (Replicate 1) 
- **esblEc_loads_b**: loads of ESBL-_E. coli_ expressed as colony forming units per day and per person, the esblEc_cfu_100ml_b variable was normalized by the flowrate and the number of inhabitants in the wwtp (Replicate 2)

### intestinal_carriage_bangladesh.csv
The file intestinal_carriage_bangladesh.csv includes the data of ESBL-_E. coli_ relative to total _E. coli_ in 67 children from Bangladesh, and have been retrieved by a cohort study conducted by Montealegre et _al.,_ 2022 (https://doi.org/10.1289/EHP11359).

Variables in the data file:
- **child_id**: number assigned to the children (1:67)
- **totEcoli_CFUs_gram_feces**: number of colony forming units of total _E. coli_ detected in 1 gram of feces 
- **esblEcoli_CFUs_gram_feces**: number of colony forming units of ESBL-_E. coli_ detected in 1 gram of feces
- **percentage_esblEc**: esblEcoli_CFUs_gram_feces/totEcoli_CFUs_gram_feces*100

### reads_counts.csv
The file reads_counts.csv derives from the file sequencing_targets.csv. Variables (columns) were added and levels were modified in excel manually. 

Variables in the data file:
- **sample**: name of the ESBL-_E. coli_ isolate in the format wwtp_month_#, where # is either 1, 2, or 3
- **month**: month in which the sample was isolated from wastewater
- **dMLA_set**: set in which sample was sequenced, either 1, 2, or 3
- **wwtp**: wastewater treatment plant from where sample was isolated
- **well**: well of the 96-well plate in which the sample undergone dMLA, each well is associated with a unique barcode (listed in the file sequencing_targets.csv)
- **X2**: target ESBL-gene family in the format of C#aT or C#bT, where # is a number and a, b refer to one of the two probe-pairs targeting the same ESBL-gene family
- **n**: number of unique molecular counts for the specific probe-pair in each sample
- **threshold**: threshold calculated for each probe target, as the sum of the mean and three times the standard deviation of unique molecules from false positive signals
- **real_n**: n minus threshold


