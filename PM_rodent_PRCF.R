# Pine marten cyclicity diagnosis_PRCF 
# Author: Jennifer Angoh
# Date: 07.02.23

# Load PM transect data 
# Set working directory
rm(list=ls()) 
setwd("C:/Marten_Vole_Codes")
PMtracks_cov_vole <- readRDS("PMtracks_cov_vole.rds") 
length(unique(PMtracks_cov_vole$Linenr))

# 1. Group transects into 60 transect groups ###################################
# Transect groups constitute new individual time-series for investigations of temporal variation in the PM population
require(sp)
require(rgdal)
library(spatstat)
library(plyr)
citation("raster")
citation("spatstat")
citation("base")

# Calculate mean of start and stop xy
coordinates <- ddply(PMtracks_cov_vole,~ Linenr, summarise, x=mean(X), y=mean(Y))
x <- coordinates$x
y <- coordinates$y
xy <- as.matrix(coordinates[,c(2,3)])
X <- as.ppp(xy, c(537575.5,702857, 6641574,6947303)) #smallest and largest values of x and y coordinates in c()
npoints(X)

# Group 600 transects into 100 groups
Z <- cut(X, seq(1,600,1), breaks=60)

# Add marks column in coordinates 
coordinates$marks <- Z$marks
length(unique(Z$marks))

# make group IDs
coordinates$groupID <- as.numeric(factor(coordinates$marks))

# combine transect ID with group ID
groups <- coordinates[,c(1,5)]

# Join groups to PMtracks_cov_vole by Linenr
library(tidyverse)
PMtracks_cov_group <- left_join(PMtracks_cov_vole, groups, by = "Linenr")


# 2. PM tracks index per year per group ########################################
# With the pooling procedure, 
# For every year, every group: Add snow-track data (Marten) and offsets (SnowDays + length)

# Subset columns 
PMtracks_PRCF <- PMtracks_cov_group[, c("groupID","Year", "Marten", "SnowDays","length")] 

# Remove rows that have no Marten data (NA)
PMtracks_PRCF <- PMtracks_PRCF[!is.na(PMtracks_PRCF$Marten),]

# Remove rows that have no SnowDays data (NA)
PMtracks_PRCF <- PMtracks_PRCF[!is.na(PMtracks_PRCF$SnowDays),]

# Create a total Marten, total SnowDays and total length per year for each groupID 
PMtracks_PRCF <- PMtracks_PRCF %>%
  group_by(groupID, Year) %>%
  mutate(total_Marten = ifelse(all(is.na(Marten)),NA_real_,sum(Marten, na.rm=TRUE)))%>%
  mutate(total_SnowDays = sum(SnowDays, na.rm=TRUE)) %>%
  mutate(total_length = sum(length, na.rm=TRUE))

# Subset dataframe with distinct groupID, total_length, total_Marten, total_SnowDays
PMtracks_PRCF <- PMtracks_PRCF %>% 
  distinct(groupID, total_Marten, total_length, total_SnowDays, .keep_all = F)

# 3. Add smallest obs. to transect groups with zero counts  ####################
#	For remaining zero counts (not to be confused with not surveyed) still remaining after grouping, 
# add the smallest observable entity possible (Turchin 2003), which in our case was 1 crossing PM track

## 118 group years had zero marten tracks, replace by 1 track
library(stringr)
PMtracks_PRCF$total_Marten <- str_replace(PMtracks_PRCF$total_Marten, "0", "1")
# Convert character to numeric in total_Marten column 
PMtracks_PRCF$total_Marten <- as.numeric(PMtracks_PRCF$total_Marten)

# Calculate PM tracks index 
colnames(PMtracks_PRCF)
PMtracks_PRCF$PMIndex <- PMtracks_PRCF$total_Marten/PMtracks_PRCF$total_length/PMtracks_PRCF$total_SnowDays

PMtracks_PRCF <- PMtracks_PRCF[, c("groupID","Year", "PMIndex")] 

#Save input for PRCF 
saveRDS(PMtracks_PRCF, file = "PMtracks_PRCF.rds")

# Find out which group has incomplete timeseries
library(officer)
library(dplyr)
library(tidyr)
# Order on group and year 
# Order by year first
PMtracks_PRCF<-PMtracks_PRCF[order(PMtracks_PRCF$Year),]
# Order by groupID next
PMtracks_PRCF<-PMtracks_PRCF[order(PMtracks_PRCF$groupID),]
PRCF_timeseries <- PMtracks_PRCF %>% 
  pivot_wider(
    names_from = Year,
    values_from = PMIndex
  ) 

#saveRDS(PRCF_timeseries, file = "PRCF_timeseries.rds")

# 4. Create function for PRCF ##################################################
# Partial Rate Correlation Function 
PRCF<- function(x1, log=TRUE){
  r1<-pacf(x1, lag.max = 12, plot = FALSE, na.action = na.omit) 
  PACF<- r1$acf
  v<-c(NULL)
  t<-length(x1)-1
  for(i in 1:t){
    if(log==TRUE) {
      vi<-log(x1[i+1]/x1[i])
      v<-c(v,vi)
    } else {
      vi <- x1[i+1] - x1[i]
      v <- c(v,vi)
    }
  }
  v<-v
  N<-c(NULL)
  for(i in 1:t){
    N<-c(N, x1[i])
  }
  N<-N
  df1 <- data.frame(v, N)
  df1 <- na.omit(df1)
  l1<-cor(df1$v, df1$N)
  PACF[1]<-l1
  prcf<-PACF
  v<-NULL
  for(i in 1:12){
    v[i]<-prcf[i]
  }
  barlett <- 2/sqrt(length(na.omit(x1)))
  names<-c("1","2","3","4","5","6","7","8", "9", "10", "11", "12")
  {barplot(v, col="darkgrey", width = 0.9, axis.lty=1, space = 0.1, names.arg=names, ylim=c(-1,1), ylab="PRCF", xlab="Lags", font.lab=2)
    abline(h=0)
    abline(h=barlett,lty=2, lwd=2)
    abline(h=-barlett, lty=2, lwd=2)
  }
  V<-v
}

# 5. Apply PRCF function to dataset ############################################
# Some of the 60 groups do not have data for all years 2003-2014

### Groups with complete time series 
# Check how many groups have complete time series (12 years) 
library(data.table)
PMtracks_PRCF <- data.table(PMtracks_PRCF)
timeseries_N <- PMtracks_PRCF[, .N, by = list(groupID)]
which(grepl(12, timeseries_N$N))
#Groups with complete time series: 2  5 15 16 20 25 26 27 44 45 54 58 59
# Use these groups to make sure that first lag is the important lag for residuals
complete_ts <- c(2,5,15,16,20,25,26,27,44,45,54,58,59)
# Only keep groups with complete time series 
PMtracks_comp_ts <- PMtracks_PRCF[PMtracks_PRCF$groupID %in% complete_ts,]

### Calculate residuals PRCF[1] of groups with complete time series (n= 13) 
# Make sure order on group and year is right before running the loop 
# Loop going through all groupID 
resi_comp_ts_lag_1 <- data.frame() #Empty dataframe for residuals
for (group in unique(PMtracks_comp_ts$groupID)) { # for each unique group 
  tmp <- PMtracks_comp_ts[which(PMtracks_comp_ts$groupID == group),] # temporary groupID for each loop
  tmptmp <-
    data.frame(lag = 1:(length(PRCF(tmp$PMIndex)) - 1), resid = NA) #fill lag column with num of lags and empty resid columnn 
  for (lag in 1:(length(PRCF(tmp$PMIndex)) - 1)) { # loop through each lag 
    tmptmp[lag, 'resid'] <- PRCF(tmp$PMIndex)[lag] # fill lag with respective resid from PRCF
  }
  tmptmp$abs_resid <- abs(tmptmp$resid) #convert to abs val. 
  output_lag_1 <-
    data.frame(resid = tmptmp[which(tmptmp$abs_resid == (abs(tmptmp$resid)[1])), 'resid'], #row for resid
               group = group, #row for groupID
               lag = tmptmp[which(tmptmp$abs_resid == (abs(tmptmp$resid)[1])), 'lag']) #row for lag of highest absolute resid 
  # Add residual, group and lag to dataframe 
  resi_comp_ts_lag_1 <- rbind(resi_comp_ts_lag_1, output_lag_1) # bind row with each of the group 
}
# Save residuals lag 1
saveRDS(resi_comp_ts_lag_1, file = "resi_comp_ts_lag_1.rds")

### Groups with incomplete time series 
# Extract longest sequence of monitoring in a row into a new time series 
# Add column with sequence groups with consecutive years
# Join cumsum column to PMtracks
PMtracks_PRCF$seq <- cumsum(c(TRUE, diff(PMtracks_PRCF$Year) != 1))

# Get sequence groups that have at least 5 years in a row time series
# First get all row number belonging to each sequence group
count_seq <- split(seq_along(PMtracks_PRCF$Year), cumsum(c(TRUE, diff(PMtracks_PRCF$Year) != 1)))
# Number of rows in each sequence group 
seq_length <- as.data.frame(lengths(count_seq))
# Add column with row number 
seq_length$seq <- seq.int(nrow(seq_length))

# Join seq_length to PMtracks_PRCF
PMtracks_incomp_ts <- left_join(PMtracks_PRCF, seq_length, by = "seq")

# Only keep sequence group that have between 4-11 years in a row
library(dplyr)
PMtracks_incomp_ts <- rename(PMtracks_incomp_ts, length = "lengths(count_seq)")
PMtracks_incomp_ts  <- subset(PMtracks_incomp_ts , length > 3 )
PMtracks_incomp_ts  <- subset(PMtracks_incomp_ts , length < 12 )

length(unique(PMtracks_incomp_ts$groupID))

# Groups 6,19,33,53,55,56,57 have 2 groups of sequence with 4+
# Remove sequence group with lower number of years (i.e., 4 or 5)
PMtracks_incomp_ts <- PMtracks_incomp_ts[!(PMtracks_incomp_ts$groupID == 6 & PMtracks_incomp_ts$length == 5 )]
PMtracks_incomp_ts <- PMtracks_incomp_ts[!(PMtracks_incomp_ts$groupID == 19 & PMtracks_incomp_ts$length == 4 )]
PMtracks_incomp_ts <- PMtracks_incomp_ts[!(PMtracks_incomp_ts$groupID == 33 & PMtracks_incomp_ts$length == 4 )]
PMtracks_incomp_ts <- PMtracks_incomp_ts[!(PMtracks_incomp_ts$groupID == 53 & PMtracks_incomp_ts$length == 4 )]
PMtracks_incomp_ts <- PMtracks_incomp_ts[!(PMtracks_incomp_ts$groupID == 55 & PMtracks_incomp_ts$length == 4 )]
PMtracks_incomp_ts <- PMtracks_incomp_ts[!(PMtracks_incomp_ts$groupID == 56 & PMtracks_incomp_ts$length == 4 )]
PMtracks_incomp_ts <- PMtracks_incomp_ts[!(PMtracks_incomp_ts$groupID == 57 & PMtracks_incomp_ts$length == 4 )]

### Calculate residuals PRCF[1] of groups with incomplete time series (n=45) 
# Make sure order on group and year is right before running the loop 
# Order on group and year 
# Order by year first
PMtracks_incomp_ts<-PMtracks_incomp_ts[order(PMtracks_incomp_ts$Year),]
# Order by groupID next
PMtracks_incomp_ts<-PMtracks_incomp_ts[order(PMtracks_incomp_ts$groupID),]

# Loop going through all groupID 
resi_incomp_ts_lag_1 <- data.frame() #Empty dataframe for residuals
for (group in unique(PMtracks_incomp_ts$groupID)) { # for each unique group 
  tmp <- PMtracks_incomp_ts[which(PMtracks_incomp_ts$groupID == group),] # temporary groupID for each loop
  tmptmp <-
    data.frame(lag = 1:(length(PRCF(tmp$PMIndex)) - 1), resid = NA) #fill lag column with num of lags and empty resid columnn 
  for (lag in 1:(length(PRCF(tmp$PMIndex)) - 1)) { # loop through each lag 
    tmptmp[lag, 'resid'] <- PRCF(tmp$PMIndex)[lag] # fill lag with respective resid from PRCF
  }
  tmptmp$abs_resid <- abs(tmptmp$resid) #convert to abs val. 
  output_lag_1 <-
    data.frame(resid = tmptmp[which(tmptmp$abs_resid == (abs(tmptmp$resid)[1])), 'resid'], #row for resid
               group = group, #row for groupID
               lag = tmptmp[which(tmptmp$abs_resid == (abs(tmptmp$resid)[1])), 'lag']) #row for lag of highest absolute resid 
  # Add residual, group and lag to dataframe 
  resi_incomp_ts_lag_1 <- rbind(resi_incomp_ts_lag_1, output_lag_1) # bind row with each of the group 
}
# Save residuals lag 1
saveRDS(resi_incomp_ts_lag_1, file = "resi_incomp_ts_lag_1.rds")

# 6. Plot PRCF residuals #######################################################
# Residuals for complete time series
for (group in unique(PMtracks_comp_ts$groupID)) {
  PRCF(PMtracks_comp_ts[which(PMtracks_comp_ts$groupID == group),]$PMIndex)
}

# Residuals for incomplete time series
for (group in unique(PMtracks_incomp_ts$groupID)) {
  PRCF(PMtracks_incomp_ts[which(PMtracks_incomp_ts$groupID == group),]$PMIndex)
}

