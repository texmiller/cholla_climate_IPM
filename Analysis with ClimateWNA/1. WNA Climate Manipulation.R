## Use this code to organize climate data for analysis
## 7/15/2018

setwd('C:/Users/Kevin Czachura/Dropbox/Kevin Czachura senior thesis/Analysis with ClimateWNA')

## read data ouput from Climate WNA
historical <- read.csv('Historical Climate WNA Data Output.csv')
future <- read.csv('Future Climate WNA Data Output.csv')

data <- rbind(historical, future)

library(data.table)

## Add new columns for averages of warm and cold periods
## Cold period runs from October Year_t to April Year_t+1

## Precipitation (Total precipitation of the period)
data$prcp_warm <- as.numeric(as.character(data$PPT05)) + 
  as.numeric(as.character(data$PPT06)) + 
  as.numeric(as.character(data$PPT07)) +
  as.numeric(as.character(data$PPT08)) + 
  as.numeric(as.character(data$PPT09))

data$prcp_cold <- as.numeric(as.character(data$PPT10)) + 
  as.numeric(as.character(data$PPT11)) + 
  as.numeric(as.character(data$PPT12))
+ as.numeric(as.character(shift(data$PPT01, n = 1, type = 'lead'))) 
+ as.numeric(as.character(shift(data$PPT02, n = 1, type = 'lead')))
+ as.numeric(as.character(shift(data$PPT03, n = 1, type = 'lead'))) 
+ as.numeric(as.character(shift(data$PPT04, n = 1, type = 'lead')))


## Mean Temperature (Average temperature of the period)
data$meantemp_warm <- (as.numeric(as.character(data$Tave05)) + 
                             as.numeric(as.character(data$Tave06)) 
                           + as.numeric(as.character(data$Tave07))
                           + as.numeric(as.character(data$Tave08)) + 
                             as.numeric(as.character(data$Tave09)))/5

data$meantemp_cold <- (as.numeric(as.character(data$Tave10)) + 
   as.numeric(as.character(data$Tave11)) 
 + as.numeric(as.character(data$Tave12))
 + as.numeric(as.character(shift(data$Tave01, n = 1, type = 'lead'))) 
 + as.numeric(as.character(shift(data$Tave02, n = 1, type = 'lead')))
 + as.numeric(as.character(shift(data$Tave03, n = 1, type = 'lead'))) 
 + as.numeric(as.character(shift(data$Tave04, n = 1, type = 'lead'))))/7


## Maximum Temperature (Maximum temperature of the period)
data$maxtemp_warm <- pmax(as.numeric(as.character(data$Tmax05)), 
                             as.numeric(as.character(data$Tmax06)),
                             as.numeric(as.character(data$Tmax07)),
                             as.numeric(as.character(data$Tmax08)),
                             as.numeric(as.character(data$Tmax09)))

data$maxtemp_cold <- pmax(as.numeric(as.character(data$Tmax10)), 
                             as.numeric(as.character(data$Tmax11)),
                             as.numeric(as.character(data$Tmax12)),
  as.numeric(as.character(shift(data$Tmax01, n = 1, type = 'lead'))), 
  as.numeric(as.character(shift(data$Tmax02, n = 1, type = 'lead'))),
  as.numeric(as.character(shift(data$Tmax03, n = 1, type = 'lead'))),
  as.numeric(as.character(shift(data$Tmax04, n = 1, type = 'lead'))))


## Minimum Temperature (Minimum Temperature of the period)
data$mintemp_warm <- pmin(as.numeric(as.character(data$Tmin05)), 
                             as.numeric(as.character(data$Tmin06)),
                             as.numeric(as.character(data$Tmin07)),
                             as.numeric(as.character(data$Tmin08)),
                             as.numeric(as.character(data$Tmin09)))

data$mintemp_cold <- pmin(as.numeric(as.character(data$Tmin10)), 
                             as.numeric(as.character(data$Tmin11)),
                             as.numeric(as.character(data$Tmin12)),
 as.numeric(as.character(shift(data$Tmin01, n = 1, type = 'lead'))), 
 as.numeric(as.character(shift(data$Tmin02, n = 1, type = 'lead'))),
 as.numeric(as.character(shift(data$Tmin03, n = 1, type = 'lead'))),
 as.numeric(as.character(shift(data$Tmin04, n = 1, type = 'lead'))))


## Create a new dataframe to hold the data we want

climate <- data.frame(data$Year,
                      data$meantemp_warm,
                      data$meantemp_cold,
                      data$maxtemp_warm,
                      data$maxtemp_cold,
                      data$mintemp_warm,
                      data$mintemp_cold,
                      data$prcp_warm,
                      data$prcp_cold)

colnames(climate) <- c('Year_t', 'Mean_Temp_Warm', 'Mean_Temp_Cold',                           'Max_Temp_Warm', 'Max_Temp_Cold',
                          'Min_Temp_Warm', 'Min_Temp_Cold',
                          'Tot_Precip_Warm', 'Tot_Precip_Cold')
climate <- subset(climate, Year_t < 2100)
                            
## write to CSV
write.csv(climate, 'Climate WNA.csv')

