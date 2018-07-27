## Script to create a figure for the lambda values
## 7/16/2018

setwd('C:/Users/Kevin Czachura/Dropbox/Kevin Czachura senior thesis/Analysis with ClimateWNA')

library(LaplacesDemon)
library(xlsx)

## read in the data
## prepare the data
climate <- read.csv('Climate WNA.csv')
demography <- read.csv('ChollaDemographyCSV.csv')
prin.comp <- read.csv('ClimateWNA PC Values.csv')
lambda <- read.csv('ClimateWNA Lambda Values.csv')

dat <- merge(climate, lambda, by = 'Year_t')
dat <- merge(dat, prin.comp, by = 'Year_t')



## Make a plot
plot(dat$Year_t, dat$lambda, type="n", xlim=c(1901, 2099), ylim=
       c(min(dat$Lambda)-0.15, max(dat$Lambda)+0.15),
     ylab = 'Population Growth Rate', xlab = 'Year t', cex.lab = 1.5)

## Backcasted years
lines(1901:2004, dat$Lambda[1:104], lty = 1, lwd = 2, col = 'red')

## Years for demography data
lines(2004:2017, dat$Lambda[104:117], lty = 1, lwd = 2, col = 'green')

## Forecasted years
lines(2018:2050, dat$Lambda[118:199], lty = 1, lwd = 2, col = 'blue')


abline(1, 0, lty = 1, col = 'grey42', lwd = 1.5)
title(outer=FALSE, adj=0.025, main="Population Growth Rate", cex.main=2,
      col="black",font=2,line=0.75)

legend('topleft', legend = c('Backcasted', 'Demography Years'),
       lty = c(1, 1, 1), lwd = c(2, 2, 2), 
       col = c('red', 'blue', 'green'), 
       bty = 'n')



