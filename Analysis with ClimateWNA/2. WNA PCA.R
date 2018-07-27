## Use the newly organized climate data for the PCA
## 7/16/2018

setwd('C:/Users/Kevin Czachura/Dropbox/Kevin Czachura senior thesis/Analysis with ClimateWNA')

climate <- read.csv('Climate WNA.csv')


## cannot log transform the continuous variables due to negative values
clim.vars <- climate[, 3:10]
## use year as the categorical variable
clim.years <- climate[, 2]

## apply a PCA with scale. and center set to TRUE
clim.pca <- prcomp(clim.vars, center = TRUE, scale. = TRUE)
## find the weight of each climate variable in the PCs
clim.pca$rotation

## write the PC values to a csv file

Year_t <- 1901:2099
clim.pca$x

pca.dat <- data.frame(Year_t, clim.pca$x)
write.csv(pca.dat, 'ClimateWNA PC Values.csv')




## ANALYZE THE RESULTS




## print method
print(clim.pca)

## plot method
plot(clim.pca, type = 'l')


## summary method
summary(clim.pca)

## predict principle components
predict(clim.pca, newdata = tail(clim.vars, 2))
