## Selecting which models to use 

## NOTE: Why are th emodel selections so strange?

# 7/16/2018

setwd('C:/Users/Kevin Czachura/Dropbox/Kevin Czachura senior thesis/Analysis with ClimateWNA')

library(xlsx)
library(LaplacesDemon)
library(AICcmodavg)
library(lme4)

## prepare the data
climate <- read.csv('Climate WNA.csv')
demography <- read.csv('ChollaDemographyCSV.csv')
prin.comp <- read.csv('ClimateWNA PC Values.csv')

## merge the data frames by year
data <- merge(climate, demography, by = 'Year_t', all = TRUE)
data <- merge(data, prin.comp, by = 'Year_t', all = TRUE)

## write function to calculate volume of cone
volume <- function(h, w, p){
  {(1/3)*pi*h*(((w + p)/2)/2)^2}
}

## add volume variables to data frame
data$volume_t <- volume(h = data$Height_t, w = data$Width_t,
                        p = data$Perp_t)
data$volume_t1 <- volume(h = data$Height_t1, w = data$Width_t1,
                         p = data$Perp_t1)

## take the log of the volumes
data$logvol_t <- log(data$volume_t)
data$logvol_t1 <- log(data$volume_t1)

## create a field for a binomial 'Flowering' category
data$Flowering <- ifelse(data$Goodbuds_t1 > 0, 1, 0)


###############################################################################


## FLOWERING MODEL SELECTION

# Model 1: logvol_t1 ~ logvol_t + (logvol_t|Year_t)
mod1 <- glmer(Flowering ~ logvol_t + (logvol_t|Year_t), data = data,
              family = 'binomial')

# Model 2: logvol_t1 ~ logvol_t + PC1 + PC2
mod2 <- glmer(Flowering ~ logvol_t + PC1 + PC2 + (logvol_t|Year_t), 
              data = data, family = 'binomial')

# Model 3: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3
mod3 <- glmer(Flowering ~ logvol_t + PC1 + PC2 + PC3 + (logvol_t|Year_t), 
              data = data, family = 'binomial')

# Model 4: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4
mod4 <- glmer(Flowering ~ logvol_t + PC1 + PC2 + PC3 + PC4
              + (logvol_t|Year_t), data = data, family = 'binomial')

# Model 5: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5
mod5 <- glmer(Flowering ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5
              + (logvol_t|Year_t), data = data, family = 'binomial')

# Model 6: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6
mod6 <- glmer(Flowering ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6
              + (logvol_t|Year_t), data = data, family = 'binomial')

# Model 7: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7
mod7 <- glmer(Flowering ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7
              + (logvol_t|Year_t), data = data, family = 'binomial')


# Model 8: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8
mod8 <- glmer(Flowering ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8
              + (logvol_t|Year_t), data = data, family = 'binomial')

# create an object for the candidate models
cand.models <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8)
mod.names <- paste("mod", 1:length(cand.models), sep = "")

# compare the AIC values
aictab(cand.set = cand.models, modnames = mod.names, sort = TRUE)

# the best candidate model has the lowest AICc score


##############################################################################


## FLOWER PRODUCTION MODEL SELECTION

# Model 1: logvol_t1 ~ logvol_t + (logvol_t|Year_t)
mod1 <- glmer(Goodbuds_t1 ~ logvol_t + (logvol_t|Year_t), data = data,
              family = 'poisson')

# Model 2: logvol_t1 ~ logvol_t + PC1 + PC2
mod2 <- glmer(Goodbuds_t1 ~ logvol_t + PC1 + PC2 + (logvol_t|Year_t), data = data, 
              family = 'poisson')

# Model 3: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3
mod3 <- glmer(Goodbuds_t1 ~ logvol_t + PC1 + PC2 + PC3 + (logvol_t|Year_t), 
              data = data, family = 'poisson')

# Model 4: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4
mod4 <- glmer(Goodbuds_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4
              + (logvol_t|Year_t), data = data, family = 'poisson')

# Model 5: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5
mod5 <- glmer(Goodbuds_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5
              + (logvol_t|Year_t), data = data, family = 'poisson')

# Model 6: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6
mod6 <- glmer(Goodbuds_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6
              + (logvol_t|Year_t), data = data, family = 'poisson')

# Model 7: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7
mod7 <- glmer(Goodbuds_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7
              + (logvol_t|Year_t), data = data, family = 'poisson')

# Model 8: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8
mod8 <- glmer(Goodbuds_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8
              + (logvol_t|Year_t), data = data, family = 'poisson')

# create an object for the candidate models
cand.models <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8)
mod.names <- paste("mod", 1:length(cand.models), sep = "")

# compare the AIC values
aictab(cand.set = cand.models, modnames = mod.names, sort = TRUE)

# the best candidate model has the lowest AICc score


################################################################################


## GROWTH MODEL SELECTION

# Model 1: logvol_t1 ~ logvol_t + (logvol_t|Year_t)
mod1 <- lmer(logvol_t1 ~ logvol_t + (logvol_t|Year_t), data = data, REML=F)

# Model 2: logvol_t1 ~ logvol_t + PC1 + PC2
mod2 <- lmer(logvol_t1 ~ logvol_t + PC1 + PC2 + (logvol_t|Year_t), data = data, 
             REML=F)

# Model 3: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3
mod3 <- lmer(logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + (logvol_t|Year_t), 
             data = data, REML=F)

# Model 4: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4
mod4 <- lmer(logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + (logvol_t|Year_t), 
             data = data, REML=F)

# Model 5: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5
mod5 <- lmer(logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5
              + (logvol_t|Year_t), data = data, REML=F)

# Model 6: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6
mod6 <- lmer(logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6
             + (logvol_t|Year_t), data = data, REML=F)

# Model 7: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7
mod7 <- lmer(logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7
             + (logvol_t|Year_t), data = data, REML=F)

# Model 8: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8
mod8 <- lmer(logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8
             + (logvol_t|Year_t), data = data, REML=F)

# create an object for the candidate models
cand.models <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8)
mod.names <- paste("mod", 1:length(cand.models), sep = "")

# compare the AIC values
aictab(cand.set = cand.models, modnames = mod.names, sort = TRUE)

# the best candidate model has the lowest AICc score


################################################################################


## SURVIVAL MODEL SELECTION

# Model 1: logvol_t1 ~ logvol_t + (logvol_t|Year_t)
mod1 <- glmer(Survival_t1 ~ logvol_t + (logvol_t|Year_t), data = data,
              family = 'binomial')

# Model 2: logvol_t1 ~ logvol_t + PC1 + PC2
mod2 <- glmer(Survival_t1 ~ logvol_t + PC1 + PC2 + (logvol_t|Year_t), 
              data = data, family = 'binomial')

# Model 3: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3
mod3 <- glmer(Survival_t1 ~ logvol_t + PC1 + PC2 + PC3 + (logvol_t|Year_t), 
              data = data, family = 'binomial')

# Model 4: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4
mod4 <- glmer(Survival_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4
              + (logvol_t|Year_t), data = data, family = 'binomial')

# Model 5: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5
mod5 <- glmer(Survival_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5
              + (logvol_t|Year_t), data = data, family = 'binomial')

# Model 6: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6
mod6 <- glmer(Survival_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6
              + (logvol_t|Year_t), data = data, family = 'binomial')

# Model 7: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7
mod7 <- glmer(Survival_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7
              + (logvol_t|Year_t), data = data, family = 'binomial')


# Model 8: logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8
mod8 <- glmer(Survival_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8
              + (logvol_t|Year_t), data = data, family = 'binomial')

# create an object for the candidate models
cand.models <- list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8)
mod.names <- paste("mod", 1:length(cand.models), sep = "")

# compare the AIC values
aictab(cand.set = cand.models, modnames = mod.names, sort = TRUE)

# the best candidate model has the lowest AICc score

