## IPM Analysis using Climate WNA


## NOTE: the models and params and function must all be manually changed
##       each time it is run with different models

## 7/16/2018

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

## subset a dataframe for new recruits and combine into one size column
recruits <- subset(data, Recruit == 1)
recruits$Year_t <- recruits$Year_t1 - 1
recruits$logvol_t[is.na(recruits$logvol_t)] <- 0
recruits$logvol_t1[is.na(recruits$logvol_t1)] <- 0
recruits$recvol <- recruits$logvol_t + recruits$logvol_t1
recruits <- subset(recruits, recvol != 0)


## Models
################################################################################

library(lme4)

# Growth Model
growth.mod <- lmer(logvol_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + PC5 + 
                     PC6 + PC7 + PC8 +
                     (logvol_t|Year_t), data = data, REML=F)

# Flower Production Model
flowerprod.mod <- glmer(Goodbuds_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + 
                          PC5
                        + PC6 + 
                          PC7 + PC8 + (logvol_t|Year_t), data = data, 
                        family = 'poisson')

# Survival Model
survival.mod <- glmer(Survival_t1 ~ logvol_t + PC1 + PC2 + PC3 + PC4 + 
                        PC5 
                      #+ PC6 + PC7 + PC8
                      + (logvol_t|Year_t), 
                      data = data, family = 'binomial')

# Flowering Model
flowering.mod <- glmer(Flowering ~ logvol_t + PC1 + PC2 + PC3 + PC4 + 
                         PC5
                       + PC6 + PC7 + 
                         #PC8 + 
                         (logvol_t|Year_t), data = data, 
                       family = 'binomial')

################################################################################

# create a vector of all of the variables needed for the IPM

params <- vector("numeric", length = 45)

params[1] <- fixef(growth.mod)[1] # growth intercept
params[2] <- fixef(growth.mod)[2] # growth size t
params[3] <- fixef(growth.mod)[3] # growth PC1
params[4] <- fixef(growth.mod)[4] # growth PC2
params[5] <- fixef(growth.mod)[5] # growth PC3
params[6] <- fixef(growth.mod)[6] # growth PC4
params[7] <- fixef(growth.mod)[7] # growth PC5
params[8] <- fixef(growth.mod)[8] # growth PC6
params[9] <- fixef(growth.mod)[9] # growth PC7
params[10] <- fixef(growth.mod)[10] # growth PC8
params[11] <- var(resid(growth.mod)) # growth variance of residuals
params[12] <- fixef(flowerprod.mod)[1] # flower production intercept
params[13] <- fixef(flowerprod.mod)[2] # flower production size t
params[14] <- fixef(flowerprod.mod)[3] # flower production PC1
params[15] <- fixef(flowerprod.mod)[4] # flower production PC2
params[16] <- fixef(flowerprod.mod)[5] # flower production PC3
params[17] <- fixef(flowerprod.mod)[6] # flower production PC4
params[18] <- fixef(flowerprod.mod)[7] # flower production PC5
params[19] <- fixef(flowerprod.mod)[8] # flower production PC6
params[20] <- fixef(flowerprod.mod)[9] # flower production PC7
params[21] <- fixef(flowerprod.mod)[10] # flower production PC8
params[22] <- fixef(survival.mod)[1] # survival intercept
params[23] <- fixef(survival.mod)[2] # survival size t
params[24] <- fixef(survival.mod)[3] # survival PC1
params[25] <- fixef(survival.mod)[4] # survival PC2
params[26] <- fixef(survival.mod)[5] # survival PC3
params[27] <- fixef(survival.mod)[6] # survival PC4
params[28] <- fixef(survival.mod)[7] # survival PC5
params[29] <- fixef(survival.mod)[8] # survival PC6
params[30] <- fixef(survival.mod)[9] # survival PC7
params[31] <- fixef(survival.mod)[10] # survival PC8
params[32] <- fixef(flowering.mod)[1] # flowering intercept
params[33] <- fixef(flowering.mod)[2] # flowering size t
params[34] <- fixef(flowering.mod)[3] # flowering PC1
params[35] <- fixef(flowering.mod)[4] # flowering PC2
params[36] <- fixef(flowering.mod)[5] # flowering PC3
params[37] <- fixef(flowering.mod)[6] # flowering PC4
params[38] <- fixef(flowering.mod)[7] # flowering PC5
params[39] <- fixef(flowering.mod)[8] # flowering PC6
params[40] <- fixef(flowering.mod)[9] # flowering PC7
params[41] <- fixef(flowering.mod)[10] # flowering PC8
params[42] <- mean(recruits$recvol) # mean recruit size
params[43] <- var(recruits$recvol) # recruit size variance
params[44] <- 0.0001 # establishment probability
params[45] <- 200 # seeds per bud


## Replace NAs with zeroes so they fall out of the equations
params[is.na(params)] <- 0



################################################################################

# write functions to collect the components into a projection kernel

# SURVIVAL OF X-SIZED PLANTS
sx <- function(x, params, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8){
  u <- exp(params[22] + params[23]*x + params[24]*pc1 + params[25]*pc2 +
             params[26]*pc3 + params[27]*pc4 + params[28]*pc5
           + params[29]*pc6 + params[30]*pc7 + params[31]*pc8)
  return(u/(1 + u));
}

# GROWTH FROM SIZE X TO Y
gxy <- function(x, y, params, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8){
  mux <- params[1] + params[2]*x + params[3]*pc1 + params[4]*pc2 +
    params[5]*pc3 + params[6]*pc4 + params[7]*pc5 + params[8]*pc6 +
    params[9]*pc7 + params[10]*pc8;
  sigmax2 <- params[11];
  sigmax <- sqrt(sigmax2);
  fac1 <- sqrt(2*pi)*sigmax; ## these lines are distribing sizes normally around the mean
  fac2 <- ((y - mux)^2)/(2*sigmax2); 
  return(exp(-fac2)/fac1); ## the function returns the probability density of the normal distribution
}

# GROWTH*SURVIVAL
# this step makes a new function that combines the two previous ones
# the pxy function will govern growth from x to y, conditioned on survival at size x
pxy <- function(x, y, params, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8){
  return(sx(x, params, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8)*gxy(x, y, 
                    params, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8))
}

# PRODUCTION OF Y-SIZED INDIVIDUALS FROM X-SIZED ADULTS
fxy <- function(x, y, params, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8){
  pr_flower <- invlogit(params[32] + params[33]*x + params[34]*pc1 + 
                          params[35]*pc2 + params[36]*pc3 + 
                          params[37]*pc4 + 
                          params[38]*pc5 + params[39]*pc6 + 
                          params[40]*pc7
                        + params[41]*pc8);
  flowerbuds <- exp(params[12] + params[13]*x + params[14]*pc1 + 
                      params[15]*pc2 + params[16]*pc3 + params[17]*pc4 + 
                      params[18]*pc5 + params[19]*pc6 + params[20]*pc7 + 
                      params[21]*pc8);
  seeds_per_bud <- params[45];
  sdlgsize.mean <- params[42];
  sdlgsize.var <- params[43];
  fac1 <- sqrt(2*pi)*sqrt(sdlgsize.var);
  fac2 <- ((y - sdlgsize.mean)^2)/(2*sdlgsize.var);
  f <- params[44]*pr_flower*flowerbuds*seeds_per_bud*(exp(-fac2)/fac1);
  return(f); ## again, it involves distributing new plants around a normal size distribution
}


################################################################################

# Put the functions together to make the projection kernel

n.big.matrix <- 50# 200 # this is the size of the matrix (nxn); value doesn't matter much as long as its big
minsize <- min(data$logvol_t, na.rm = T)	# smallest log(size) in the data set
maxsize <- max(data$logvol_t, na.rm = T)	# largest log(size) in the data set

# now the "big" matrix (= discrete version of the continuous kernel, with very fine slices)
bigmatrix <- function(n, params, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8){
  L <- minsize; U <- maxsize; # these are the upper and lower integration limits
  b <- L + c(0:n)*(U - L)/n;
  y <- 0.5*(b[1:n] + b[2:(n + 1)]);
  # these steps set the boundary points (b) and mesh points (y); basically the skeleton of the matrix
  # now construct the matrix
  I <- diag(n);
  ## the 'P' and 'B' steps are where all the demographic info comes in
  P <- t(outer(y, y, pxy, params = params, pc1 = pc1, pc2 = pc2, pc3 = pc3,
               pc4 = pc4, pc5 = pc5, pc6 = pc6, pc7 = pc7, pc8 = pc8))	
  B <- t(outer(y, y, fxy, params = params, pc1 = pc1, pc2 = pc2, pc3 = pc3,
               pc4 = pc4, pc5 = pc5, pc6 = pc6, pc7 = pc7, pc8 = pc8))
  M <- P + B
  K <- M;
  M <- (U - L)*M/n;
  P <- (U - L)*P/n;
  B <- (U - L)*B/n;
  return(list(matrix = M, kernel = K, meshpts = y, Pmatrix = P,
              Bmatrix = B, Imatrix = I));				 
  ## note that this function is returning a bunch of things
  ## we care mostly about the kernel K and the discrete matrix M
}



## Get the lambda values
pc.1 <- prin.comp$PC1
pc.2 <- prin.comp$PC2
pc.3 <- prin.comp$PC3
pc.4 <- prin.comp$PC4
pc.5 <- prin.comp$PC5
pc.6 <- prin.comp$PC6
pc.7 <- prin.comp$PC7
pc.8 <- prin.comp$PC8


myfunction <- function(n, params, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8){
  Re(eigen(bigmatrix(n, params, pc1, pc2, pc3, pc4, pc5, pc6, pc7, 
                     pc8)$matrix)$values[1])
}

lambda <- vector(mode = 'numeric', length = 199)

for (i in 1:199){
  lambda[i] <- myfunction(n.big.matrix, params, pc.1[i], pc.2[i], pc.3[i],
                          pc.4[i], pc.5[i], pc.6[i], pc.7[i], pc.8[i])
}

lambda.yr <- c(1901:2099)

pop.growth <- data.frame(lambda.yr, lambda)
colnames(pop.growth) <- c('Year_t', 'Lambda')

# write to csv
write.csv(pop.growth, 'ClimateWNA lambda values.csv')


