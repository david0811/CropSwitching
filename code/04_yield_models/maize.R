setwd(getwd())

library(tidyr)
library(dplyr)
library(rjags)
library(R2jags)
library(parallel)
library(coda)

source('./00_useful_R_funcs.R')

# Load yield data (used for all models!)
usda <- read.csv('../data/historical/yields/maize_all.csv')
usda <- subset(usda, select=c('fips','year','yield', 'state'))
usda <- filter(usda, year <= 2013) # Match GCM period
usda <- filter(usda, year >= 1950) # Match GCM period
usda <- merge(usda,  usda %>% count(fips)) # Filter to counties with >=50% coverage
usda <- filter(usda, n >= 30) # Filter to counties with >=50% coverage
usda$log_yield <- log(usda$yield)

# Climate data
loca.predictors <- read.csv('../data/historical/predictors/Livneh/Livneh_maize_H21_predictors.csv')
loca.covariates <- read.csv('../data/historical/predictors/Livneh/Livneh_covariates_all.csv')

# Remove counties with NA irrigation fraction (not producing maize circa 2010-2015 when dataset was made)
# Around 300 counties
# Need to fix later
loca.covariates <- drop_na(loca.covariates, maize_irr_frac)

# Merge
df <- merge(loca.predictors, usda)
df <- df[df$fips %in% loca.covariates$fips,]
df$year2 <- df$year **2

# Save fips codes
fips <- as.data.frame(unique(df$fips))
colnames(fips) <- 'fips'
write.csv(fips, file='../data/historical/crop_fips_historical/maize_fips.csv',
          row.names = FALSE)

# Input data for models
Y <- df %>% pivot_wider(names_from = fips, values_from = log_yield, id_cols = year)
GDD <- df %>% pivot_wider(names_from = fips, values_from = GDD, id_cols = year, values_fill=-1)
EDD <- df %>% pivot_wider(names_from = fips, values_from = EDD, id_cols = year, values_fill=-1.)
P <- df %>% pivot_wider(names_from = fips, values_from = prcp, id_cols = year, values_fill=-1.)
P2 <- df %>% pivot_wider(names_from = fips, values_from = prcp2, id_cols = year, values_fill=-1.)
SM <- df %>% pivot_wider(names_from = fips, values_from = SM_mean, id_cols = year, values_fill=-1.)
SM2 <- df %>% pivot_wider(names_from = fips, values_from = SM_mean2, id_cols = year, values_fill=-1.)
SM_min <- df %>% pivot_wider(names_from = fips, values_from = SM_week_min, id_cols = year, values_fill=-1.)
SM_max <- df %>% pivot_wider(names_from = fips, values_from = SM_week_max, id_cols = year, values_fill=-1.)

# Get covariates for modeled counties
bios <- c('bio_1', 'bio_3', 'bio_4', 'bio_12')
covariates <- subset(loca.covariates, select=append(bios, c('maize_irr_frac', 'fips')))
covariates <- covariates[covariates$fips %in% df$fips,]
covariates <- subset(covariates, select=append(bios, 'maize_irr_frac'))

# Dims
nyrs <- nrow(Y)
ncounties <- ncol(Y)-1
ncovs = ncol(covariates)

# Some arrangement # 
Technology = as.matrix(Y[,1]) - 1949
Y = as.matrix(Y[,1:ncounties+1])
GDD = as.matrix(GDD[,1:ncounties+1])
EDD = as.matrix(EDD[,1:ncounties+1])
P = as.matrix(P[,1:ncounties+1])
P2 = as.matrix(P2[,1:ncounties+1])
SM = as.matrix(SM[,1:ncounties+1])
SM2 = as.matrix(SM2[,1:ncounties+1])
SM_max = as.matrix(SM_max[,1:ncounties+1])
SM_min = as.matrix(SM_min[,1:ncounties+1])

########################### 
## (GDD, EDD, P ,P2)
########################### 
model.file <- "./04_yield_models/bugs_files/SR_09.txt"

jags.data <- list('Technology' = Technology,
                  'GDD' = GDD,
                  'EDD' = EDD,
                  'Y' = Y,
                  'P' = P,
                  'P2' = P2,
                  'covariates' = covariates,
                  'ncounties' = ncounties,
                  'nyrs'= nyrs,
                  'ncovs' = ncovs)

jags.params <- c('tau.y',
                 'alpha','beta1','beta2','beta3','beta4','beta5',
                 'tau.alpha','tau.beta1','tau.beta2','tau.beta3','tau.beta4','tau.beta5',
                 'a0','a1','a2','a3','a4','a5',
                 'b0','b1','b2','b3','b4','b5'
                 )

jags.inits <- NULL

how_many_samples(n_iter=10000, n_burnin=5000, n_thin=1, n_chains=3)

# Fit the model using run.jags -- parallel #
start_time <- Sys.time()
SR_09_model  <- jags.parallel(model.file = model.file, 
                                data = jags.data, 
                                parameters.to.save = jags.params,
                                inits = jags.inits,
                                n.iter=10000L,
                                n.burnin=5000L,
                                n.thin = 1L,
                                jags.module = c("glm","dic"),
                                n.chains = 3)
end_time <- Sys.time()
end_time - start_time

# Store
saveRDS(SR_09_model, file = "../data/simulated/yield_models/maize_SR_09.rds")

########################### 
## (GDD, EDD, SM, SM2)
########################### 
model.file <- "./04_yield_models/bugs_files/SM_ave.txt"

jags.data <- list('Technology' = Technology,
                  'GDD' = GDD,
                  'EDD' = EDD,
                  'Y' = Y,
                  'SM' = SM,
                  'SM2' = SM2,
                  'covariates' = covariates,
                  'ncounties' = ncounties,
                  'nyrs'= nyrs,
                  'ncovs' = ncovs)

jags.params <- c('tau.y',
                 'alpha','beta1','beta2','beta3','beta4','beta5',
                 'tau.alpha','tau.beta1','tau.beta2','tau.beta3','tau.beta4','tau.beta5',
                 'a0','a1','a2','a3','a4','a5',
                 'b0','b1','b2','b3','b4','b5'
)

jags.inits <- NULL

# Fit the model using run.jags -- parallel #
start_time <- Sys.time()
SM_ave_model  <- jags.parallel(model.file = model.file, 
                               data = jags.data, 
                               parameters.to.save = jags.params,
                               inits = jags.inits,
                               n.iter=10000L,
                               n.burnin=5000L,
                               n.thin = 1L,
                               jags.module = c("glm","dic"),
                               n.chains = 3)
end_time <- Sys.time()
end_time - start_time

# Store
saveRDS(SM_ave_model, file = "../data/simulated/yield_models/maize_SM_ave.rds")

###################################################### 
## (GDD, EDD, SM_weekly_min, SM_weekly_max)
###################################################### 
model.file <- "./04_yield_models/bugs_files/SM_week.txt"

jags.data <- list('Technology' = Technology,
                  'GDD' = GDD,
                  'EDD' = EDD,
                  'Y' = Y,
                  'SM_min' = SM_min,
                  'SM_max' = SM_max,
                  'covariates' = covariates,
                  'ncounties' = ncounties,
                  'nyrs'= nyrs,
                  'ncovs' = ncovs)

jags.params <- c('tau.y',
                 'alpha','beta1','beta2','beta3','beta4','beta5',
                 'tau.alpha','tau.beta1','tau.beta2','tau.beta3','tau.beta4','tau.beta5',
                 'a0','a1','a2','a3','a4','a5',
                 'b0','b1','b2','b3','b4','b5'
)

jags.inits <- NULL

# Fit the model using run.jags -- parallel #
start_time <- Sys.time()
SM_week_model  <- jags.parallel(model.file = model.file, 
                                data = jags.data, 
                                parameters.to.save = jags.params,
                                inits = jags.inits,
                                n.iter = 10000L,
                                n.burnin = 5000L,
                                n.thin = 1L,
                                jags.module = c("glm","dic"),
                                n.chains = 3)
end_time <- Sys.time()
end_time - start_time

# Store
saveRDS(SM_week_model, file = "../data/simulated/yield_models/maize_SM_week.rds")
