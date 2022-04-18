setwd(getwd())

#library(ggplot2)
#library(tidyr)
#library(dplyr)
library(foreach)
library(doParallel)
library(data.table)

#############################################
# Function to get yield projections
#############################################
yield_project_NOTREND <- function(crop, yield.model, county, npar, nthin, nsim, ncores, models, fips){
  # Thinning of parameters
  l <-seq(1, npar, nthin)
  
  # Read JAGS model
  yield_model <- readRDS(paste('../data/simulated/yield_models/', crop, '_', yield.model, '.rds', sep=''))
  tau.y <- yield_model$BUGSoutput$sims.list$tau.y
  
  # Need to use level 1 or level 2 parameters
  if (is.element(county, as.matrix(fips))){
    alpha <- yield_model$BUGSoutput$sims.list$alpha
    beta1 <- yield_model$BUGSoutput$sims.list$beta1
    beta2 <- yield_model$BUGSoutput$sims.list$beta2
    beta3 <- yield_model$BUGSoutput$sims.list$beta3
    beta4 <- yield_model$BUGSoutput$sims.list$beta4
    beta5 <- yield_model$BUGSoutput$sims.list$beta5
    rm(yield_model) 
  
    # Get list of means, sampling from every nthin-th parameter set
    county.id <- which(fips == county)
    Params <- matrix(c(beta1[l,county.id], beta2[l,county.id], beta3[l,county.id],
                       beta4[l,county.id], beta5[l,county.id]), nrow=npar/nthin, ncol=5)
    alpha <- alpha[l,county.id]
    
  } else{
    a0 <- yield_model$BUGSoutput$sims.list$a0
    a1 <- yield_model$BUGSoutput$sims.list$a1
    a2 <- yield_model$BUGSoutput$sims.list$a2
    a3 <- yield_model$BUGSoutput$sims.list$a3
    a4 <- yield_model$BUGSoutput$sims.list$a4
    a5 <- yield_model$BUGSoutput$sims.list$a5
    
    b0 <- yield_model$BUGSoutput$sims.list$b0
    b1 <- yield_model$BUGSoutput$sims.list$b1
    b2 <- yield_model$BUGSoutput$sims.list$b2
    b3 <- yield_model$BUGSoutput$sims.list$b3
    b4 <- yield_model$BUGSoutput$sims.list$b4
    b5 <- yield_model$BUGSoutput$sims.list$b5
    
    tau.alpha <- yield_model$BUGSoutput$sims.list$tau.alpha
    tau.beta1 <- yield_model$BUGSoutput$sims.list$tau.beta1
    tau.beta2 <- yield_model$BUGSoutput$sims.list$tau.beta2
    tau.beta3 <- yield_model$BUGSoutput$sims.list$tau.beta3
    tau.beta4 <- yield_model$BUGSoutput$sims.list$tau.beta4
    tau.beta5 <- yield_model$BUGSoutput$sims.list$tau.beta5
    rm(yield_model) 
    
    # Covariates
    covariates <- read.csv('../data/historical/predictors/Livneh/Livneh_covariates_all.csv')
    covariates <- covariates[covariates$fips == countyfips,]
    bios <- c('bio_1', 'bio_3', 'bio_4', 'bio_12')
    
    # if no irrigation data use average of other crops in county
    if (is.na(subset(covariates, select=paste(crop, '_irr_frac', sep="")))[1]){
      covariates$irr_frac <- mean(as.matrix(covariates[,2:9]), na.rm = TRUE)
      covariates <- subset(covariates, select=append(bios, 'irr_frac'))
    } else {
      covariates <- subset(covariates, select=append(bios, paste(crop, '_irr_frac', sep="")))  
    }
    
    # Get list of means, sampling from every nthin-th parameter set
    Params <- matrix(c(rnorm(length(l), a1[l] + b1[l,] %*% t(covariates), 1/tau.beta1[l,]), 
                       rnorm(length(l), a2[l] + b2[l,] %*% t(covariates), 1/tau.beta2[l,]),
                       rnorm(length(l), a3[l] + b3[l,] %*% t(covariates), 1/tau.beta3[l,]),
                       rnorm(length(l), a4[l] + b4[l,] %*% t(covariates), 1/tau.beta4[l,]), 
                       rnorm(length(l), a5[l] + b5[l,] %*% t(covariates), 1/tau.beta5[l,])), 
                     nrow=length(l), ncol=5)
    alpha <- rnorm(length(l), a0[l] + b0[l,] %*% t(covariates), 1/tau.alpha[l,])
  }
  
  # Get predictors
  if (yield.model == 'SM_week') {
    predictors <- c('SM_week_min', 'SM_week_max')
  }
  else if (yield.model == 'SR_09') {
    predictors <- c('prcp', 'prcp2')
  }
  else if (yield.model == 'SM_ave') {
    predictors <- c('SM_mean', 'SM_mean2')
  }
  
  # Set up cluster
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  # Get projections
  Ypred <- foreach(i=1:nmodels, .combine='c') %:%
    foreach(k=1:nyrs) %dopar% {
      # Subset climate data for county, year, predictors
      if (crop == 'sorghum'){
        clim <- read.csv(paste('../data/future/predictors/soy/', models[i], 
                               '_soy_H21_RCP45_30-50.csv', sep=""))
      }
      else{
        clim <- read.csv(paste('../data/future/predictors/', crop, '/', models[i], '_', 
                               crop, '_H21_RCP45_30-50.csv', sep=""))
      }
      clim <- clim[clim$fips == county,]
      clim <- subset(clim, select=c('fips', 'year', 'GDD', 'EDD', predictors))
      clim <- clim[k, 2:6]
      clim$year <- 2020 - 1949
      
      # Simulate future yield for given year and county, nsim times for each parameter set
      Ymean <- Params %*% t(clim) + alpha
      log_yield_sim <- rnorm(nsim * npar/nthin, Ymean, 1/tau.y[l])
      
      # Make dataframe and append
      Y <- data.frame(log_yield_sim = log_yield_sim,
                      fips = county,
                      year = 2029 + k,
                      clim_model = i,
                      param = seq(1,npar/nthin))
      # 'return'
      Y 
    }
  
  # Stop cluster
  stopCluster(cl)
  
  # Save single dataframe
  Ypred.df <- rbindlist(Ypred)
  path <- paste('../data/future/yields/', crop, '/', crop, '_', county,
                '_RCP45_30-50_', yield.model,'_NOTREND.csv', sep='')
  fwrite(Ypred.df, path)
}

# Climate models
models <- c('ACCESS1-0','ACCESS1-3','CCSM4','CESM1-BGC','CESM1-CAM5','CMCC-CMS',
            'CMCC-CM','CNRM-CM5','CSIRO-Mk3-6-0','CanESM2','EC-EARTH','FGOALS-g2',
            'GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-H','GISS-E2-R','HadGEM2-AO',
            'HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','MIROC-ESM-CHEM',
            'MIROC-ESM','MIROC5','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M',
            'inmcm4')

# Dimensions for single county
nyrs <- 1 + 2050 - 2030
nmodels <- length(models)
npar <- 15000
nthin <- 15
nsim <- 10
# How many simulations per year?
log(npar * nsim / nthin, base = 10)
# How many total simulations?
log(nyrs * nmodels * npar * nsim / nthin, base = 10)

# Historical counties
maize_fips <- read.csv('../data/historical/crop_fips_historical/maize_fips.csv')
soy_fips <- read.csv('../data/historical/crop_fips_historical/soy_fips.csv')
sorghum_fips <- read.csv('../data/historical/crop_fips_historical/sorghum_fips.csv')
cotton_fips <- read.csv('../data/historical/crop_fips_historical/cotton_fips.csv')
springwheat_fips <- read.csv('../data/historical/crop_fips_historical/springwheat_fips.csv')
barley_fips <- read.csv('../data/historical/crop_fips_historical/barley_fips.csv')

## Choose county
countyfips <- 17019
start_time <- Sys.time()
yield_project_NOTREND('maize', 'SM_week', countyfips, npar, nthin, nsim, 14, models, maize_fips)
yield_project_NOTREND('maize', 'SR_09', countyfips, npar, nthin, nsim, 14, models, maize_fips)
yield_project_NOTREND('maize', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, maize_fips)

yield_project_NOTREND('soy', 'SM_week', countyfips, npar, nthin, nsim, 14, models, soy_fips)
yield_project_NOTREND('soy', 'SR_09', countyfips, npar, nthin, nsim, 14, models, soy_fips)
yield_project_NOTREND('soy', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, soy_fips)

yield_project_NOTREND('sorghum', 'SM_week', countyfips, npar, nthin, nsim, 14, models, sorghum_fips)
yield_project_NOTREND('sorghum', 'SR_09', countyfips, npar, nthin, nsim, 14, models, sorghum_fips)
yield_project_NOTREND('sorghum', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, sorghum_fips)

yield_project_NOTREND('cotton', 'SM_week', countyfips, npar, nthin, nsim, 14, models, cotton_fips)
yield_project_NOTREND('cotton', 'SR_09', countyfips, npar, nthin, nsim, 14, models, cotton_fips)
yield_project_NOTREND('cotton', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, cotton_fips)

yield_project_NOTREND('barley', 'SM_week', countyfips, npar, nthin, nsim, 14, models, barley_fips)
yield_project_NOTREND('barley', 'SR_09', countyfips, npar, nthin, nsim, 14, models, barley_fips)
yield_project_NOTREND('barley', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, barley_fips)

yield_project_NOTREND('springwheat', 'SM_week', countyfips, npar, nthin, nsim, 14, models, springwheat_fips)
yield_project_NOTREND('springwheat', 'SR_09', countyfips, npar, nthin, nsim, 14, models, springwheat_fips)
yield_project_NOTREND('springwheat', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, springwheat_fips)

print(Sys.time() - start_time)

countyfips <- 40039
start_time <- Sys.time()
yield_project_NOTREND('maize', 'SM_week', countyfips, npar, nthin, nsim, 14, models, maize_fips)
yield_project_NOTREND('maize', 'SR_09', countyfips, npar, nthin, nsim, 14, models, maize_fips)
yield_project_NOTREND('maize', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, maize_fips)

yield_project_NOTREND('soy', 'SM_week', countyfips, npar, nthin, nsim, 14, models, soy_fips)
yield_project_NOTREND('soy', 'SR_09', countyfips, npar, nthin, nsim, 14, models, soy_fips)
yield_project_NOTREND('soy', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, soy_fips)

yield_project_NOTREND('sorghum', 'SM_week', countyfips, npar, nthin, nsim, 14, models, sorghum_fips)
yield_project_NOTREND('sorghum', 'SR_09', countyfips, npar, nthin, nsim, 14, models, sorghum_fips)
yield_project_NOTREND('sorghum', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, sorghum_fips)

yield_project_NOTREND('cotton', 'SM_week', countyfips, npar, nthin, nsim, 14, models, cotton_fips)
yield_project_NOTREND('cotton', 'SR_09', countyfips, npar, nthin, nsim, 14, models, cotton_fips)
yield_project_NOTREND('cotton', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, cotton_fips)

yield_project_NOTREND('barley', 'SM_week', countyfips, npar, nthin, nsim, 14, models, barley_fips)
yield_project_NOTREND('barley', 'SR_09', countyfips, npar, nthin, nsim, 14, models, barley_fips)
yield_project_NOTREND('barley', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, barley_fips)

yield_project_NOTREND('springwheat', 'SM_week', countyfips, npar, nthin, nsim, 14, models, springwheat_fips)
yield_project_NOTREND('springwheat', 'SR_09', countyfips, npar, nthin, nsim, 14, models, springwheat_fips)
yield_project_NOTREND('springwheat', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, springwheat_fips)

print(Sys.time() - start_time)

countyfips <- 38071
start_time <- Sys.time()
yield_project_NOTREND('maize', 'SM_week', countyfips, npar, nthin, nsim, 14, models, maize_fips)
yield_project_NOTREND('maize', 'SR_09', countyfips, npar, nthin, nsim, 14, models, maize_fips)
yield_project_NOTREND('maize', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, maize_fips)

yield_project_NOTREND('soy', 'SM_week', countyfips, npar, nthin, nsim, 14, models, soy_fips)
yield_project_NOTREND('soy', 'SR_09', countyfips, npar, nthin, nsim, 14, models, soy_fips)
yield_project_NOTREND('soy', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, soy_fips)

yield_project_NOTREND('sorghum', 'SM_week', countyfips, npar, nthin, nsim, 14, models, sorghum_fips)
yield_project_NOTREND('sorghum', 'SR_09', countyfips, npar, nthin, nsim, 14, models, sorghum_fips)
yield_project_NOTREND('sorghum', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, sorghum_fips)

yield_project_NOTREND('cotton', 'SM_week', countyfips, npar, nthin, nsim, 14, models, cotton_fips)
yield_project_NOTREND('cotton', 'SR_09', countyfips, npar, nthin, nsim, 14, models, cotton_fips)
yield_project_NOTREND('cotton', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, cotton_fips)

yield_project_NOTREND('barley', 'SM_week', countyfips, npar, nthin, nsim, 14, models, barley_fips)
yield_project_NOTREND('barley', 'SR_09', countyfips, npar, nthin, nsim, 14, models, barley_fips)
yield_project_NOTREND('barley', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, barley_fips)

yield_project_NOTREND('springwheat', 'SM_week', countyfips, npar, nthin, nsim, 14, models, springwheat_fips)
yield_project_NOTREND('springwheat', 'SR_09', countyfips, npar, nthin, nsim, 14, models, springwheat_fips)
yield_project_NOTREND('springwheat', 'SM_ave', countyfips, npar, nthin, nsim, 14, models, springwheat_fips)

print(Sys.time() - start_time)