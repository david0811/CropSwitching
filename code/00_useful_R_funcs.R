########################################################
### Cross validation function returning average RMSE
########################################################
crossvalidate <- function(data, k, fml){
  # 'data' is the training set with the ".folds" column
  # 'k' is the number of folds we have
  # 'model' is a string describing a linear regression model formula
  
  # Initialize empty list for recording performances
  performances <- c()
  
  # One iteration per fold
  for (fold in 1:k){
    
    # Create training set for this iteration
    # Subset all the datapoints where .folds does not match the current fold
    training_set <- data[data$.folds != fold,]
    
    # Create test set for this iteration
    # Subset all the datapoints where .folds matches the current fold
    testing_set <- data[data$.folds == fold,]
    
    ## Train model
    model <- feols(fml, training_set)
    
    ## Test model
    
    # Predict the dependent variable in the testing_set with the trained model
    predicted <- predict(model, testing_set)
    
    # Get the Root Mean Square Error between the predicted and the observed
    RMSE <- sqrt(mean((predicted - testing_set$log_yield)^2))
    
    # Add the RMSE to the performance list
    performances[fold] <- RMSE
    
    
  }
  
  # Return the mean of the recorded RMSEs
  return(c('RMSE' = mean(performances)))
  
}

########################################################
### Assess convergence function via Rhat and Neff
########################################################
assess_convergence <- function(jags_model, n_params){
  # Rhat
  rhat_all <- jags_model$BUGSoutput$summary[,8]
  rhat_all_summary <- c(min(rhat_all),
                        quantile(rhat_all,c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)),
                        max(rhat_all))
  print('Rhat all:')
  print(rhat_all_summary)
  
  # Neff
  neff_all <- effectiveSize(as.mcmc(jags_model))
  neff_all_summary <- c(min(neff_all),
                        quantile(neff_all,c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)),
                        max(neff_all))
  print('Neff all:')
  print(neff_all_summary)
  
  if (n_params != -1) {
    # Rhat
    rhat_params <- jags_model$BUGSoutput$summary[,8][1:n_params]
    rhat_params_summary <- c(min(rhat_params),
                             quantile(rhat_params,c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)),
                             max(rhat_params))
    print('Rhat parameters only:')
    print(rhat_params_summary)
    
    # Neff
    neff_params <- effectiveSize(as.mcmc(jags_model))
    neff_params_summary <- c(min(neff_params),
                             quantile(neff_params,c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)),
                             max(neff_params))
    print('Neff parameters only:')
    print(neff_params_summary)
  }
}

how_many_samples <- function(n_iter, n_burnin, n_thin, n_chains){
  return ((n_iter-n_burnin) / n_thin * n_chains)
}

#############################################
# Function to get yield projections
#############################################

yield_project <- function(crop, yield.model, county, county.id, nthin, nsim, ncores, models){
  # Read JAGS model
  yield_model <- readRDS(paste('../data/simulated/yield_models/', crop, '_', yield.model, '.rds', sep=''))
  beta1 <- yield_model$BUGSoutput$sims.list$beta1
  beta2 <- yield_model$BUGSoutput$sims.list$beta2
  beta3 <- yield_model$BUGSoutput$sims.list$beta3
  beta4 <- yield_model$BUGSoutput$sims.list$beta4
  beta5 <- yield_model$BUGSoutput$sims.list$beta5
  alpha <- yield_model$BUGSoutput$sims.list$alpha
  tau.y <- yield_model$BUGSoutput$sims.list$tau.y
  rm(yield_model)
  
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
  
  # Thinning of parameters
  npar <- dim(beta1)[1]
  l <-seq(1, npar, nthin)
  
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
      
      # Get list of means, sampling from every nthin-th parameter set
      Ymean <- matrix(c(beta1[l,county.id], beta2[l,county.id], beta3[l,county.id],
                        beta4[l,county.id], beta5[l,county.id]), nrow=npar/nthin, ncol=5)
      Ymean <- Ymean %*% t(clim) + alpha[l,county.id]
      
      # Simulate future yield for given year and county, nsim times for each parameter set
      log_yield_sim <-rnorm(nsim * npar/nthin, Ymean, 1/tau.y[l,county.id])
      
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
                '_RCP45_30-50_', yield.model,'_TREND.csv', sep='')
  fwrite(Ypred.df, path)
}