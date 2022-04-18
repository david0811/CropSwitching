setwd(getwd())

library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
#library(reshape2)
library(rjags)
library(R2jags)
library(parallel)
library(coda)
library(mnormt)
source('./00_useful_R_funcs.R')

###############################################
#### Read and tidy data 
###############################################
barley <- read.csv('../data/historical/prices/barley_ppi.csv')
barley$log_barley <- log(barley$value_ppi)

maize <- read.csv('../data/historical/prices/maize_ppi.csv')
maize$log_maize  <- log(maize$value_ppi)

soy <- read.csv('../data/historical/prices/soy_ppi.csv')
soy$log_soy  <- log(soy$value_ppi)

cotton <- read.csv('../data/historical/prices/cotton_ppi.csv')
cotton$log_cotton  <- log(cotton$value_ppi)

rice <- read.csv('../data/historical/prices/rice_ppi.csv')
rice$log_rice  <- log(rice$value_ppi)

wheat <- read.csv('../data/historical/prices/wheat_ppi.csv')
wheat$log_wheat  <- log(wheat$value_ppi)

sorghum <- read.csv('../data/historical/prices/sorghum_ppi.csv')
sorghum$log_sorghum  <- log(sorghum$value_ppi)

prices <- subset(soy, select=c('date', 'log_soy')) %>%
                merge(subset(maize, select=c('date', 'log_maize'))) %>%
                merge(subset(cotton, select=c('date', 'log_cotton'))) %>%
                merge(subset(rice, select=c('date', 'log_rice'))) %>%
                merge(subset(sorghum, select=c('date', 'log_sorghum'))) %>%
                merge(subset(wheat, select=c('date', 'log_wheat'))) %>%
                merge(subset(barley, select=c('date', 'log_barley'))
)

prices <- prices %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
prices$year <- format(prices$date, format="%Y")
prices <- subset(prices, select=-c(date))
prices$year <- as.numeric(as.character(prices$year))

# check plot
df <- as.data.frame(prices)
df <- pivot_longer(df, cols=c('log_maize', 'log_soy', 'log_cotton', 'log_barley',
                        'log_sorghum', 'log_rice', 'log_wheat'), names_to = 'crops')
ggplot(df, aes(year,value)) + geom_line(aes(colour=crops)) +
  ylab('log(Price) [2018$]') + xlab('Year')

# De-trend log prices
prices <- prices[order(prices$year),]
barley_lm <- lm(prices$log_barley ~ prices$year)
maize_lm <- lm(prices$log_maize ~ prices$year)
soy_lm <- lm(prices$log_soy ~ prices$year)
cotton_lm <- lm(prices$log_cotton ~ prices$year)
wheat_lm <- lm(prices$log_wheat ~ prices$year)
sorghum_lm <- lm(prices$log_sorghum ~ prices$year)
rice_lm <- lm(prices$log_rice~ prices$year)

prices_fit <- list(prices$year,
                   residuals(soy_lm),
                   residuals(maize_lm),
                   residuals(cotton_lm),
                   residuals(rice_lm),
                   residuals(sorghum_lm),
                   residuals(wheat_lm),
                   residuals(barley_lm))

crops_date <- c('year', 'soy', 'maize', 'cotton', 'rice', 'sorghum', 'wheat', 'barley')
crops <- c('soy', 'maize', 'cotton', 'rice', 'sorghum', 'wheat', 'barley')
names(prices_fit) <- crops_date
prices_fit <- as.data.frame(prices_fit)

# Plot
df_fit <- pivot_longer(prices_fit, cols=all_of(crops), names_to = 'crops')
ggplot(df_fit, aes(year,value)) + geom_line(aes(colour = crops)) +
 ylab('log(Price) [2018$]') + xlab('Year')

# For model
prices_fit <- subset(prices_fit, select=-c(year))

# ###############################################
# #### Bayesian VAR model: manual
# ###############################################
# # Model
# bvar1_model <- function () {
#   # Initial conditions
#   z[1,1:p] ~ dmnorm(rep(0,p), Sigma.inv)
#   zrep[1,1:p] ~ dmnorm(rep(0,p), Sigma.inv)
#   zfut[1,1:p] <- z[nyrs,]
#   
#   # Loop through years
#   for (t in 2:nyrs){
#     # Define Beta vector of coefficients
#     for (pp in 0:(p-1)){
#       Beta[t, pp+1] <- inprod(rho[(pp*p + 1) : (pp*p + p)], z[t-1,1:p])
#       Beta.rep[t, pp+1] <- inprod(rho[(pp*p + 1) : (pp*p + p)], zrep[t-1,1:p])
#     }
#     #  Multi-variate normal likelihood
#     z[t,1:p] ~ dmnorm(Beta[t,], Sigma.inv)
#     zrep[t,1:p] ~ dmnorm(Beta.rep[t,], Sigma.inv)
#   }
#   # Future prediction
#   for (t in 2:(npred+1)){
#     for (pp in 0:(p-1)){
#       Beta.fut[t, pp+1] <- inprod(rho[(pp*p + 1) : (pp*p + p)], zfut[t-1,1:p])
#     }
#     zfut[t,1:p] ~ dmnorm(Beta.fut[t,], Sigma.inv)
#   }
#   
#   # Priors
#   rho ~ dmnorm(mu0, Prec)
#   Sigma.inv ~ dscaled.wishart(rep(0.5, p), df)
#   
#   # Record covariance matrix
#   Sigma <- inverse(Sigma.inv)
# }
# 
# # Dimensions
# p <- ncol(prices_fit)
# nyrs <- nrow(prices_fit)
# npred <- 40
# 
# jags.data <- list('z' = as.matrix(prices_fit),
#                   'Prec' = diag(rep(0.001, p**2)),
#                   'mu0' = rep(0, p**2),
#                   'df' = 2,
#                   'p' = p,
#                   'npred' = npred,
#                   'nyrs'= nyrs
#                   )
# 
# jags.params <- c('rho', 'zrep', 'zfut', 'Sigma')
# 
# jags.inits <- NULL
# 
# how_many_samples(n_iter=100000, n_burnin=10000, n_thin=5, n_chains=3)
# 
# # Fit the model using run.jags -- parallel. Takes around 10 mins#
# start.time <- Sys.time()
# bayes.mod.fit  <- jags.parallel(model.file = bvar1_model, 
#                                 data = jags.data, 
#                                 parameters.to.save = jags.params,
#                                 inits = jags.inits,
#                                 n.iter = 100000L,
#                                 n.burnin = 10000L,
#                                 n.thin = 5L,
#                                 n.chains = 3)
# print(Sys.time() - start.time)
# 
# ## Assess convergence
# assess_convergence(bayes.mod.fit, 2*p**2 + 1)
# 
# # Save model
# saveRDS(bayes.mod.fit, '../data/simulated/prices/BVAR_allcrops.rds')

# Load model
bayes.mod.fit <- readRDS('../data/simulated/prices/BVAR_allcrops.rds')

# Attach
attach.jags(bayes.mod.fit)

# Dimensions
p <- ncol(prices_fit)
nyrs <- nrow(prices_fit)
npred <- 40

## Store predictions
for (ncrop in 1:length(crops)) {
  zfuture.df <- as.data.frame(t(zfut[,,ncrop]))
  colnames(zfuture.df) <-  seq(1, dim(zfuture.df)[2])
  cols_to_pivot <- colnames(zfuture.df)
  zfuture.df$year <- seq(2020, 2020 + npred)
  zfuture.df <- pivot_longer(zfuture.df, cols=all_of(cols_to_pivot), names_to='param')
  colnames(zfuture.df) <- c('year', 'param', 'log_price')

  # Assuming no trend: add fluctuations to avg of last 20 years
  last20_avg <- prices %>% filter(year >= 2000) %>% 
    select(paste('log_', crops[ncrop], sep='')) %>% lapply(mean)

  zfuture.df$log_price_NOTREND <- zfuture.df$log_price + last20_avg[[1]]

  # Assuming continuation of historical trend in log-prices
  lm <- get(paste(crops[ncrop], '_lm', sep=''))
  zfuture.df$log_price_TREND <- lm$coefficients[1] + lm$coefficients[2]*zfuture.df$year
  zfuture.df$log_price_TREND <- zfuture.df$log_price_TREND + zfuture.df$log_price

  # Subset to climate data time period
  zfuture.df <- filter(zfuture.df, year >= 2030 & year <= 2050)

  # Save
  fwrite(zfuture.df, 
         paste('../data/future/prices/', crops[ncrop], '_30-50_BVAR.csv', sep=''))
}

############### Analysis
# Covariance matrix
vcov_mean <- apply(Sigma, c(2,3), mean)
vcov_mean

# Check correlation
zrep_corr <- mapply(cor, as.data.frame(t(zrep[,,1])),
                    as.data.frame(t(zrep[,,2])))

zfut_corr <- mapply(cor, as.data.frame(t(zfut[,,1])),
                    as.data.frame(t(zfut[,,2])))

ggplot() + 
  geom_density(data=as.data.frame(zrep_corr), aes(zrep_corr), color='gray') +
  geom_density(data=as.data.frame(zfut_corr), aes(zfut_corr), color='black') +
  geom_vline(aes(xintercept=cor(prices_fit[,1], prices_fit[,2])), color='blue')

# Plot future
plot_prediction <- function(ncrop){
  # Historical data
  crop_hist <- as.data.frame(prices[ncrop])
  colnames(crop_hist) <- c('value')
  crop_hist$date <- prices$year
  crop_hist$value <- exp(crop_hist$value)
  
  # Get predictions
  zfuture.df <- as.data.frame(t(zfut[,,ncrop]))
  zfuture.df$date <- seq(2020, 2020 + npred)
  zfuture.df <- pivot_longer(zfuture.df,
                           cols=colnames(zfuture.df)[1:(length(colnames(zfuture.df))-1)])
  zfuture.df$value <- exp(log(crop_hist$value[nyrs]) + zfuture.df$value - zfuture.df$value[1])

  # Calculate quantiles for credible intervals
  zfuture.df.q05 <- zfuture.df %>% group_by(date) %>% group_modify(~ {
    quantile(.x$value, probs = c(0.05, 0.95)) %>%
    tibble::enframe(name = "prob", value = "value")
  })

  zfuture.df.q05 <- pivot_wider(zfuture.df.q05, names_from = prob, values_from = value)

  zfuture.df.q01 <- zfuture.df %>% group_by(date) %>% group_modify(~ {
    quantile(.x$value, probs = c(0.01, 0.99)) %>%
    tibble::enframe(name = "prob", value = "value")
  })
  zfuture.df.q01 <- pivot_wider(zfuture.df.q01, names_from = prob, values_from = value)

  zfuture.df.q50 <- zfuture.df %>% group_by(date) %>% group_modify(~ {
    quantile(.x$value, probs = c(0.5)) %>%
    tibble::enframe(name = "prob", value = "value")
  })
  
  # Get one random sample
  # zfuture.df <- as.data.frame(zfut[sample(1:dim(rho)[1], 1),,ncrop])
  zfuture.df <- as.data.frame(zfut[10000,,ncrop])
  zfuture.df$date <- seq(2020, 2020 + npred)
  zfuture.df <- pivot_longer(zfuture.df,
                           cols=colnames(zfuture.df)[1:(length(colnames(zfuture.df))-1)])
  zfuture.df$value <- exp(log(crop_hist$value[nyrs]) + zfuture.df$value - zfuture.df$value[1])

  ggplot() + 
    geom_line(data=crop_hist, aes(x=date, y=value), lwd=0.8, color='black') +
    geom_ribbon(data=zfuture.df.q01, aes(x=date, ymin=`1%`, ymax=`99%`), alpha=0.7, fill='#8ACFFF') +
    geom_ribbon(data=zfuture.df.q05, aes(x=date, ymin=`5%`, ymax=`95%`), alpha=0.9, fill='#8ACFFF') +
    geom_line(data=zfuture.df.q50, aes(x=date, y=value), lwd=0.8, color='white') +
    geom_line(data=zfuture.df, aes(x=date, y=value, group=name), lwd=0.8, color='#fc8d62') +
    xlab('Year') + ylab(expression('Price per bushel [$'[2018]*']')) +
    scale_x_continuous(breaks = seq(1920, 2060, by = 10)) +
    ggtitle(paste(toupper(substr(crops[ncrop], 1, 1)), substr(crops[ncrop], 2, nchar(crops[ncrop])), sep=""))
}

plot_prediction(1)
plot_prediction(2)
plot_prediction(3)
plot_prediction(4)
plot_prediction(5)
plot_prediction(6)
plot_prediction(7)

###############################################
#### Bayesian VAR model: BVAR package
####
library(BVAR)
library(BVARverse)
###############################################
# Hyperprior
mn <- bv_minnesota(lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 0.0001, max = 5),
                   alpha = bv_alpha(mode = 2, sd = 0.25, min = 1, max = 3),
                   psi = bv_psi(),
                   var = 1e7)

priors <- bv_priors(hyper = "full", mn = mn)

# MH settings
mh <- bv_metropolis(scale_hess = 0.05,
                    adjust_acc = TRUE, acc_lower = 0.25, acc_upper = 0.45)

# Run 3 chains and check Gelman-Rubin convergence
n_draw <- 250000L
n_burn <- 100000L
cl <- makeCluster(3L)
runs <- par_bvar(cl, n_runs = 3L,
                 data = prices_fit, lags = 1, n_draw = n_draw, n_burn = n_burn, n_thin = 10L,
                 priors = priors, mh = mh)
stopCluster(cl)

runs_mcmc <- as.mcmc(runs, vars = crops)
gelman.diag(runs_mcmc, autoburnin=FALSE, confidence = 0.99)

# Run single model
run <- bvar(prices_fit, lags = 1, n_draw = n_draw, n_burn = n_burn, n_thin = 10L,
            priors = priors, mh = mh, verbose = TRUE)
run_mcmc <- as.mcmc(run, vars=crops)
effectiveSize(run_mcmc)

summary(run)
plot(run)

run.res <- augment(run, conf_bands = 0.025)
run.res

# Plot all
bv_ggplot(run, type='density', vars=crops) #, chains=runs)

# Predictions
run.predict <- predict(run, horizon = 20, n_thin=1L, conf_bands = c(0.01, 0.05))

bv_ggplot(run.predict, vars = c('maize'), t_back=100L, col = "#4455f2")

# Individual predictions: these look way too narrow!
ncrop <- 2
crop.predict <- data.frame(t(run.predict[['fcast']][,,ncrop]))
crop.predict$date <- seq(as.Date("2021-12-31"), as.Date("2040-12-31"), by="years")
crop.predict <- melt(crop.predict, id.vars = 'date')

crop.predict.df.q5 <- crop.predict %>% group_by(date) %>% group_modify(~ {
  quantile(.x$value, probs = c(0.05, 0.95)) %>%
    tibble::enframe(name = "prob", value = "value")
})

crop.predict.df.q01 <- crop.predict %>% group_by(date) %>% group_modify(~ {
  quantile(.x$value, probs = c(0.01, 0.99)) %>%
    tibble::enframe(name = "prob", value = "value")
})

crop.predict.df.q50 <- crop.predict %>% group_by(date) %>% group_modify(~ {
  quantile(.x$value, probs = c(0.5)) %>%
    tibble::enframe(name = "prob", value = "value")
})

prices_plot <- prices_fit
prices_plot$date <- prices$date
ggplot() + 
  geom_line(data=crop.predict, aes(x=date, y=value, group = variable), alpha=0.3, color='darkgray') +
  geom_line(data=prices_plot, aes_string(x="date", y=crops[ncrop]), color='red') +
  geom_line(data=crop.predict.df.q5, aes(x=date, y=value, group=prob), alpha=1, color='blue') +
  geom_line(data=crop.predict.df.q01, aes(x=date, y=value, group=prob), alpha=1, color='darkblue') +
  geom_line(data=crop.predict.df.q50, aes(x=date, y=value, group=prob), alpha=1, color='black') +
  ylim(-1,0.6) + ylab(crops[ncrop])

###############################################
#### Frequentist VAR model
####
library(vars)
###############################################
prices <- prices[order(prices$date),]
prices_fit <- subset(prices, select=-c(date))

# Get order selection 
VARselect(prices_fit, type='none')

# Fit model
var1 = VAR(prices_fit, p=1, type="both")
summary(var1)
roots(var1) # want to be <1
# plot residual ACF
acf(residuals(var1)[,7])
# testing serial correlation
var1.serial <- serial.test(var1, lags.pt=16, type="PT.asymptotic")
var1.serial
plot(var1.serial , names = 'log_maize')
# testing heteroscedasticity
var1.arch <- arch.test(var1, lags.multi = 5, multivariate.only=TRUE)
var1.arch
# test normality
var1.norm <- normality.test(var1, multivariate.only = TRUE)
var1.norm
# look at resids
resid.test <- residuals(var1)[,6]
qqnorm(resid.test)
qqline(resid.test)

# predict
fanchart(predict(var1,n.ahead = 20, ci = 0.95))
