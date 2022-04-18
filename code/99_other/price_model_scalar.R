setwd(getwd())

library(ggplot2)
library(tidyr)
library(dplyr)
library(rjags)
library(R2jags)
library(runjags)
library(parallel)
library(coda)
library(MCMCpack)
library(reshape2)
library(mnormt)

###################################
##### Maize
###################################

# Read data
maize <- read.csv('../data/processed/prices/maize_ppi.csv')
maize$log_maize  <- log(maize$value_ppi)
maize <- maize %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
maize <- maize[order(maize$date),]
maize_lm <- lm(maize$log_maize ~ maize$date)

# Model
ar1_model <- function () {
  z[1] ~ dnorm(0, tausq)
  zrep[1] ~ dnorm(0, tausq)
  zfuture[1] = z[length(z)] 
  
  for(t in 2:length(z)){
    z[t] ~ dnorm(rho * z[t-1], tausq / (1 - rho^2))
    
    zrep[t] ~ dnorm(rho * zrep[t-1], tausq / (1 - rho^2))
  }
  for(t in 2:npred){
    zfuture[t] ~ dnorm(rho * zfuture[t-1], tausq / (1 - rho^2))
  }
  
  zrep.min <- min(zrep)
  zrep.sd <- sd(zrep)
  zrep.max <- max(zrep)
  
  rho ~ dunif(-1, 1)
  tausq ~ dgamma(0.001, 0.001) # check this!!
  
  sigmasq <- 1 / tausq
}

# JAGS model
npred <- 20
jags.data <- list('z' = residuals(maize_lm),
                  'npred' = npred)

jags.params <- c('rho','sigmasq', 'zrep.min', 'zrep.sd', 'zrep.max', 'zfuture')

jags.inits <- NULL

#=============#
# using jags  #
#=============#

# Fit the model using run.jags -- parallel #
bayes.mod.fit  <- jags.parallel(model.file = ar1_model, 
                                data = jags.data, 
                                parameters.to.save = jags.params,
                                inits = jags.inits,
                                n.iter=10000,
                                n.chains = 3)


all_convergence <- bayes.mod.fit$BUGSoutput$summary[,8]
rhat_summary <- c(min(all_convergence),quantile(all_convergence,c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)),max(all_convergence),mean(all_convergence))
rhat_summary

# Output
bayes.mod.fit
attach.jags(bayes.mod.fit)

# Check zrep
ggplot() + 
  geom_density(data=as.data.frame(zrep.sd), aes(V1)) + 
  geom_vline(aes(xintercept=sd(resid(maize_lm))), color='blue')

# Plot future
zfuture.df <- as.data.frame(t(zfuture))
maize_zfut <- zfuture.df
zfuture.df$date <- seq(as.Date("2021-12-31"), as.Date("2040-12-31"), by="years")
zfuture.df <- melt(zfuture.df, id.vars = 'date')

maize_resid <- as.data.frame(residuals(maize_lm))
colnames(maize_resid) <- c('maize_resid')
maize_resid$date <- maize$date

zfuture.df.q05 <- zfuture.df %>% group_by(date) %>% group_modify(~ {
  quantile(.x$value, probs = c(0.05, 0.95)) %>%
    tibble::enframe(name = "prob", value = "value")
})

zfuture.df.q01 <- zfuture.df %>% group_by(date) %>% group_modify(~ {
  quantile(.x$value, probs = c(0.01, 0.99)) %>%
    tibble::enframe(name = "prob", value = "value")
})

zfuture.df.q50 <- zfuture.df %>% group_by(date) %>% group_modify(~ {
  quantile(.x$value, probs = c(0.5)) %>%
    tibble::enframe(name = "prob", value = "value")
})

ggplot() + 
  geom_line(data=zfuture.df, aes(x=date, y=value, group=variable), alpha=0.5, color='darkgray') +
  geom_line(data=maize_resid, aes(x=date, y=maize_resid), color='red') +
  geom_line(data=zfuture.df.q05, aes(x=date, y=value, group=prob), alpha=1, color='blue') +
  geom_line(data=zfuture.df.q01, aes(x=date, y=value, group=prob), alpha=1, color='darkblue') +
  geom_line(data=zfuture.df.q50, aes(x=date, y=value, group=prob), alpha=1, color='black') +
  ylim(-1,1) + ylab("maize")

###################################
##### Soy
###################################

# Read data
soy <- read.csv('../data/processed/prices/soy_ppi.csv')
soy$log_soy  <- log(soy$value_ppi)
soy <- soy %>% mutate(date = as.Date(date, format = "%Y-%m-%d"))
soy <- soy[order(soy$date),]
soy_lm <- lm(soy$log_soy ~ soy$date)

# Model
ar1_model <- function () {
  z[1] ~ dnorm(0, tausq)
  zrep[1] ~ dnorm(0, tausq)
  zfuture[1] = z[length(z)] 
  
  for(t in 2:length(z)){
    z[t] ~ dnorm(rho * z[t-1], tausq / (1 - rho^2))
    
    zrep[t] ~ dnorm(rho * zrep[t-1], tausq / (1 - rho^2))
  }
  for(t in 2:npred){
    zfuture[t] ~ dnorm(rho * zfuture[t-1], tausq / (1 - rho^2))
  }
  
  zrep.min <- min(zrep)
  zrep.sd <- sd(zrep)
  zrep.max <- max(zrep)
  
  rho ~ dunif(-1, 1)
  tausq ~ dgamma(0.001, 0.001) # check this!!
  
  sigmasq <- 1 / tausq
}

# JAGS model
npred <- 20
jags.data <- list('z' = residuals(soy_lm),
                  'npred' = npred)

jags.params <- c('rho','sigmasq', 'zrep.min', 'zrep.sd', 'zrep.max', 'zfuture')

jags.inits <- NULL

#=============#
# using jags  #
#=============#

# Fit the model using run.jags -- parallel #
bayes.mod.fit  <- jags.parallel(model.file = ar1_model, 
                                data = jags.data, 
                                parameters.to.save = jags.params,
                                inits = jags.inits,
                                n.iter=10000,
                                n.chains = 3)


all_convergence <- bayes.mod.fit$BUGSoutput$summary[,8]
rhat_summary <- c(min(all_convergence),quantile(all_convergence,c(0.01,0.1,0.25,0.5,0.75,0.9,0.99)),max(all_convergence),mean(all_convergence))
rhat_summary

# Output
bayes.mod.fit
attach.jags(bayes.mod.fit)

# Check zrep
ggplot() + 
  geom_density(data=as.data.frame(zrep.sd), aes(V1)) + 
  geom_vline(aes(xintercept=sd(resid(soy_lm))), color='blue')

# Plot future
zfuture.df <- as.data.frame(t(zfuture))
soy_zfut <- zfuture.df
zfuture.df$date <- seq(as.Date("2021-12-31"), as.Date("2040-12-31"), by="years")
zfuture.df <- melt(zfuture.df, id.vars = 'date')

soy_resid <- as.data.frame(residuals(soy_lm))
colnames(soy_resid) <- c('soy_resid')
soy_resid$date <- soy$date

zfuture.df.q05 <- zfuture.df %>% group_by(date) %>% group_modify(~ {
  quantile(.x$value, probs = c(0.05, 0.95)) %>%
    tibble::enframe(name = "prob", value = "value")
})

zfuture.df.q01 <- zfuture.df %>% group_by(date) %>% group_modify(~ {
  quantile(.x$value, probs = c(0.01, 0.99)) %>%
    tibble::enframe(name = "prob", value = "value")
})

zfuture.df.q50 <- zfuture.df %>% group_by(date) %>% group_modify(~ {
  quantile(.x$value, probs = c(0.5)) %>%
    tibble::enframe(name = "prob", value = "value")
})

ggplot() + 
  geom_line(data=zfuture.df, aes(x=date, y=value, group=variable), alpha=0.5, color='darkgray') +
  geom_line(data=soy_resid, aes(x=date, y=soy_resid), color='red') +
  geom_line(data=zfuture.df.q05, aes(x=date, y=value, group=prob), alpha=1, color='blue') +
  geom_line(data=zfuture.df.q01, aes(x=date, y=value, group=prob), alpha=1, color='darkblue') +
  geom_line(data=zfuture.df.q50, aes(x=date, y=value, group=prob), alpha=1, color='black') +
  ylim(-1,1) + ylab("soy")

# Check correlation
zfut_corr <- mapply(cor, maize_zfut, soy_zfut)

ggplot() + 
  geom_density(data=as.data.frame(zfut_corr), aes(zfut_corr))