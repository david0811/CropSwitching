setwd(getwd())

library(dplyr)
library(fixest)
library(ggplot2)
library(groupdata2)
source('./00_useful_R_funcs.R') 

# Load yield data (used for all models!)
usda <- read.csv('../data/historical/yields/maize_all.csv')
usda <- subset(usda, select=c('fips','year','yield', 'state'))
usda <- filter(usda, year <= 2013) # Match GCM period
usda <- filter(usda, year >= 1950) # Match GCM period
usda <- merge(usda,  usda %>% count(fips)) # Filter to counties with >=50% coverage
usda <- filter(usda, n >= 30) # Filter to counties with >=50% coverage
usda$log_yield <- log(usda$yield)

########################### 
## Haqiqi model
########################### 
# Climate data
loca.obs <- read.csv('../data/historical/predictors/Livneh/Livneh_maize_H21_predictors.csv')

################### 1950-, county trends
# Merge
df <- merge(loca.obs, usda)
df$year2 <- df$year **2

# Models (Haqiqi)
fx.mod1 <- feols(log_yield ~ GDD + EDD + prcp + prcp2 | 
                   fips[year, year2] + fips, df, cluster = "state")
fx.mod2 <- feols(log_yield ~ GDD + EDD + SM_mean + SM_mean2 |
                   fips[year, year2] + fips, df, cluster = 'state')
fx.mod3 <- feols(log_yield ~ GDD + EDD + SM_week_max + SM_week_min | 
                   fips[year, year2] + fips, df, cluster = "state")
fx.mod4 <- feols(log_yield ~ GDD + EDD_SM_75_below + EDD_SM_25_75_below + 
                   EDD_SM_0_25_norm + EDD_SM_25_75_above + EDD_SM_75_above |
                   fips[year, year2] + fips, df, cluster = "state")
fx.mod5 <- feols(log_yield ~ GDD + prcp + prcp2 + EDD_SM_75_below + EDD_SM_25_75_below + 
                   EDD_SM_0_25_norm + EDD_SM_25_75_above + EDD_SM_75_above | 
                   fips[year, year2] + fips, df, cluster = "state")
fx.mod6 <- feols(log_yield ~ GDD + SM_mean + SM_mean2 + EDD_SM_75_below + EDD_SM_25_75_below + 
                   EDD_SM_0_25_norm + EDD_SM_25_75_above + EDD_SM_75_above | 
                   fips[year, year2] + fips, df, cluster = "state")

# Results table
resH1 <- etable(fx.mod1, fx.mod2, fx.mod3, fx.mod4, fx.mod5, fx.mod6, fitstat =~ . + aic + bic + rmse)
resH1

# Cross validation
nfolds = 10
train_set <- fold(df, k = nfolds, cat_col = 'year')
train_set <- train_set %>% arrange(.folds)

crossvalidate(train_set, k = nfolds, fml = formula(fx.mod1))
crossvalidate(train_set, k = nfolds, fml = formula(fx.mod2))
crossvalidate(train_set, k = nfolds, fml = formula(fx.mod3))
crossvalidate(train_set, k = nfolds, fml = formula(fx.mod4))
crossvalidate(train_set, k = nfolds, fml = formula(fx.mod5))
crossvalidate(train_set, k = nfolds, fml = formula(fx.mod6))


################### 1981-. state trends
usda <- filter(usda, year >= 1980)

# Merge
df <- merge(loca.obs, usda)
df$year2 <- df$year **2

# Models (Haqiqi)
fx.mod1 <- feols(log_yield ~ GDD + EDD + prcp + prcp2 | 
                   state[year, year2] + fips, df, cluster = "state")
fx.mod2 <- feols(log_yield ~ GDD + EDD + SM_mean + SM_mean2 |
                   state[year, year2] + fips, df, cluster = 'state')
fx.mod3 <- feols(log_yield ~ GDD + EDD + SM_week_max + SM_week_min | 
                   state[year, year2] + fips, df, cluster = "state")
fx.mod4 <- feols(log_yield ~ GDD + EDD_SM_75_below + EDD_SM_25_75_below + 
                   EDD_SM_0_25_norm + EDD_SM_25_75_above + EDD_SM_75_above |
                   state[year, year2] + fips, df, cluster = "state")
fx.mod5 <- feols(log_yield ~ GDD + prcp + prcp2 + EDD_SM_75_below + EDD_SM_25_75_below + 
                   EDD_SM_0_25_norm + EDD_SM_25_75_above + EDD_SM_75_above | 
                   state[year, year2] + fips, df, cluster = "state")
fx.mod6 <- feols(log_yield ~ GDD + SM_mean + SM_mean2 + EDD_SM_75_below + EDD_SM_25_75_below + 
                   EDD_SM_0_25_norm + EDD_SM_25_75_above + EDD_SM_75_above | 
                   state[year, year2] + fips, df, cluster = "state")

# Results table
resH2 <- etable(fx.mod1, fx.mod2, fx.mod3, fx.mod4, fx.mod5, fx.mod6, fitstat =~ . + aic + bic + rmse)
resH2

########################### 
## Intraseasonal model
########################### 
# Load yield data
usda <- read.csv('../data/historical/yields/maize_all.csv')
usda <- subset(usda, select=c('fips','year','yield', 'state'))
usda <- filter(usda, year <= 2013) # Match GCM period
usda <- filter(usda, year >= 1950) # Match GCM period
usda <- merge(usda,  usda %>% count(fips)) # Filter to counties with >=50% coverage
usda <- filter(usda, n >= 30) # Filter to counties with >=50% coverage
usda$log_yield <- log(usda$yield)

# Climate data
loca.obs <- read.csv('../data/historical/predictors/Livneh/Livneh_maize_intra_predictors.csv')

################### 1950-, county trends
# Merge
df <- merge(loca.obs, usda)
df$year2 <- df$year **2

# Models
fx.mod0 <- feols(log_yield ~ GDD_4 + EDD_4 + 
                   GDD_5 + EDD_5 + 
                   GDD_6 + EDD_6 + 
                   GDD_7 + EDD_7 + 
                   GDD_8 + EDD_8 + 
                   GDD_9 + EDD_9 | 
                   fips[year] + fips, df, cluster = "state")

fx.mod1 <- feols(log_yield ~ GDD_4 + EDD_4 + prcp_4 + prcp2_4 + 
                   GDD_5 + EDD_5 + prcp_5 + prcp2_5 + 
                   GDD_6 + EDD_6 + prcp_6 + prcp2_6 + 
                   GDD_7 + EDD_7 + prcp_7 + prcp2_7 + 
                   GDD_8 + EDD_8 + prcp_8 + prcp2_8 + 
                   GDD_9 + EDD_9 + prcp_9 + prcp2_9 | 
                   fips[year] + fips, df, cluster = "state")

fx.mod2 <- feols(log_yield ~ GDD_4 + EDD_4 + SM_mean_4 + SM_mean2_4 + 
                   GDD_5 + EDD_5 + SM_mean_5 + SM_mean2_5 + 
                   GDD_6 + EDD_6 + SM_mean_6 + SM_mean2_6 + 
                   GDD_7 + EDD_7 + SM_mean_7 + SM_mean2_7 + 
                   GDD_8 + EDD_8 + SM_mean_8 + SM_mean2_8 + 
                   GDD_9 + EDD_9 + SM_mean_9 + SM_mean2_9 | 
                   fips[year] + fips, df, cluster = "state")

fx.mod3 <- feols(log_yield ~ GDD_4 + EDD_4 + SM_week_min_4 + SM_week_max_4 +
                   GDD_5 + EDD_5 + SM_week_min_5 + SM_week_max_5 +
                   GDD_6 + EDD_6 + SM_week_min_6 + SM_week_max_6 +
                   GDD_7 + EDD_7 + SM_week_min_7 + SM_week_max_7 +
                   GDD_8 + EDD_8 + SM_week_min_8 + SM_week_max_8 +
                   GDD_9 + EDD_9 + SM_week_min_9 + SM_week_max_9 | 
                   fips[year] + fips, df, cluster = "state")

# Results table
resI1 <- etable(fx.mod0, fx.mod1, fx.mod2, fx.mod3, fitstat =~ . + aic + bic + rmse)
resI1

########################### 1980-, county trends
usda <- filter(usda, year >= 1980)

# Merge
df <- merge(loca.obs, usda)
df$year2 <- df$year **2

# Models
fx.mod0 <- feols(log_yield ~ GDD_4 + EDD_4 + 
                   GDD_5 + EDD_5 + 
                   GDD_6 + EDD_6 + 
                   GDD_7 + EDD_7 + 
                   GDD_8 + EDD_8 + 
                   GDD_9 + EDD_9 | 
                   state[year, year2] + fips, df, cluster = "state")

fx.mod1 <- feols(log_yield ~ GDD_4 + EDD_4 + prcp_4 + prcp2_4 + 
                   GDD_5 + EDD_5 + prcp_5 + prcp2_5 + 
                   GDD_6 + EDD_6 + prcp_6 + prcp2_6 + 
                   GDD_7 + EDD_7 + prcp_7 + prcp2_7 + 
                   GDD_8 + EDD_8 + prcp_8 + prcp2_8 + 
                   GDD_9 + EDD_9 + prcp_9 + prcp2_9 | 
                   state[year, year2] + fips, df, cluster = "state")

fx.mod2 <- feols(log_yield ~ GDD_4 + EDD_4 + SM_mean_4 + SM_mean2_4 + 
                   GDD_5 + EDD_5 + SM_mean_5 + SM_mean2_5 + 
                   GDD_6 + EDD_6 + SM_mean_6 + SM_mean2_6 + 
                   GDD_7 + EDD_7 + SM_mean_7 + SM_mean2_7 + 
                   GDD_8 + EDD_8 + SM_mean_8 + SM_mean2_8 + 
                   GDD_9 + EDD_9 + SM_mean_9 + SM_mean2_9 | 
                   state[year, year2] + fips, df, cluster = "state")

fx.mod3 <- feols(log_yield ~ GDD_4 + EDD_4 + SM_week_min_4 + SM_week_max_4 +
                   GDD_5 + EDD_5 + SM_week_min_5 + SM_week_max_5 +
                   GDD_6 + EDD_6 + SM_week_min_6 + SM_week_max_6 +
                   GDD_7 + EDD_7 + SM_week_min_7 + SM_week_max_7 +
                   GDD_8 + EDD_8 + SM_week_min_8 + SM_week_max_8 +
                   GDD_9 + EDD_9 + SM_week_min_9 + SM_week_max_9 | 
                   state[year, year2] + fips, df, cluster = "state")

# Results table
resI2 <- etable(fx.mod0, fx.mod1, fx.mod2, fx.mod3, fitstat =~ . + aic + bic + rmse)
resI2

# Cross Validation
nfolds <- 10
train_set <- fold(df, k = nfolds, cat_col = 'year')
train_set <- train_set %>% arrange(.folds)

crossvalidate(train_set, k = nfolds, fml = formula(fx.mod1))
crossvalidate(train_set, k = nfolds, fml = formula(fx.mod2))
crossvalidate(train_set, k = nfolds, fml = formula(fx.mod3))

