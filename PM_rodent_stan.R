# Models for pine marten population growth rate + pine marten abundance
# Author: Jennifer Angoh
# Date: 12.06.2023

rm(list=ls()) 
setwd("C:/Marten_Vole_Codes")

# Pine marten population growth rate ###########################################
# Load data
Growth_recov_vole_abund <- readRDS("Growth_recov_vole_abund.rds")
hist(Growth_recov_vole_abund$Growth, breaks=100)
length(unique(Growth_recov_vole_abund$groupID))

# MODEL_pine marten growth rate
# Independent variables: vole_abundance; snow; elevation; elev*vole
# Dependent variable: Growth (pine marten)
# Random effect: groupID
# Family = gaussian

### Fit model using rstanarm 
library("rstanarm")
PM_GR_gp <- stan_glmer(Growth ~ 1 + reRoAbund + reSnow + reElev + reElev*reRoAbund + (1|groupID),
                           data = Growth_recov_vole_abund,
                           family = gaussian,
                           prior = normal(0,2), #prior for all slope, here prior pretty flat
                           prior_intercept = normal(0,2),#prior for intercept
                           adapt_delta=0.99,
                           iter = 6000)

saveRDS(PM_GR_gp, file = "PM_GR_gp.rds")

### Model checks 
summary(PM_GR_gp)
## Check convergence of algorism in MCMC diagnostics with n_eff and Rhat:
# n_eff: effective number of drawn, the bigger the better, 
# and should be bigger than sample size,if too small, need to double number of iteration and try again
# if n_eff too small it means thta there may be autocorrelation in the samples
# Rhat: needs to be equal to 1, if it is higher the iterations should be increased, 
# if the chains have converged to the same estimate, then the variance within chain 
# would be equal to variance across chain

## Check convergence, all the chain lines should be mixing well
library(bayesplot) #check help page for different plots that can be ploted to visualise the data
orange_scheme <- c("#ffebcc", "#ffcc80",
                   "#ffad33", "#e68a00",
                   "#995c00", "#663d00")
color_scheme_set(orange_scheme)
plot(PM_GR_gp, plotfun="trace", pars ="beta") #only the post warmup for rstanarm

## Check mean_PPD and this value should be similar to mean of observed data (e.g. mean(Growth))
Growth_recov_vole_abund$Growth |> mean() 
posterior_predict(PM_GR_gp) |> mean() 

## Summary for posterior estimates 80% CRI 
summary_stats <- function(posterior) {
  x <- (posterior)  # log-odds -> probabilities
  t(apply(x, 2, quantile, probs = c(0.1, 0.5, 0.9))) 
}
# as.matrix extracts the posterior draws
estimates <- summary_stats(as.matrix(PM_GR_gp))  

## Posterior predictive check
# give info on the fit of the model 
# generate data from the model after beliefs are updated with observed data 
pp_check (PM_GR_gp) 
# orangge lines are samples from the generate data 
# black line is observed data

### Plot posterior estimates for each parameter
# Plot estimate and associate distribution (and variance)
library(bayesplot)
p <- plot(PM_GR_gp, pars ="beta", prob = 0.8, prob_outer = 0.95) + theme_classic()
p + scale_x_continuous(breaks=c(-1,0.0,1)) 


# Pine marten abundance ########################################################
# Load data
PMtracks_recov_vole_noNA <- readRDS("PMtracks_recov_vole_noNA.rds") 
hist(PMtracks_recov_vole_noNA$PM_index, breaks=100)
length(unique(PMtracks_recov_vole_noNA$Linenr))
sum(PMtracks_recov_vole_noNA$Marten,na.rm=TRUE)

# MODEL_pine marten abundance
# Independent variables: mature spruce; agriculture; snow; vole_abundance; elevation; elevation*vole; 
# Dependent variable: Marten (number of pine marten tracks observed per transect line)
# Random effect: year
# Offsets: days since the last snowfall; track length 
# Family = negative binomial
PM_Abund <- stan_glmer.nb ( Marten ~ 1 + reSpruce + reAgri + reSnow + reRoAbund + reElev + reElev*reRoAbund + (1|Year) + offset(log(length * SnowDays)),
                              data = PMtracks_recov_vole_noNA, 
                              prior = normal(0,2), #prior for all slope, here prior pretty flat
                              prior_intercept = normal(0,2),#prior for intercept
                              adapt_delta=0.99,
                              iter = 6000)

saveRDS(PM_Abund, file = "PM_Abund.rds")

## Model checks 
summary(PM_Abund)
## Check convergence of algorism in MCMC diagnostics with n_eff and Rhat:
# n_eff: effective number of drawn, the bigger the better, 
# and should be bigger than sample size,if too small, need to double number of iteration and try again
# if n_eff too small it means that there may be autocorrelation in the samples
# Rhat: needs to be equal to 1, if it is higher the iterations should be increased, 
# if the chains have converged to the same estimate, then the variance within chain 
# would be equal to variance across chain

## Check convergence, all the chain lines should be mixing well
library(bayesplot) #check help page for different plots that can be ploted to visualise the data
color_scheme_set("blue")
plot(PM_Abund, plotfun="trace", pars ="beta") #only the post warmup for rstanarm

## Check mean_PPD and this value should be similar to mean of observed data (e.g. mean(Marten))
PMtracks_recov_vole_noNA$Marten |> mean()
posterior_predict(PM_Abund) |> mean()

## Summary for posterior estimates 80% CRI
summary_stats <- function(posterior) {
  x <- (posterior)  # log-odds -> probabilities
  t(apply(x, 2, quantile, probs = c(0.1, 0.5, 0.9))) 
}
# as.matrix extracts the posterior draws
estimates <- summary_stats(as.matrix(PM_Abund))  

## Posterior predictive check
# give info on the fit of the model 
# generate data from the model after beliefs are updated with observed data 
pp_check (PM_Abund) 
# blue lines are samples from the generate data
# black line is observed data

### Plot posterior estimates for each parameter
# Plot estimate and associate distribution (and variance)
library(bayesplot)
q <- plot(PM_Abund, pars ="beta", prob = 0.8, prob_outer = 0.95) + theme_classic() 
q + scale_x_continuous(breaks=c(-1,0.0,1))

