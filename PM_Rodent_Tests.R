# Tests for pine marten_microtine rodent study
# Author: Jennifer Angoh
# Date: 29.03.2024

rm(list=ls()) 
setwd("C:/Marten_Vole_Codes")

################################################################################
## Microtine rodent index validation #################################

# Load data
library(readxl)
rodent_peaks <- read_excel("C:/Marten_Vole_codes/rodent_peaks.xlsx")

#Test for correlation ##########################################################
# correlation plot 
correlation <- cor(
  rodent_peaks[,
                    c("Engerdal",
                      "Gutulia")], 
  y= NULL, 
  use="complete.obs",
  method = "pearson"
)

# Correlation chart for temporal level
corr.data.loc <- data.frame(rodent_peaks$Engerdal, 
                            rodent_peaks$Gutulia)

names(corr.data.loc)[1] <- "Engerdal"
names(corr.data.loc)[2] <- "Gutulia"

# Check correlation and significance of interaction between variables
library(Hmisc)
rcorr.m1 <- rcorr(as.matrix(corr.data.loc))
rcorr.m1$r
rcorr.m1$P 


################################################################################
## Zero inflated NB vs. NB: LOOIC comparison ###################################

# Pine marten abundance models in brms 
# Load data
PMtracks_recov_vole_noNA <- readRDS("PMtracks_recov_vole_noNA_290324.rds")

# Negative binomial model ######################################################
library("brms")
PM_Abund_brms_nb_dist <- brms::brm(
  formula = Marten ~ 1 + reSpruce + reAgri + reSnow + reRoAbund + reElev + reElev*reRoAbund + (1|Year) + offset(log(length * SnowDays)),
  data = PMtracks_recov_vole_noNA,
  family = brmsfamily("negbinomial"),
  prior = c(
    prior(normal(0, 2), class = "b"),   # Prior for all slope parameters
    prior(normal(0, 2), class = "Intercept")  # Prior for the intercept
  ),
  control = list(adapt_delta = 0.99),
  iter = 6000
)


saveRDS(PM_Abund_brms_nb, file = "PM_Abund_brms_nb_290324.rds")
PM_Abund_brms_nb <- readRDS("PM_Abund_brms_nb_290324.rds")

plot(PM_Abund_brms_nb)

#plot posterior estimates
library(bayesplot)

# Extract posterior samples from the model
posterior_samples <- as.matrix(PM_Abund_brms_nb)

# Get the parameter names
param_names <- colnames(posterior_samples)

# Filter parameter names starting with "b_"
beta_param_names <- param_names[grep("^b_", param_names)]

# Plot posterior intervals for regression coefficients
p <- mcmc_intervals(
  posterior_samples, 
  pars = beta_param_names, 
  prob = 0.8, # Adjust this value as needed
  prob_outer = FALSE,
)

# Print the plot
print(p)


# Zero inflated negative binomial model ########################################
library("brms")
PM_Abund_brms_zero <- brms::brm(
  formula = Marten ~ 1 + reSpruce + reAgri + reSnow + reRoAbund + reElev + reElev*reRoAbund + (1|Year) + offset(log(length * SnowDays)),
  data = PMtracks_recov_vole_noNA,
  family = brmsfamily("zero_inflated_negbinomial"),
  prior = c(
    prior(normal(0, 2), class = "b"),   # Prior for all slope parameters
    prior(normal(0, 2), class = "Intercept")  # Prior for the intercept
  ),
  control = list(adapt_delta = 0.99),
  iter = 6000
)

saveRDS(PM_Abund_brms_zero, file = "PM_Abund_brms_zero_290324.rds")
PM_Abund_brms_zero <- readRDS("PM_Abund_brms_zero_290324.rds")


plot(PM_Abund_brms_zero)

#plot posterior estimates
library(bayesplot)

# Extract posterior samples from the model
posterior_samples_zero <- as.matrix(PM_Abund_brms_zero)

# Get the parameter names
param_names_zero <- colnames(posterior_samples_zero)

# Filter parameter names starting with "b_"
beta_param_names_zero <- param_names[grep("^b_", param_names_zero)]

# Plot posterior intervals for regression coefficients
q <- mcmc_intervals(
  posterior_samples_zero, 
  pars = beta_param_names_zero, 
  prob = 0.8, # Adjust this value as needed
  prob_outer = FALSE,
)

# Print the plot
print(q)


# LOOIC Comparison #############################################################
# Load the loo package
library(loo)

# Compute the LOOIC for each model
loo_zero <- loo(PM_Abund_brms_zero)
loo_zero
loo_negbinomial <- loo(PM_Abund_brms_nb)
loo_negbinomial

# Compare the models using the LOOIC
compare_models <- loo::loo_compare(loo_zero, loo_negbinomial)
print(compare_models)



################################################################################
## Spatial Autocorrelation ######################################################
library(brms)

## Adding spatially correlated error term to model #############################
# Calculate an adjacency matrix using the coordinates
# First create a spatial coordinates matrix 
coords_all <- cbind(PMtracks_recov_vole_noNA$X, PMtracks_recov_vole_noNA$Y)

library(sp)
#crs <- CRS("+init=epsg:32632")
# Define the CRS of your original coordinates
original_crs <- CRS("+proj=utm +zone=32 +ellps=WGS84")

# Create a SpatialPoints object with your original coordinates
original_coords <- SpatialPoints(coords_all, proj4string = original_crs)

# Define the target CRS for decimal degrees (WGS84)
target_crs <- CRS("+proj=longlat +datum=WGS84")

# Transform the coordinates to decimal degrees
decimal_coords <- spTransform(original_coords, target_crs)

# Adjacency matrix - dist
library(geosphere)
dist <- distm(decimal_coords, decimal_coords, fun = distVincentyEllipsoid)

#?distm
# head(dist)
# Name the columns and rows
rownames(dist) = PMtracks_recov_vole_noNA$Linenr
colnames(dist) = PMtracks_recov_vole_noNA$Linenr
#head(dist)


# within 100km2 circle; neighbours = 1
radius <- sqrt(100000000 / pi)
radius #5641.896m radius

#Change the distance that are higher than 5641.896m to 0 
# Change distances greater than 2820.95 to 0
dist[dist > 5641.896] <- 0
#Check
max_vector <- as.vector(dist)
max(max_vector)

#Change the distance than are not zero to 1
dist[dist!= 0] <- 1
#Check
max_vector <- as.vector(dist)
max(max_vector)

#Negative binomial with adjacency matrix #######################################
library("brms")
PM_Abund_brms_nb_dist <- brms::brm(
  formula = Marten ~ 1 + reSpruce + reAgri + reSnow + reRoAbund + reElev + reElev*reRoAbund + (1|Year) + offset(log(length * SnowDays)) + car(dist, gr = Linenr, type = "icar"),
  data = PMtracks_recov_vole_noNA,
  data2 = list(dist = dist),
  family = brmsfamily("negbinomial"),
  prior = c(
    prior(normal(0, 2), class = "b"),   # Prior for all slope parameters
    prior(normal(0, 2), class = "Intercept")  # Prior for the intercept
  ),
  control = list(adapt_delta = 0.99),
  iter = 6000
)

#?car()

saveRDS(PM_Abund_brms_nb_dist, file = "PM_Abund_brms_nb_290324_dist.rds")



# LOOIC comparison between nb model with and without adjacency matrix ##########
library(loo)

PM_Abund_brms_nb_dist <- readRDS("PM_Abund_brms_nb_290324_dist.rds") # model with adjacency matrix
PM_Abund_brms_nb <- readRDS("PM_Abund_brms_nb_290324.rds") # model without adjacency matrix

# Compute the LOOIC for each model
loo_dist <- loo(PM_Abund_brms_nb_dist)
loo_dist
loo_negbinomial <- loo(PM_Abund_brms_nb)
loo_negbinomial


# Compare the models using the LOOIC
compare_models <- loo::loo_compare(loo_dist, loo_negbinomial)
print(compare_models)
# if the difference is more than 4 between elpd of differnet models, then you want to look at SE of the elpd
# -12.7/5.8 = 2.19
# The difference is small wrt to the standard error of the difference (< 4 times)
# In our case we have elpd_diff > 2 * se_diff, so we cannot be confident that the model with spatial error term is better.


# EXTRA 
# testSpatialAutocorrelation
library(DHARMa)
# Load data
PMtracks_recov_vole_noNA <- readRDS("PMtracks_recov_vole_noNA_290324.rds")
PM_Abund_brms_nb <- readRDS("PM_Abund_brms_nb_290324.rds")

# Generate posterior predictive samples for each year subgroup
ppd_by_year <- lapply(unique(PMtracks_recov_vole_noNA$Year), function(year) {
  newdata <- subset(PMtracks_recov_vole_noNA, Year == year)
  posterior_predict(PM_Abund_brms_nb, newdata = newdata, ndraws = 10000, summary = FALSE)
})

#only for year 3 (4th position in ppd_by_year) - change year to test for spatial autocorrelation in other years
ppd_yr <- ppd_by_year[4] 
ppd_yr_matrix <- do.call(cbind, ppd_yr)
ppd_yr_matrix <- t(ppd_yr_matrix)

# Filter the counts column based on the year group selected
PMtracks_recov_vole_noNA_yr <- PMtracks_recov_vole_noNA[PMtracks_recov_vole_noNA$Year == 3, "Marten"] 
PMtracks_recov_vole_noNA_X <- PMtracks_recov_vole_noNA[PMtracks_recov_vole_noNA$Year == 3, "X"] 
PMtracks_recov_vole_noNA_Y <- PMtracks_recov_vole_noNA[PMtracks_recov_vole_noNA$Year == 3, "Y"] 

res <- createDHARMa(simulatedResponse = ppd_yr_matrix, 
                    fittedPredictedResponse = apply(ppd_yr_matrix, 1, median), 
                    observedResponse = PMtracks_recov_vole_noNA_yr, 
                    integerResponse = T)

#D_sp <- recalculateResiduals(res, group = PMtracks_recov_vole_noNA_yr1$Linenr)
testSpatialAutocorrelation(res, x= PMtracks_recov_vole_noNA_X, y = PMtracks_recov_vole_noNA_Y)

# DHARMa Moran's I test for distance-based autocorrelation, p-value = 3.056e-06