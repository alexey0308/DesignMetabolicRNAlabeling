library(pulseR)
source("experiment/funs.R")
source("experiment/models.R")

## Here we use the fitted model to estimate the 95% confidence intervals (CIs)
## using the pulseR package (the profile likelihood CIs).
## The procedure is performed for both, mu and d parameters.
## The CIs are saved to the cis.rds file as a list with an individual
## table with 2 columns (min, max values) for every parameter.

## read the fitting results
rdsDir <- "experiment/rds"
path <- getPath(rdsDir)
opts <- readRDS(path$opts)
par  <- readRDS(path$fit)
pd   <- readRDS(path$pd)

## compute profile likelihood CI
## start optimization 5 times from randomly perturbed parameter values for every gene
opts$replicates <- 5
cis <- lapply(c(d = "d", mu = "mu"), function(parName)
  ciGene(parName, seq_along(par$mu), pd, par, options = opts))
saveRDS(cis, path$cis)
