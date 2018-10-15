## create a data set with several time points
##
## estimate confidence intervals for the degradation rates for different
## sets of time points used
##

set.seed(239)
source("simulations/sim.R")

library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(cowplot)

model <- getPureModel()

## Instead of making many replications I decided to implement
## everything in on run: we created `geneReplicates` number
## of identical gene in the count, hence replicateNumber = 1.
## 20 genes X 10 geneReplicates X 3 fractions X 3 points + 0hr total
simulationParams <- list(
  replicateNumber = 1,
  geneNumber      = 20,
  geneReplicates  = 10,
  expressionLevel = 10000,
  overdispersion  = 1e8, # huge size-parameteter --> minimal overdispersion
  timePoints      = c(0, 2, 4, 8)
)

batches <- letters[1:simulationParams$replicateNumber]

model$conditions <- createConditions(
  batches   = batches,
  fractions = names(model$formIndexes),
  times     = simulationParams$timePoints)
model$conditions <- model$conditions %>%
  filter((time == 0 & fraction == "total_fraction") |
         (time > 0  & fraction != "total_fraction"))

## generate parameters
## d is logarithimically equidistant
ranges <- list(d = 1/c(.1,200))
params <- list()
params$d <- exp(
  seq(length.out = simulationParams$geneNumber,
    from = log(ranges$d[1]),
    to = log(ranges$d[2]))
)
params$d <- rep(params$d, simulationParams$geneReplicates)
## we use logarithm of the expression, so mean read counts
## in the total sample is exp(mu)
params$mu <- rep(log(simulationParams$expressionLevel), length(params$d))
params$size <- simulationParams$overdispersion


## to mimick the SLAMseq read counts, all normalisation factors are equal 1,
## then labelled + unlabelled = total
## the normalisation factors are not fitted, because in this case we assume as
## labelled and unlabelled counts originate from the same sequencing sample,
## the simulation is done assuming the same sequencing depth.
normFactors <- list(
  total_fraction = 1,
  pull_down      = 1,
  flow_through   = 1
)


allNormFactors <- multiplyList(normFactors, model$conditions$fraction)

## manually generate: use internal pulseR function for creation of test data
counts <- pulseR:::generateTestDataFrom(
  model$formulas,
  model$formIndexes,
  allNormFactors,
  params,
  model$conditions)

## create sets of time points
z <- simulationParams$timePoints[-1]
times <- c(as.list(z),  list(z))

fitToSubset <- function(timePoints) {
  cat(timePoints, "\n")
  i <- model$conditions$time %in% c(0, timePoints)
  pd <- PulseData(
    counts[,i,drop = FALSE],
    model$conditions[i,],
    model$formulas,
    model$formIndexes,
    groups = ~ fraction + time )

  ## fit only gene parameters
  opts <- pulseR::setFittingOptions(verbose = "silent")
  opts <- pulseR::setTolerance(params = .01, options = opts)
  opts$tolerance$logLik <- .5
  opts <- pulseR::setBoundaries(
    list(d = c(1e-5,50),
      mu = (log(c(1,1e5))),
      size = c(1e7,1e8)))
  opts$replicates <- 3
  opts$cores <- 4

  ## this is very important: we do not fit normalisation in this simulations
  opts$fixedNorms <- TRUE

  fit <- fitModel(
    pulseData = pd,
    par       = params,
    options   = opts
  )
  cat("successfully fitted\n")

  cis <- lapply(
    c(#mu = "mu",
      d = "d"),
    function(parName) {
      pulseR::ciGene(
        parName=parName,
        geneIndexes=seq_along(fit$mu),
        pd = pd,
        par = fit,
        options = opts)
    })
  list(
    pd = pd,
    params = params,
    fit = fit,
    opts = opts,
    simulationParams = simulationParams,
    cis = cis)
}

sims <- lapply(times, fitToSubset)

timeNames <- unlist(lapply(times, paste, collapse = " | "))
names(sims) <- timeNames
saveRDS(sims, file = "simulations/rds/ci-for-subsets.rds")
