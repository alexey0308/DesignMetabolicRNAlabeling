library(pulseR)
source("experiment/funs.R")
source("experiment/models.R")

## Here we fit the kinetic model to the read counts
## For every gene there are two parameters:
##   - degradation rate d
##   - logarithm of expression level
##       in comparison to the main text, we work here with the logarithm of
##       the expression level, i.e. the mean read count in the total fraction
##       is exp(mu), and in the unlabelled one, it is exp(mu - d*t).
##       It is done to keep the parameters in the same scale during
##       the optimization and does not affect the likelihood function values.
##       In the final supplementary table with the parameter values, we use
##       the exp(mu) values, which represent expected read count, not the logarithm,
##       so the reader may see the expression level in usual units.
## In addition, there is the size parameter of the NB distribution
## (the dnbinom function).
## The normalization factors are shared between samples from the same
## fraction and time point. In addition, the DESeq-like normalization is applied
## within the group to account for different sequencing depth (implemented
## in the pulseR package).
##
## pulseR v1.0.2 was used

## fitting options: stopping criteria and parameter boundaries for optimization
tolerance <- list(params = 0.01,
                  normFactors = 0.01,
                  logLik = 0.1)
bounds <- list(mu = log(c(1e-1, 1e8)),
               d = c(1e-3, 6),
               size = c(1,1e3),
               normFactors = c(1e-6,20))
cores <- 4
filterExressionLevel <- 50

## input/output options
rdsDir <- "experiment/rds/"
dir.create(rdsDir, showWarnings = FALSE, recursive = TRUE)
rdsFiles <- getPath(rdsDir)

data <- getPulseData(rdsDir)
## total samples at t > 0 hr were biotinylated and labeled RNA
## was removed during library preparation by the RiboZero kit
data <- biotinsAsFlowThrough(data)

## use genes only with higher expression to reduce the data set
## using geometric mean of [counts + 1] from the total samples
highGenes <- whichHigh(
  1 + data$genes[, data$conditions$condition == "total_fraction"],
  filterExressionLevel)
counts <- data$genes[highGenes,]

## - specify the model: all three fraction types, no contamination
## - create a PulseData object
## - initialize first guess of parameters for optimization
## - save fitted values and the PusleData object to .rds files
fractions <- c(
  "total_fraction",
  "flow_through",
  "pull_down")
isContaminated <- FALSE
model <- createPulseModel(fractions, isContaminated)

pd <- PulseData(
  counts = counts,
  conditions = data$conditions,
  formulas = model$formulas,
  formulaIndexes = model$formIndexes,
  groups = ~ condition + time
)
saveRDS(rdsFiles$pd, object = pd)

opts <- setOpts(bounds, tolerance, rdsFiles$fit)
saveRDS(opts, rdsFiles$opts)

par <- initParameters(list(),
                      geneParams = c("mu", "d"),
                      pulseData = pd,
                      opts)
par$mu <- log(.1 + counts[, 1])
par$size <- 1e1
par$normFactors <- lapply(par$normFactors,
                          function(x) {x[1] <- 10; x})

par <- fitModel(pd, par, opts)

saveRDS(par, rdsFiles$fit)

