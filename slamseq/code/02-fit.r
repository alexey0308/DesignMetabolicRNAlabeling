library(pulseR)
library(purrr)

rdsDir <- "slamseq/rds"
rawCounts <- readRDS(file.path(rdsDir, "counts.rds"))
samples   <- readRDS(file.path(rdsDir, "samples.rds"))
samples <- samples[c(4,1:3)]

totals <- map(unique(samples$id), function(.id) {
  i <- samples$desc[samples$id == .id]
  apply(rawCounts[,i],1,sum)
})
totals <- do.call(cbind, totals)

norms <- pulseR:::findDeseqFactorsSingle(totals)
norms <- rep(norms, 2)
names(norms) <- samples$desc

## create a list with options for the pulseR fitting
## (parameter boundaries, stopping criteria - by changes in logLik)
makeOpts <- function() {
  opts <- setFittingOptions(verbose = "verbose")
  opts <- setBoundaries(options = opts,
    list(
      mu1 = c(log(.1),log(1e6)),
      mu2 = c(log(.1),log(1e6)),
      d = c(1e-4,5),
      size = c(1, 1e5)
    )
  )
  opts <- setTolerance(
    params = 10, shared = 10, logLik = .01, options = opts)
  ## we do not fit normalisation coefficients, because they
  ## are derived from the DESeq-like normalisation
  opts$fixedNorms <- TRUE
  opts$replicates <- 3
  opts
}

## define the kinetic model
makeFormulae <- function() {
  list(
    formulas = MeanFormulas(
      labelled =  exp(mu2) * exp(-d*chase_time),
      unlabelled = exp(mu1) + exp(mu2) * (1 - exp(-d*chase_time))),
    formulaIndexes = list(
      unlabelled = "unlabelled",
      labelled = "labelled"
    ))
}

## prepare a PulseData object for the pulseR package
makePD <- function(counts, samples) {
  forms <- makeFormulae()
  pd <- PulseData(
    counts[, samples$desc],
    samples,
    forms$formulas,
    forms$formulaIndexes,
    groups=~chase_time+fraction
  )
  pd$depthNormalisation <- norms[samples$desc]
  pd
}

## first guess for the fitting
## I used a fit for the whole data set as a starting point
## for all other subset of samples
init <- function(counts) {
  fit <- list(
    mu1 = log(counts[,1]/10),
    mu2 = log(counts[,1]/10),
    d = rep(.5, length(counts[,1])),
    size = 1e5)
  ## fit <- readRDS("slamseq/rds/fit.rds")
  fit
}


## run pulseR
getFit <- function(counts, samples) {
  pd <- makePD(counts, samples)
  opts <- makeOpts()
  initPars <- init(counts)
  fit <- fitModel(pd, initPars, opts)
  list(fit=fit, pd=pd, opts=opts)
}

timeSets <- list(
  c(0, .5),
  c(0, .5, 1),
  c(0, 6),
  c(0, 12),
  unique(samples$chase_time),
  c(0, 3),
  c(0, 24)
)

fitTimePoints <- function(counts, samples, tp) {
  tpsamples <- samples[samples$chase_time %in% tp,]
  tpcounts <- counts[, tpsamples$desc]
  res <- getFit(tpcounts, tpsamples)
  saveRDS(
    res,
    sprintf("slamseq/rds/fit-%s.rds",
      paste(tp, collapse ="-")))
  res
}

res <- map(
  timeSets,
  fitTimePoints,
  counts=rawCounts, samples=samples)


x <- map(res, c("fit","d")) %>%
  do.call(what = cbind)

plot(as.data.frame(x))

res <- map(dir("slamseq/rds", "fit-0"), ~readRDS(file.path("rds", .x)))

## compute confidence intervals for d
ciss <- map(
  res,
  function(.x) {
    .x$opts$cores <- 3
    .x$opts$replicates <- 2
    cis <- ciGene("d",
      par=.x$fit,
      geneIndexes = seq_along(.x$fit$d),
      pd = .x$pd,
      options = .x$opts)
    tp <- unique(.x$pd$conditions$chase_time)
    saveRDS(cis,
      sprintf("slamseq/rds/cis-%s.rds",
        paste(tp, collapse ="-")))
    cis
  })
