library(pulseR)

## return formulas and indexes
getPureModel <- function() {
  model <- list()
  model$formulas <- MeanFormulas(
    total = exp(mu),
    labelled = exp(mu) * (1 - exp(-d * time)),
    unlabelled = exp(mu - d * time))
  model$formIndexes <- list(
    total_fraction = "total",
    flow_through = "unlabelled",
    pull_down = "labelled")
  model
}

createConditions <-
  function(fractions = c("total_fraction", "pull_down", "flow_through"),
           times = c(0,2,4,8),
           batches = c("a", "b")) {
    expand.grid(
      fraction = fractions,
      time = times,
      batch = batches
    )
  }
