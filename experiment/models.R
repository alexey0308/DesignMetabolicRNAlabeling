##' A generic helper to generate different model specifications
##' (fractions used, pure/contaminated).
##'
##' Returns formulas and formulaIndexes to use in `PulseData` function.
##'
createPulseModel <- function(fractions, contaminated) {
  formulas <- MeanFormulas(
    total      = exp(mu),
    labelled   = exp(mu) * (1 - exp(-d * time)),
    unlabelled = exp(mu - d * time)
  )
  formIndexes  <- list(
    total_fraction = "total",
    flow_through   = c("unlabelled", "labelled"),
    pull_down      = c("labelled", "unlabelled")
  )[fractions]
  if (!contaminated) {
    formIndexes <- lapply(formIndexes, `[[`, 1)
  }
  usedFormulas <- unique(unlist(formIndexes))
  list(formulas    = formulas[usedFormulas],
       formIndexes = formIndexes)
}
