##' a helper to get all rds result files
getPath <- function(rdsDir) {
  x <- c("pd", "fit", "opts", "cis")
  path <- lapply(setNames(x, x),
                    function(y)
                      file.path(rdsDir, paste0(y, ".rds")))
  path
}

##' 1) read conditions and read counts
##' 2) ensure the order of samples in the count table and sample
##'    sheet to be the same
##' 3) treat the total_fraction as the reference factor level
getPulseData <- function(rdsDir) {
  genes <- readRDS(file.path(rdsDir, "genes.rds"))
  conditions <- readRDS(file.path(rdsDir, "conditions.rds"))
  genes <- genes[conditions$sample]

  conditions <- conditions[c("fraction", "time")]
  names(conditions)[1] <- "condition"
  conditions$condition <- factor(conditions$condition)
  conditions$condition <-
    relevel(conditions$condition, "total_fraction")
  list(conditions = conditions,
       genes = genes)
}

##' get a logical vector to filter out low expressed genes
##' on the basis of the geometric mean
whichHigh <- function(x, level) {
  apply(x, 1, function(y) exp(mean(log(y)))) > level
}

##' The total fractions were biotinylated as well.
##' In result, the total fractions at t > 0 hr behave as flow-through ones,
##' because the RiboZero removed the labelled molecules together with the
##' ribosomal RNA (the protocol uses streptavidin beads and biotinylated oligos
##' against the ribo-RNA). The `conditions.rds` is already processed by this
##' function, I keep it here for the record.
##'
##' Returns a data frame with fixed labels of the samples.
biotinsAsFlowThrough <- function(data) {
  unlabelledIndex <-  (data$conditions$condition == "total_fraction" &
                       data$conditions$time > 0)

  data$conditions$condition <- as.character(data$conditions$condition)
  data$conditions$condition[unlabelledIndex] <- "flow_through"
  data$conditions$condition <- factor(data$conditions$condition)
  data$conditions$condition <- relevel(data$conditions$condition, "total_fraction")
  data
}

##' a helper to read fit results and to make a result table.
##' returns a list with
##'     result$pd  - PulseData object
##'     result$fit - fitting results
##'     result$tab - fitting results with confidence intervals as a data.frame
readFits <- function(d) {
  models <- c("pure", "contaminated")
  res <- readRDS(file.path(d, "fit.rds"))
  pd <- readRDS(file.path(d, "pd.rds"))
  fits <- res[c("mu", "d")]
  if (file.exists(file.path(d, "cis.rds"))) {
    cis <- readRDS(file.path(d, "cis.rds"))
    cis <- lapply(cis,
                  function(x) {
                    setNames(as.data.frame(x), c("min", "max"))
                  })
    cis <- do.call(cbind, cis)
    fits <- cbind(as.data.frame(fits), cis)
  }
  list(fit = res,
       pd = pd,
       tab = as.data.frame(fits))
}

##' helper to make the list of fitting options
setOpts <- function(bounds, tolerance, resultRDSpath) {
  opts <- setBoundaries(bounds, normFactors = bounds$normFactors)

  opts <- setTolerance(
    params = tolerance$params,
    normFactors = tolerance$normFactors,
    options = opts
  )
  opts$tolerance$logLik <- tolerance$logLik
  opts$resultRDS <- resultRDSpath
  opts$verbose <- "verbose"
  opts$cores <- cores
  opts
}
