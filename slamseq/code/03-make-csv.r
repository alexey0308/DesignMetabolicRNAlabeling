## Just read the rds files from fiting results and save as csv

library(purrr)
library(cowplot)
library(stringr)
library(dplyr)

extractTimesFromNames <- function(files) {
  timeSets <- str_match_all(files, "-([0-9\\.].+)\\.rds") %>%
    map(2)
  timeSets
}

points2label <- function(times)
  paste(unique(times), collapse = "-")

rdsDir <- "slamseq/rds"

res <- map(dir(rdsDir, "fit-0"), ~readRDS(file.path(rdsDir, .x)))
timeSets <- extractTimesFromNames(dir(rdsDir, "fit-0"))
names(res) <- timeSets

ref <- readxl::read_excel(
  "./slamseq/data/pulse_chase/nmeth.4435-S4.xls")
ref$id <- paste(ref$Chromosome, ref$Start, ref$End, sep = "-")
ref <- ref %>%
  group_by(id) %>% filter(row_number() == 1) %>% ungroup

stopifnot(identical(timeSets,
  extractTimesFromNames(dir(rdsDir, "cis-0"))))

ci <- map(dir(rdsDir, "cis-0"),
  ~readRDS(file.path(rdsDir, .x)))
names(ci) <- timeSets

head(ci[[1]])

dat <- map2(res, ci, function(.x ,cis) {
  x <- data.frame(
    d = .x$fit$d,
    tau =  1 / .x$fit$d,
    tau_all = 1 / allPointsFit$fit$d,
    mu1 = exp(.x$fit$mu1),
    mu2 = exp(.x$fit$mu2),
    points = points2label(.x$pd$conditions$chase_time)
  )
  colnames(cis) <- c("d.min", "d.max")
  rownames(x) <- rownames(.x$pd$counts)
  x <- cbind(x, cis)
  x
})

csvDir <- "slamseq/csv"
dir.create(csvDir, showWarnings = FALSE)

iwalk(dat,
  ~ write.csv(.x, file = file.path(csvDir, paste0(.y,".csv")))
)
