## create supplementary tables with the gene counts and
## the estimations for the parameters (gene name is an ENSEMBL gene_id's).
## The model was fitted to a subset of genes prefiltered by expression level,
## hence parameter values are available only for this subset.

r <- readRDS("experiment/rds/genes.rds")
conditions <- readRDS("experiment/rds/conditions.rds")

sampleNames <- do.call(paste, c(conditions[,-1], list(sep = "-")))
names(r) <- sampleNames
r <- r[grep("^ENS", rownames(r)),]
r <- cbind(data.frame(gene_id = rownames(r)), r)

write.table(x = r,
            file = "supplementary/raw-counts.csv",
            sep = ",",
            row.names=FALSE)


source("experiment/funs.R")
library(tidyverse)
dat <- readFits("experiment/rds")
dat <- dat$tab
dat <- rownames_to_column(as.data.frame(dat), "gene_id")

dat <- dat %>%
  mutate_at(vars(starts_with("mu")), exp)
dat <- dat[, c(1:3, 6:7, 4:5)]
dat <- dat[grep("^ENS", dat$gene_id),]

write.table(x = dat,
            file = "supplementary/fit-values.csv",
            sep = ",",
            row.names=FALSE)
