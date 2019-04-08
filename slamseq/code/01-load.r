## Here the csv and Excel files from [Herzog et. al]
## are read into R objects and  saved as .rds files

library(cowplot)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(readxl)

rdsDir <- "slamseq/rds/pulse_chase"
dir.create(rdsDir, showWarnings = FALSE)


countDir <- "slamseq/data/pulse_chase/"
countFiles <- dir(countDir, "tsv")

selectedTx <- readxl::read_excel(
  file.path(countDir, "nmeth.4435-S4.xls"))
selectedTx <- selectedTx %>%
  mutate(id=paste(Chromosome, Start, End, sep="-"))


counts <- map(
  file.path(countDir, countFiles),
  read.table, head = TRUE)

## all names are the same
all(map_lgl(
  map(counts, "Name"),
  ~ identical(.x, y = counts[[1]]$Name)))

for(i in seq_along(counts)) {
  counts[[i]] <- counts[[i]] %>%
    mutate(id=paste(Chromosome, Start, End, sep="-"))
}
all(map_lgl(
  map(counts, "id"),
  ~ identical(.x, y = counts[[1]]$id)))

info <- strsplit(countFiles, "_")
## all same length
unique(map(info,length))

info <- data.frame(
  id = map_chr(info, 1),
  desc = map_chr(info, 5)
)

info$chase_time <-
  str_match(info$desc, "-([0-9\\.]+)h-")[,2] %>% as.numeric()
names(counts) <- info$desc

## there are duplicates
counts[[1]] %>%
  group_by(id) %>%
  filter(row_number() > 1) %>% dim

## take the first from the duplicates
for(i in seq_along(counts)) {
    counts[[i]] <- counts[[i]] %>%
        filter(id %in% selectedTx$id) %>%
        group_by(id) %>%
        filter(row_number() == 1) %>% ungroup()
}

totalCounts <- map(counts, "ReadCount") %>% bind_cols
labelledCounts <- map(counts, "TcReadCount") %>% bind_cols
unlabelledCounts <- totalCounts - labelledCounts
all(unlabelledCounts >= 0)

names(unlabelledCounts) <- paste0("U_", names(unlabelledCounts))
names(labelledCounts) <- paste0("L_", names(labelledCounts))

samples <- rbind(
  info %>% mutate(
    desc = paste0("U_", desc),
    fraction = "unlabelled"),
  info %>% mutate(
    desc = paste0("L_", desc),
    fraction = "labelled")
) %>% filter(!is.na(chase_time))
samples

##
heatmap(cor(totalCounts), symm = T)

rawCounts <- cbind(labelledCounts, unlabelledCounts)

rownames(rawCounts) <- as.character(counts[[1]]$id)

expressionLevel <- apply(totalCounts[1:6], 2, mean)
qplot(x = log10(expressionLevel + 1), geom = "density")

saveRDS(rawCounts, file.path(rdsDir, "counts.rds"))
saveRDS(samples, file.path(rdsDir, "samples.rds"))
saveRDS(selectedTx, file.path(rdsDir, "selectedTx.rds"))
