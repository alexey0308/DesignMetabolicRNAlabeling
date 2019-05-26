library(purrr)
library(cowplot)
library(stringr)
library(dplyr)

extractTimesFromNames <- function(files) {
  timeSets <- str_match_all(files, "-([0-9\\.].+)\\.rds") %>%
    map(2)
  timeSets
}

naToVal <- function(x, val) {
  x[is.na(x)] <- val
  x
}

naToBorders <- function(x, vals) {
  for(i in 1:2)
    x[,i] <- naToVal(x[,i], vals[i])
  x
}

widths <- function(x)
    x[,2] - x[,1]
relWidth <- function(x,d)
    widths(x)/d

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

allPointsFit <- res[["0-0.5-1-3-6-12-24"]]
allPointsCIs <- ci[["0-0.5-1-3-6-12-24"]]

joinRefWithEstimation <- function(ref, x, xcis) {
  d <- data.frame(
    id = rownames(x$pd$counts),
    d = x$fit$d,
    mu1 = x$fit$mu1,
    mu2 = x$fit$mu2,
    stringsAsFactors = FALSE)
  ref <- ref[match(d$id, ref$id),]
  d <- cbind(d, ref)
  d <- cbind(d, xcis)
  d
}


plotCIOrdered <- function(fit, cis) {
  o <- order(fit$d)
  cis <- naToBorders(cis, c(1e-4,5))
  d <- data.frame(
    d = fit$d[o],
    min=cis[o,1],
    max=cis[o,2],
    i = seq_along(cis[,1])
  )
  ggplot(data = d) +
    geom_linerange(aes(
      ymin = min,
      ymax = max,
      x = i
    ),
    colour = "dodgerblue") +
    geom_point(aes(
      x = i, y = d),
      colour = "grey90") +
    scale_y_log10() +
    xlab("order by rate (from slow to fast)") +
    ylab(expression(paste(delta, ",  ", hr^-1))) +
    theme(axis.title.x = element_text(size = 12))

}

dash2comma <- function(x)
  str_replace_all(x, "-", ", ")

plotFigure3 <- function(res, ci) {
  setsToPlot <- c("0-0.5", "0-3", "0-6", "0-12")
  qs <- map(setsToPlot,
            ~ plotCIOrdered(res[[.x]]$fit, ci[[.x]]) +
              ggtitle(paste(dash2comma(.x), "hr")))
  q <- plot_grid(plotlist=qs)
  q
}

qRelCIOrdered <- plotFigure3(res, ci)

figureDir <- "figure/slamseq"
dir.create(figureDir, recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path(figureDir, "ci.pdf"),
  qRelCIOrdered,
  width = 174,
  height = 160,
  units = "mm",
  device = cairo_pdf)
ggsave(
  file.path(figureDir, "ci.png"),
  qRelCIOrdered,
  width = 174,
  height = 130,
  units = "mm",
  dpi = 400)

qCITwoPointsTogether <- plotCIOrdered(
  res[["0-0.5-1"]]$fit,
  ci[["0-0.5-1"]]
) + ggtitle("0, 0.5, 1 hr")
ggsave(
  file.path("figure/supplementary/ci-0-0.5-1.pdf"),
  qCITwoPointsTogether,
  width = 174/2,
  height = 80,
  units = "mm",
  device = cairo_pdf)
ggsave(
  file.path("figure/supplementary/ci-0-0.5-1.png"),
  qCITwoPointsTogether,
  width = 174/2,
  height = 80,
  units = "mm",
  dpi = 400)

dat <- map2(res, ci, function(.x ,cis) {
  data.frame(
    tau =  1 / .x$fit$d,
    tau_0 = 1 / allPointsFit$fit$d,
    mu1 = .x$fit$mu1,
    mu2 = .x$fit$mu2,
    relDeltaNA = relWidth(cis, .x$fit$d),
    relDeltaNARef = relWidth(cis, allPointsFit$fit$d),
    relDelta = relWidth(naToBorders(cis,c(1e-5,1e2)),
       .x$fit$d),
    relDeltaRef = relWidth(naToBorders(cis,c(1e-5,1e2)),
      allPointsFit$fit$d),
    points = points2label(.x$pd$conditions$chase_time)
  )
})
dat <- bind_rows(dat)

mu2Quantiles <- quantile(allPointsFit$fit$mu2, c(.4,.6))
setsRelCIvsTau <- c(
  "0-0.5-1-3-6-12-24",
  "0-12",
  "0-6"
)

qRelCIvsTau <- ggplot(
  data = dat %>% filter(
    points %in% setsRelCIvsTau,
    mu2 > mu2Quantiles[1],
    mu2 < mu2Quantiles[2]
  )) +
  geom_point(
    aes(x = tau,
      y = relDeltaNA,
      colour = points),
    size = .3,
    ##shape = ".",
    alpha = .6
  ) +
  xlab(expression(paste(tau, ", ",hr))) +
  ylab("relative CI (95%) ") +
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 5))) +
  annotation_logticks(sides = 'l') +
  coord_cartesian(
    xlim = c(1, 20),
    ylim = c(.2,2)
  ) +
  scale_color_brewer(
    palette = "Set1",
    breaks = setsRelCIvsTau,
    labels = paste(dash2comma(setsRelCIvsTau), " hr")
  ) +
  scale_y_log10() +
  scale_x_log10() +
  theme(legend.position = c(0.03, 1.00),
    legend.justification = c(0, 1),
    legend.text = element_text(margin = margin(l = 5)),
    legend.title = element_blank())


I_slam_mu0 <- function(pars) {
  eval(expression(
    ## I_slam_labeled_delta_delta
    k*mu2*t^2*exp(-d*t) /
      (mu2*exp(-d*t)+k) +
      ## I_slam_unlabeled_delta_delta
      (k*mu2^2*t^2*exp(-2*d*t)) /
      ((mu2*(1-exp(-d*t))+mu1)*(mu2*(1-exp(-d*t))+mu1+k)))
  , pars)
}

getOptTime <- function(k, mu2, mu1=2e2) {
  f <- function(x) {
    pars <- list(d = 1, t = x, mu2 = mu2, k = k, mu1 = mu1)
    I_slam_mu0(pars)
  }
  optimize(f, c(1e-2,20), maximum = TRUE)
}

optRatio <- getOptTime(allPointsFit$fit$size,
  mu2 = exp(median(allPointsFit$fit$mu2)),
  mu1 = exp(median(allPointsFit$fit$mu1)))


exp(median(allPointsFit$fit$mu2))
exp(median(allPointsFit$fit$mu1))
(summary(1/allPointsFit$fit$d))
2.9/(summary(allPointsFit$fit$d))
allPointsFit$fit$size

getIvsTime <- function(pars, times) {
  Is <- map_dbl(times, function(x) {
    pars$t <- x
    I_slam_mu0(pars)
  })
  data.frame(time = times, I = Is)
}

x <- getIvsTime(
  list(d = 1,
    mu2 = median(exp(allPointsFit$fit$mu2)),
    mu1 = median(exp(allPointsFit$fit$mu1)),
    k = allPointsFit$fit$size),
  seq(.1,10,.1))

qInformationVsTime <- qplot(
  data = x,
  x = time,
  y = I,
  geom = "line",
  log='xy',
  xlab = expression(paste(t, "/", tau)),
  ylab = expression(paste(I[paste(delta, delta)] %.% delta^2))) +
  theme(axis.title = element_text(family = "serif"))

qInformationVsTime <-
qInformationVsTime +
  annotate("point",
    x = optRatio$maximum,
    y = optRatio$objective, shape = 1, size = 3,
  colour = "dodgerblue") +
  annotate("segment",
    x = optRatio$maximum,
    xend = optRatio$maximum,
    y = 0,
    yend = optRatio$objective,  size = .5,
    colour = "dodgerblue") +
  scale_x_continuous(
    trans = "log",
    breaks = c(.1,.5,1, optRatio$maximum, 5, 10),
    labels = c("0.1","0.5", "1",
      format(optRatio$maximum, digits=2), "5", "10"))

optimalTimeText <- format(optRatio$maximum, digits=2)

qRelCIvsTauAnnotated <-
  qRelCIvsTau +
  annotate("segment",
      x = c(6,12)/optRatio$maximum,
      xend = c(6,12)/optRatio$maximum,
      y = c(.1,.1),
      yend = c(.35,.35),
      arrow = arrow(length=unit(.05, "npc"),
        angle = 15)
  ) +
  scale_x_continuous(
    trans = "log",
    breaks = c(.1,.5,1,
      6/optRatio$maximum,
      12/optRatio$maximum,
      10),
    labels = c("0.1","0.5", "1",
      "6/2.9",
      "12/2.9",
      "10"))
qRelCIvsTauAnnotated


## plot figure 4
q <- plot_grid(
  qInformationVsTime,
  qRelCIvsTauAnnotated,
  rel_widths = c(1,1),
  nrow = 1,
  labels = c("A", "B")
)
ggsave(
  file.path(figureDir, "ci-time.pdf"),
  q,
  width = 174,
  height = 80,
  units = "mm",
  device = cairo_pdf)
ggsave(
  file.path(figureDir, "ci-time.png"),
  q,
  width = 174,
  height = 80,
  units = "mm",
  dpi = 400)


## correlation between single-point and all-point estimates

plotSupplCorrelation <- function() {
  dCor <- map(res, ~.x$fit$d) %>%
    bind_cols() %>%
    cor( method = "spearman")

  ii <- c(
    "0-0.5",
    "0-3",
    "0-6",
    "0-12",
    "0-24"
    #  "0-0.5-1-3-6-12-24"
  )
  allPointsIndex <- "0-0.5-1-3-6-12-24"

  qCor <- qplot(
    y = dCor[allPointsIndex,ii],
    x = factor(ii, levels = ii),
    ) +
    geom_segment(aes(
      x = factor(ii, levels = ii),
      y = 0,
      xend = factor(ii, levels = ii),
      yend = dCor[allPointsIndex,ii])) +
    xlab("points") +
    ylim(c(0,1)) +
    ylab("correlation") +
    ggtitle(expression(
      paste( delta[subset], " vs ", delta["all points"])
    ))
  qCor
}

qCor <- plotSupplCorrelation()
ggsave(
  file.path("figure/supplementary/correlation.pdf"),
  qCor,
  width = 174/2,
  height = 80,
  units = "mm",
  device = cairo_pdf)
