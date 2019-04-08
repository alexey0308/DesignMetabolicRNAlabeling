library(cowplot)
library(dplyr)
library(purrr)

source(file.path("figure", "funs.R"))

## Fisher information dd term vs time for d = 1, mu = 1
FIMvsTime <- function() {
  dat <- lapply(c(I_pulse, I_chase, I_slam), function(f) {
    t <- exp(seq(log(.01), log(10), length.out = 1000))
    data.frame(time = t,
               I = f(list(t = t, d = 1, mu = 1), NB = FALSE))
  })
  names(dat) <- c("labeled", "unlabeled", "both")
  dat <- dat %>% bind_rows(.id = "fraction")
  dat$fraction <- factor(dat$fraction, c("labeled", "unlabeled", "both"))

  mu <- 10000

  q <- qplot(
    data   = dat,
    x      = log2(time),
    y = I,
    geom   = "line",
    size   = I(1.5),
    colour = fraction,
    log    = 'y',
    xlab   = expression(paste(log[2],"(",t ,"/" , tau, ")")),
    ylab   =
      expression(paste(I[paste(delta, delta)]))) +
    coord_cartesian(xlim = c(-4, 4), ylim = c(.01, 2)) +
    scale_color_brewer(palette = "Set1", name = "") +
    annotation_logticks(sides = "l") +
    theme(
      legend.position = c(0.1,1.15),
      legend.justification = c(0,1),
      axis.title.y = element_text(family = "serif"))
  q
}

## Confidence interval width (relative) of the degradation rates
## with a different set of time points included in the simulation.
CIvsSet <- function() {
  sims <- readRDS(
    file = file.path("simulations","rds", "ci-for-subsets.rds"))
  ## make table with cis
  cis <- map(sims, c("cis", "d"))
  cis <- lapply(cis, data.frame)
  cis <- bind_rows(cis, .id = "set")
  names(cis)[2:3] <- c("min", "max")
  cis$trueValue <- sims[[1]]$params$d
  cis$set <- factor(cis$set)
  cis$set <- factor(cis$set, levels(cis$set)[order(nchar(levels(cis$set)))])
  levels(cis$set)[4] <- gsub("\\|", "+", levels(cis$set)[4])
  pal <- scales::brewer_pal(pal="GnBu")(5)[-1]
  q <- ggplot() +
    geom_smooth(
      data = cis,
      aes(
        x = trueValue,
        y = (max - min) / trueValue,
        colour = set
      ),
      se = FALSE,
      size = 1.5
    ) +
    scale_x_log10() +
    annotation_logticks(sides = "l") +
    scale_y_log10() +
    xlab(expression(paste(delta, ", ", hr^-1))) +
    ylab("relative 95% CI") +
    theme(legend.position = c(0.29,1.05),
          legend.justification = c(0,1),
          legend.key = element_rect("white"),
          legend.title = element_text(size = 12)
          ) +
    guides(colour = guide_legend(reverse = TRUE)) +
    #scale_colour_brewer(palette = "GnBu", name = "time points")
    scale_colour_manual(values = pal)
  q
}


## Expression level saturation effect
saturationPlot <- function() {
  mus <- 10^seq(0,5,.01)
  qplot(
    x   = mus,
    y   = sqrt(I_slam_inv(list(d = 1, t = 1, mu = mus, k = 100))),
    ylab = expression(sd(delta)/delta),
    xlab = expression(mu),
    geom = "line",
    size = I(1)
  ) +
    geom_line(aes(
      x   = mus,
      y   = sqrt(I_slam_inv(list(d = 1, t = 1, mu = mus, k = 1e100)))),
      size = 1,
      linetype = 2
      ) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
}
saturationPlot()

##' effect of overdispersion on the dd diagonal term of the inverse of the FIM
inverseInfo <- function() {
  mu <- 100
  I.vs.k <- lapply(c(1, 10,1e2,1e8), function(k) {
    d <- 1
    ts <- exp(seq(log(1e-3), log(15), .01))
    mus <- mu
    data.frame(
      labeled = I_pulse(list(mu = mus, k = k, d = d, t = ts)),
      unlabeled = I_chase(list(mu = mus, k = k, d = d, t = ts)),
      both = I_slam(list(mu = mus, k = k, d = d, t = ts)),
      inverse = I_slam_inv(list(mu = mus, k = k, d = d, t = ts)),
      t = ts,
      d = d,
      mu = mus,
      k = k
    )
  })
  I.vs.k <- bind_rows(I.vs.k)
  pal <- scales::brewer_pal(pal="GnBu")(5)[-1]

  q <-
    I.vs.k %>%
    qplot(data = .,
          x = t,
          y = sqrt(inverse),
          colour = factor(k),
          log='xy',
          geom = "line",
          size = I(1.0),
          group = rev(k)) +
    xlab(expression(t/tau)) +
    ylab(expression(sd(delta)/delta)) +
    coord_cartesian(ylim = c(1e-1, 20)) +
   scale_colour_manual(
     values = rev(pal),
     name = "k",
     labels =
        scales::trans_format("log10", scales::math_format(10^.x))(unique(I.vs.k$k))
   ) +
    theme(strip.background = element_rect("white"),
          legend.position = c(.2,0.9),
          legend.justification = c(0,1),
          legend.title = element_text(margin = rep(unit(-.7, "cm"),4))
          )
  q
}
inverseInfo()


figure2 <- plot_grid(
  FIMvsTime(),
  CIvsSet(),
  inverseInfo(),
  saturationPlot(),
  labels = c("A", "B", "C", "D"),
  align = "hv"
)
figure2

ggsave(
  filename = file.path("figure", "2", "figure.pdf"),
  plot = figure2,
  width = 174,
  height = 160,
  device = cairo_pdf,
  units = "mm")
