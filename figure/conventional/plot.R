## figure 6

library(cowplot)
library(dplyr)
library(purrr)

source(file.path("figure", "funs.R"))
source(file.path("experiment", "funs.R"))

rdsDir <- "experiment/rds"
r <- readFits(rdsDir)

## calculate I_dd for a given condition and time
## for all the genes, which are described with the `fit`
getIFromData <- function(fit, normCoef, t, condition) {
  if (condition == "total_fraction") return (NULL)
  fun <-c("pull_down" = I_pulse, "flow_through" = I_chase)[condition]
  res <- data.frame(
    I = fun[[1]](list(
      mu = exp(fit$mu)* normCoef,
      d = fit$d,
      k = fit$size,
      t = t
    )),
    t = t,
    condition = c("pull_down" = "labeled",
                  "flow_through" = "unlabeled")[condition],
    stringsAsFactors = FALSE
    )
  res$id <- rownames(res)
  res$d <- fit$d
  res$mu <- fit$mu
  res
}

## use internal pulseR function to get the normalisation factors
## calculate the normalisation coefficient for every sample
norms <- pulseR:::getNorms(r$pd, r$fit$normFactors)
norms <- apply(norms, 2, max)
## apply getIFromData to every sample:
## i.e. for every set of norm, time, condition used in the experiment.
## construct a table with I_dd for every gene+sample
o <- pmap(list(normCoef = as.list(norms),
               t = as.list(r$pd$conditions$time),
               condition = as.character(r$pd$conditions$condition)),
          getIFromData,
          fit = r$fit)
o <- bind_rows(o)

plotFIForGenes <- function(o) {
  q <- ggplot(data = o) +
    geom_point(
      aes(
        x = d,
        y = d^2*I,
        colour = condition,
        ),
      stroke = 0,
      size = .18
    ) +
    facet_wrap(~ t,
               scales = "free_y",
               labeller = as_labeller(function(x) paste(x, " hr"))) +
    ylab(
     expression(paste(I[paste(delta, delta)] %.% delta^2))) +
    xlab(expression(paste("d, ", hr^-1))) +
    scale_x_log10() +
    scale_y_continuous(trans = "sqrt") +
    scale_color_brewer(palette = "Set1", name = "") +
    guides(colour = guide_legend(override.aes = aes(size = 1))) +
    annotation_logticks(sides = 'b') +
    geom_line(data = o %>% filter(d < .71),
              aes(x = d,y = t^2*r$fit$size*d^2),
              linetype = "42") +
    geom_line(data = o %>% filter(d < 2.3),
              linetype = "42",
              aes(x = d,y = t^2*r$fit$size*d^2 * exp(-2*d*t) / (1 - exp(-d*t))^2)) +
    theme(strip.background = element_rect(fill = "white"),
          legend.position = c(.91,1),
          legend.justification = c(1,1),
          axis.title.y = element_text(family = "serif"),
          text = element_text(size = 10),
          axis.text = element_text(size = 10))
  q
}

fiForGenes <- plotFIForGenes(o)


## a helper to keep same depth in simulation:
## if in the total sample the depth is sum(exp(mu)),
## estimate the coefficient for the pull-down and flow-through
## which result in the same total expected read count in the
## corresponding sample
calculateNormCoeff <- function(t, fit) {
  unlab <- exp(fit$mu) * exp(-fit$d*t)
  lab <- exp(fit$mu) * (1 - exp(-fit$d*t))
  list(
    pull_down = sum(exp(fit$mu)) / sum(lab),
    flow_through = sum(exp(fit$mu)) / sum(unlab)
  )
}

## add genes to the data set with the expected read count in the
## total fraction of 1000 reads. The degradation rates for slow and fast genes
## correspond to the 1e-3 and 1 - 1e-3 = .999 quantiles
plotFastSlow <- function() {
  fit <- r$fit
  fi <- lapply(c(fit$size, 1e10), function(k) {
    fit$mu <- c(
      set_names(rep(log(1e2), 3), c("FAST", "SLOW", "MEDIAN")),
      fit$mu)
    fit$d <- c(
      set_names(quantile(fit$d, c(1 - 1e-3, 1e-3, .5), type = 1),
                c("FAST", "SLOW", "MEDIAN")),
      fit$d)
    fit$size <- k
    fi <- lapply(
      exp(seq(log(1e-3), log(1000), length.out = 100)),
      function(t) {
        coefs <- calculateNormCoeff(t, fit)
        fi <- imap(coefs,
                   ~ getIFromData(fit, normCoef = .x, t = t, condition = .y)[1:3,])
        fi <- bind_rows(fi)
        fi
      })
    fi <- bind_rows(fi)
    fi$k <- fit$size
    fi
  })

  fi <- bind_rows(fi)
  fi$overdispersion <- "high"
  fi$overdispersion[fi$k > 1e4] <- "low"
  fi <- filter(fi, id != "MEDIAN")

  q <- ggplot(data = fi %>% filter(I*d^2 > 1e-1, d*t > 1e-3)) +
    geom_line(aes(x = t,
                  y = d^2*I,
                  colour =  condition,
                  linetype = overdispersion)) +
    scale_color_brewer(palette = "Set1", name = "") +
    scale_x_log10() +
      scale_y_log10() +
    facet_wrap(~id,
               scales = "free_x"
               ) +
    ylab(
      expression(paste(I[paste(delta, delta)] %.% delta^2))) +
    annotation_logticks(sides = 'b') +
    xlab("time, hr") +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    theme(strip.background = element_rect(fill = "white"),
          text = element_text(size = 10),
          axis.title.y = element_text(family = "serif"),
          axis.text = element_text(size = 10))
  q
}

fastSlow <- plotFastSlow()

q <- plot_grid(
  fiForGenes,
  fastSlow,
  nrow = 2,
  labels = c("A", "B")
)
q

ggsave(
  filename=file.path("figure", "conventional", "figure.png"),
  plot = q,
  dpi = 600,
  width = 174,
  height = 130,
  units = "mm")
