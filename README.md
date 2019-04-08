## Supplementary code

> On the optimal design of metabolic RNA labeling experiments.  
> Alexey Uvarovskii, Isabel Naarmann-de Vries, Christoph Dieterich  
> bioRxiv 428862; doi: [https://doi.org/10.1101/428862](https://doi.org/10.1101/428862)

## wxMaxima code for derivations

[.mac file](maxima.mac)

## Requirements

```{r}
pkgs <- c(
  "tidyverse",
  "cowplot"
)

install.packages("devtools")
library(devtools)

install_github("dieterich-lab/pulseR", subdir="pkg")
install_cran(pkgs)

```

## Run

```{bash}

## compute confidence interval widths
## in simulations for the figure 2
make simulate-time-point-sets

## create figure 2
make figure-2

## fit the model to experimental data
make fit-experiment

## estimate confidence intervals for the parameters
make compute-conf-int

## create figure 6
make figure-6

## run slamseq analysis
make prepare-slamseq-data
make fit-slamseq

```
