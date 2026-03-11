# metapgm

The goal of `metapgm` is to provide tools for fitting Penalized Gaussian Mixture (PGM) models for meta-analysis. Unlike standard random-effects models that assume a single normal distribution of true effects, `metapgm` can flexibly capture complex, non-normal, and multimodal distributions. It includes functions for standard PGM models as well as conditional models with meta-regression capabilities.

## Installation

You can install the development version of metapgm from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("daihesui/metapgm")
```

## Example: Capturing a Bimodal Distribution

To see the power of the PGM approach, let's simulate a meta-analysis where there are two distinct sub-populations of effect sizes (e.g., studies using two very different interventions), creating a bimodal distribution. 

Standard models would fail to capture this nuance, but `metapgm` handles it effortlessly:

``` r
library(metapgm)

# 1. Simulate bimodal meta-analysis data
set.seed(42)

# Sub-population A (True effect centered around 0.2)
yi_A <- rnorm(70, mean = 0.2, sd = 0.1)
vi_A <- runif(70, min = 0.01, max = 0.03)

# Sub-population B (True effect centered around 0.8)
yi_B <- rnorm(30, mean = 0.8, sd = 0.1)
vi_B <- runif(30, min = 0.01, max = 0.03)

# Combine into a single dataset
yi_bimodal <- c(yi_A, yi_B)
vi_bimodal <- c(vi_A, vi_B)

# 2. Fit the Penalized Gaussian Mixture model
fit <- rma.pgm(yi = yi_bimodal, vi = vi_bimodal)

# 3. View the summary statistics
summary(fit)
```

The true power of the model is best seen visually. The `plot()` function overlays the PGM estimated density over the histogram of observed effect sizes, cleanly identifying the two underlying peaks:

``` r
# Plot the estimated density
plot(fit)
```
