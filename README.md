The Flexible Gumbel Distribution
================

## Introduction

This is the GitHub repository for the paper titled “The Flexible Gumbel
Distribution: A New Model for Inference About the Mode” (**add the link
later**). In this paper, we introduce a new unimodal distribution
indexed by the mode. This distribution is a mixture of Gumbel
distributions for the maximum and minimum values, which we refer to as
the flexible Gumbel distribution.

The probability density function (pdf) of the flexible Gumbel
distribution is given by:

$$
\begin{aligned}
f(y)= & w \times \frac{1}{\sigma_1} \exp \left\{-\frac{x-\theta}{\sigma_1}-\exp \left(-\frac{x-\theta}{\sigma_1}\right)\right\}+ \\
& (1-w) \times \frac{1}{\sigma_2} \exp \left\{\frac{x-\theta}{\sigma_2}-\exp \left(\frac{x-\theta}{\sigma_2}\right)\right\},
\end{aligned}
$$

where $w \in[0,1]$ represents the mixing proportion parameter. We denote
that $Y$ follows the distribution specified by the pdf as
$Y \sim \operatorname{FG}\left(\theta, \sigma_1, \sigma_2, w\right)$.

## Sampling from the Flexible Gumbel Distribution

Sampling from the flexible Gumbel distribution is straightforward due to
its mixture distribution nature with two components. To generate $n$
independent and identically distributed (i.i.d.) samples from the
flexible Gumbel distribution, you can utilize the provided functions in
this GitHub repository. Below is an example of sampling 1000 i.i.d.
samples from $\operatorname{FG}(0,1,1,0.3)$.

``` r
source("./density.R")

set.seed(100)

# Set the sample size
n <- 1000

# Define parameters
w1 <- 0.2
loc <- 0
scale1 <- 1
scale2 <- 1

# Sample from the mixture distribution
Z <- rbinom(n=n, size=1, prob=w1)
X1 <- rgumbel_RNG(n=n, loc, scale1)
X2 <- 2 * loc - rgumbel_RNG(n=n, loc, scale2)
Y <- Z * X1 + (1 - Z) * X2

# Plot histogram and density
hist(Y, breaks=50, main="Distribution of FG Samples", freq=FALSE)
points(density(Y), type="l", col=rgb(0.8, 0, 0, 1))

# Add density curve of the true pdf
xaxis <- seq(min(Y)*1.2,max(Y)*1.2,length.out = 1000)
yaxis <- mixture_gumbel_pdf(x = xaxis, w1 = w1, loc = loc, 
                            scale1 = scale1, scale2 = scale2, log = FALSE)
lines(xaxis, yaxis, type="l", col=rgb(0, 0.8, 0, 1), lty=2)

# Add legend
legend("topright", legend=c("Sample Density", "True Density"), 
       col=c("red", "green"), lty=c(1, 2), cex=0.8)
```

<img src="README_files/figure-gfm/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />

## Frequentist Inference

The Frequentist inference of the flexible Gumbel (FG) distribution
includes an implementation of maximum likelihood estimation (MLE) and
the associated sandwich variance estimation to quantify the uncertainty
of the MLE. We utilized the Expectation–Conditional Maximization (ECM)
algorithm to calculate the MLE estimation of the FG distribution. An
example of how to use the ECM algorithm program is provided below.

``` r
source("./EM_v1.0.R")
estimation <- em(Y, w1, loc, scale1, scale2)
print(round(estimation,3))
```

    ##     w1    loc scale1 scale2 loglikelihood
    ## 1 0.23 -0.002  1.032  0.977     -1658.263

``` r
SE <- sandwichvar_se(Y, estimation$w1, estimation$loc, estimation$scale1, estimation$scale2)
print(round(SE,3))
```

    ##   w1_se loc_se scale1_se scale2_se
    ## 1 0.065  0.066     0.116     0.035

## Bayesian Inference

This project includes an implementation of Markov Chain Monte Carlo
(MCMC) for Bayesian inference of the FG distribution. We employ the
Metropolis-Hastings within Gibbs algorithm to draw samples from the
posterior distribution. Below is an example demonstrating how to conduct
Bayesian inference for the FG distribution using our provided functions.

``` r
source("./Bayes_MH_v4.0.r")

# Set initial parameter values
w1_initial <- 0.3
loc_initial <- 0
scale1_initial <- 1
scale2_initial <- 1

# Tune the proposal distribution parameters
loc_tune <- 0.1
scale1_tune <- 0.2
scale2_tune <- 0.1

# Run MCMC sampling
output <- MH(data = Y,
             w1_initial = w1_initial,
             loc_initial = loc_initial,
             scale1_initial = scale1_initial,
             scale2_initial = scale2_initial,
             loc_tune = loc_tune,
             scale1_tune = scale1_tune,
             scale2_tune = scale2_tune,
             random_seed = 100,
             burn_in = 5000,
             iteration = 10000,
             thin = 1)

estimation <- data.frame(t(round(summary(output)[[2]][, 3], 3)))
print(estimation)
```

    ##     w1    loc scale1 scale2
    ## 1 0.23 -0.003  1.036  0.981

``` r
SE <- round(summary(output)[[1]][, 2], 3)
print(SE)
```

    ##     w1    loc scale1 scale2 
    ##  0.063  0.065  0.135  0.037

## R shiny app for the Flexible Gumbel Distribution

We have deployed an [R Shiny
app](https://qingyang.shinyapps.io/gumbel_mixture/) for the FG
distribution. This interactive tool allows users to explore various
aspects of the FG distribution, including:

- Generating PDF plots and cumulative density plots for the FG
  distribution with different location/mode, scale, and weight
  parameters.
- Performing calculations and visualizations to verify the
  identifiability of the FG distribution.
- Calculating third and fourth central moments (skewness and kurtosis)
  of the FG distribution.
- Conducting both Frequentist and Bayesian inferences using uploaded CSV
  files.

The [R Shiny app](https://qingyang.shinyapps.io/gumbel_mixture/)
provides a user-friendly interface for conducting analyses and
visualizations related to the FG distribution. Feel free to explore and
utilize the app for your statistical needs.

## `STAN` and `JAGS` Programs

The repository includes `STAN` and `JAGS` programs for implementing the
Flexible Gumbel (FG) distribution:

- The `STAN` program can be accessed
  [here](https://github.com/rh8liuqy/flexible_Gumbel/tree/main/table2/model_fitting.stan).
- The `JAGS` program can be accessed
  [here](https://github.com/rh8liuqy/flexible_Gumbel/tree/main/table2/Gumbel_Mixture_JAGS.BUG).

These programs are useful for conducting Bayesian inference and fitting
the FG distribution to data using the `STAN` and `JAGS` software
platforms, respectively. Feel free to utilize these programs for your
statistical analyses.