# VOI-Gaussian-Approximation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3263875.svg)](https://doi.org/10.5281/zenodo.3263875)

R code for the Gaussian approximation (GA) method to compute the **expected value of sample information (EVSI)** from a single probabilistic sensitivity analysis (PSA), as described in:

> Jalal H, Alarid-Escudero F. A Gaussian Approximation Approach for Value of Information Analysis. *Medical Decision Making*. 2018;38(2):174-188. [doi:10.1177/0272989X17715627](https://doi.org/10.1177/0272989X17715627)

## The method in brief

EVSI tells you how much a proposed study of size *n* is worth before you run it. Standard approaches need nested Monte Carlo simulation, which is slow for complex models. The GA approach needs only the PSA you already have:

1. Regress each strategy's opportunity loss on the parameter(s) the study would inform, using a spline metamodel (a generalized additive model, GAM).
2. Approximate the preposterior of those parameters by shrinking each basis function of the metamodel toward its mean by a factor of sqrt(*n* / (*n* + *n0*)), where *n0* is the effective prior sample size.
3. EVSI is the expected maximum of the resulting expected conditional losses.

The same metamodel also gives the expected value of partial perfect information (EVPPI), which is the limit of EVSI as *n* grows.

## Use it with your own PSA

```r
library(mgcv)
library(matrixStats)
source("R/GA_functions.R")

# nmb: matrix of net monetary benefits (one column per strategy)
# theta: data frame of PSA parameter samples
loss <- nmb - nmb[, which.max(colMeans(nmb))]

# One metamodel per strategy on the parameter of interest
fits <- lapply(seq_len(ncol(loss)),
               function(j) gam(loss[, j] ~ s(my_param), data = theta))

# EVSI for a study of n = 100 when the prior is worth n0 = 10 observations
loss_tilde <- sapply(fits, function(f) predict.ga(f, n = 100, n0 = 10))
evsi <- mean(rowMaxs(loss_tilde))
```

For two or more parameters, use `s(a) + s(b) + ti(a, b)` and pass `n` and `n0` as named vectors, e.g. `n = c(a = 100, b = 1000)`. Appendix C of the paper (and `analysis/03_appendixC_n0_summary_statistic.R`) shows how to estimate *n0* when it is not available directly from a conjugate prior.

## Reproduce the paper

Open `VOI-Gaussian-Approximation.Rproj` (or set the working directory to the repository root) and run:

```r
source("analysis/00_run_all.R")
```

| Script | Reproduces |
|---|---|
| `analysis/01_box1_evsi_single_parameter.R` | Box 1: step-by-step EVSI for one parameter |
| `analysis/02_appendixB_evppi_evsi.R` | Appendix B: EVPI, EVPPI, and EVSI for one and two parameters, balanced and unbalanced designs |
| `analysis/03_appendixC_n0_summary_statistic.R` | Appendix C: estimating *n0* with a summary statistic (beta-binomial, gamma-exponential, gamma-Poisson) |
| `analysis/04_appendixC_n0_mcmc.R` | Appendix C: estimating *n0* with MCMC (beta-binomial, normal-Weibull); requires JAGS |
| `tests/test_GA_functions.R` | Checks that the code reproduces the published results and behaves correctly at the limits |

Expected EVPI, EVPPI, and EVSI values are in [`results/expected_results.csv`](results/expected_results.csv). The *n0* estimates in Appendix C are Monte Carlo estimates, so they vary slightly with the random seed.

## Repository structure

```
R/                   GA functions (maintained version)
analysis/            Scripts that reproduce Box 1 and Appendices B and C
data/                Example PSA dataset (see data/README.md)
results/             Expected results
tests/               Reproducibility and correctness checks
MDM_2017_Appendix/   Original files as published with the paper (unchanged)
renv.lock, renv/     Recorded package versions (see Reproducible environment)
```

## Requirements

- R (tested with R 4.3.3 and 4.5.2)
- [mgcv](https://cran.r-project.org/package=mgcv) and [matrixStats](https://cran.r-project.org/package=matrixStats)
- For the MCMC examples only: [JAGS](https://mcmc-jags.sourceforge.io) and [R2jags](https://cran.r-project.org/package=R2jags)

The paper used mgcv 1.8-17. All EVPI, EVPPI, and EVSI values are identical with mgcv 1.8-17, 1.9-1, and 1.9-3.

### Reproducible environment (renv)

The exact package versions used to check the results (R 4.5.2, mgcv 1.9-3, matrixStats 1.5.0, R2jags 0.8-9, JAGS 4.3.1) are recorded in [`renv.lock`](renv.lock) with [renv](https://rstudio.github.io/renv/). When you open the project, renv activates automatically. To install the recorded versions into a project-specific library, run:

```r
renv::restore()
```

renv does not install R or JAGS; install those separately. If you only want to use `R/GA_functions.R` in your own project, you do not need renv: copy the file and install mgcv and matrixStats.

## Errata and changes

See [NEWS.md](NEWS.md). In short:

- **Normal-Weibull MCMC example (Appendix C):** the published JAGS model referred to an undefined variable `v`, so it did not run. Corrected in `analysis/04_appendixC_n0_mcmc.R`.
- **`predict.ga()`:** `n` and `n0` can now be named by parameter, so results no longer depend on the order of terms in the GAM formula. Results for the models in the paper are unchanged.

## How to cite

If you use this code, please cite the paper above and this software (see [CITATION.cff](CITATION.cff), or use "Cite this repository" on GitHub).

## License

MIT License. See [LICENSE](LICENSE).

## Authors

[Fernando Alarid-Escudero](https://profiles.stanford.edu/fernando-alarid-escudero) (Stanford University) and [Hawre Jalal](https://github.com/hjalal). Questions and bug reports: please [open an issue](../../issues).

Part of [FerAl Lab](https://github.com/FerAlLab-git).
