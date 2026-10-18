# VOI-Gaussian-Approximation 2.0.0

Repository reorganized for reproducibility. The original files published with the paper are kept unchanged in `MDM_2017_Appendix/`.

## New

- `R/GA_functions.R`: maintained version of the GA functions, with documentation.
- `analysis/`: scripts that reproduce Box 1 and Appendices B and C, runnable from the repository root, with random seeds for the Monte Carlo examples.
- `tests/test_GA_functions.R`: checks that the maintained functions reproduce the published results, that EVSI tends to EVPPI as the sample size grows and to zero when it is zero, and that the data reproduce Table D.2 of the supplement.
- `results/expected_results.csv`: reference EVPI, EVPPI, and EVSI values.
- `CITATION.cff`, `data/README.md` (data dictionary).

## Fixes

- **Normal-Weibull MCMC example (Appendix C).** The JAGS model used `lambda[i] <- 1/pow(theta[i], v)`, but `v` was never defined, so the model failed to compile. It now uses the Weibull shape (`shape`). The same error appears in the code printed in the online supplement.
- **`predict.ga()` with several parameters.** `n` and `n0` were matched to parameters by position, and the single-parameter smooths and the tensor-product margins could use different positions when the terms in the formula were not in the same order. `n` and `n0` can now be named vectors matched by parameter name. Unnamed vectors keep the original behavior (order in which parameters first appear in the formula). Results for the models in the paper are identical.
- Tensor-product smooths are now detected by class (`te()`, `ti()`), and unsupported models (parametric terms other than the intercept, `t2()` smooths) give an informative error.

## Removed

- `README.html` (outdated rendering of the README) and an empty `Icon` file.

# VOI-Gaussian-Approximation 1.0.0

- Code and data accompanying Jalal H, Alarid-Escudero F. A Gaussian Approximation Approach for Value of Information Analysis. *Med Decis Making*. 2018;38(2):174-188. Archived on Zenodo: [doi:10.5281/zenodo.3263876](https://doi.org/10.5281/zenodo.3263876).
