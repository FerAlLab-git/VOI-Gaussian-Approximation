# Data

`psa_demo_GA.csv` holds 10,000 probabilistic sensitivity analysis (PSA) samples from the Markov model of a hypothetical genetic disorder ("syndrome X") described in Appendix D of the paper's supplement. Three strategies are compared: A, B, and C (standard care). Willingness to pay is $50,000 per QALY; costs and benefits are discounted at 5% per year.

| Column | Description | Distribution |
|---|---|---|
| `MeanVisitsA` | Mean number of hospital visits, intervention A | Gamma(shape 10, scale 0.1), mean 1.0 |
| `MeanVisitsB` | Mean number of hospital visits, intervention B | Gamma(shape 20, scale 0.1), mean 2.0 |
| `ProbFailingA` | Annual probability of failing A | Beta(2, 8), mean 0.2 |
| `ProbFailingB` | Annual probability of failing B | Beta(3, 7), mean 0.3 |
| `nmbA` | Net monetary benefit of A ($) | Model output |
| `nmbB` | Net monetary benefit of B ($) | Model output |
| `nmmC` | Net monetary benefit of C ($). The column name is a typo for `nmbC`, kept so the file matches the published version. C has no uncertain parameters, so its NMB is constant. | Model output |

The effective prior sample size is *n0* = 10 for all four parameters (Appendix D).

Mean NMB: A $445,179; B $449,543; C $445,305. B is optimal with current information and is the best strategy in 4,586 of 10,000 samples (Table D.2).

This file is identical to `MDM_2017_Appendix/psa_demo_GA.csv`.
