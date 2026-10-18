#**************************************************************************
#* Appendix C: Estimating the effective prior sample size (n0) with the
#* summary-statistic method:  n0_hat = n * (Var(S) / Var(theta) - 1)
#*
#* Jalal H, Alarid-Escudero F. A Gaussian Approximation Approach for Value
#* of Information Analysis. Med Decis Making. 2018;38(2):174-188.
#*
#* Original versions: MDM_2017_Appendix/Appendix_EVSI_GA_*.R
#* Results are Monte Carlo estimates; a seed is set for reproducibility, and
#* values vary around the paper's reported estimates with other seeds.
#**************************************************************************

set.seed(2017)

#### Beta-binomial ####
## Prior: theta ~ beta(12, 8)   -> true n0 = a + b = 20
## Data:  x ~ binomial(theta, n = 50), summary statistic S = x / n
a <- 12; b <- 8
nSample <- 1000
n <- 50
theta <- rbeta(nSample, a, b)
S <- rbinom(nSample, n, theta) / n
n0.hat.beta.binomial <- n * (var(S) / var(theta) - 1)
n0.hat.beta.binomial   # paper: 21.93

#### Gamma-exponential ####
## Prior: lambda ~ gamma(shape = 20, scale = 0.1)   -> true n0 = a = 20
## Data:  x ~ exponential(lambda), summary statistic S = 1 / mean(x)
a <- 20; b <- 0.1
nSample <- 10000
n <- 100
theta <- rgamma(nSample, shape = a, scale = b)
S <- vapply(theta, function(l) 1 / mean(rexp(n, l)), numeric(1))
n0.hat.gamma.exponential <- n * (var(S) / var(theta) - 1)
n0.hat.gamma.exponential   # paper: 25

#### Gamma-Poisson ####
## Prior: lambda ~ gamma(shape = 50, scale = 0.01)   -> true n0 = 1 / b = 100
## Data:  x ~ Poisson(lambda), summary statistic S = mean(x)
a <- 50; b <- 0.01
nSample <- 10000
n <- 50
theta <- rgamma(nSample, shape = a, scale = b)
S <- vapply(theta, function(l) mean(rpois(n, l)), numeric(1))
n0.hat.gamma.poisson <- n * (var(S) / var(theta) - 1)
n0.hat.gamma.poisson   # paper: 104
