#**************************************************************************
#* Appendix C: Estimating the effective prior sample size (n0) with MCMC:
#*   n0_hat = n * (Var(theta) / Var_X(phi) - 1)
#* where phi is the posterior mean of theta for each simulated dataset.
#*
#* Jalal H, Alarid-Escudero F. A Gaussian Approximation Approach for Value
#* of Information Analysis. Med Decis Making. 2018;38(2):174-188.
#*
#* Requires JAGS (https://mcmc-jags.sourceforge.io) and the R2jags package.
#* Takes about 3 minutes. Results are Monte Carlo estimates; a seed is set.
#*
#* Original versions: MDM_2017_Appendix/Appendix_EVSI_GA_*_mcmc.R
#* Erratum: in the published normal-Weibull code, the JAGS model used an
#* undefined variable `v` (lambda[i] <- 1/pow(theta[i], v)), so the model
#* failed to compile. It is corrected below to use the Weibull shape.
#**************************************************************************

library(R2jags)
set.seed(2017)

#### Beta-binomial ####
## Prior: p ~ beta(12, 8)   -> true n0 = a + b = 20
## Data:  x ~ binomial(p, n = 50)
a <- 12; b <- 8
n <- 50
nSample <- 2000
p <- rbeta(nSample, a, b)
x <- matrix(rbinom(nSample, n, p), nSample, 1)

bb.model <- function() {
  for (i in 1:nSample) {
    p[i] ~ dbeta(a, b)       # prior
    x[i, 1] ~ dbin(p[i], n)  # likelihood
  }
}
bb.sim <- jags(data = list(x = x, n = n, nSample = nSample, a = a, b = b),
               inits = NULL, model.file = bb.model,
               parameters.to.save = "p",
               n.chains = 1, n.iter = 11000, n.burnin = 1000, n.thin = 1,
               progress.bar = "none")
prepost <- colMeans(bb.sim$BUGSoutput$sims.list$p)  # posterior mean per dataset
n0.hat.bb.mcmc <- n * (var(p) / var(prepost) - 1)
n0.hat.bb.mcmc   # paper: 21.17

#### Normal-Weibull ####
## Prior: theta ~ normal(1, sd = 0.2)
## Data:  x ~ Weibull(shape = 1, scale = theta), i.e., JAGS dweib(shape, 1/theta^shape)
shape <- 1
n <- 100        # observations per simulated study
m <- 1000       # prior samples
q <- 2000       # posterior samples kept
n.burnin <- 200
theta <- rnorm(m, 1, 0.2)
x <- t(vapply(theta, function(th) rweibull(n, shape, th), numeric(n)))

nw.model <- function() {
  tau <- 1 / pow(0.2, 2)
  for (i in 1:m) {
    theta[i] ~ dnorm(1, tau)                  # prior
    lambda[i] <- 1 / pow(theta[i], shape)     # corrected: was pow(theta[i], v)
    for (j in 1:n) {
      x[i, j] ~ dweib(shape, lambda[i])       # likelihood
    }
  }
}
nw.sim <- jags(data = list(x = x, n = n, m = m, shape = shape),
               inits = NULL, model.file = nw.model,
               parameters.to.save = "theta",
               n.chains = 1, n.iter = q + n.burnin, n.burnin = n.burnin, n.thin = 1,
               progress.bar = "none")
prepost <- colMeans(nw.sim$BUGSoutput$sims.list$theta)
n0.hat.nw.mcmc <- n * (var(theta) / var(prepost) - 1)
n0.hat.nw.mcmc   # paper: about 25
