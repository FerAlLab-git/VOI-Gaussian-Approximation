#**************************************************************************
#* Box 1: EVSI for a single parameter with the Gaussian approximation (GA)
#*
#* Jalal H, Alarid-Escudero F. A Gaussian Approximation Approach for Value
#* of Information Analysis. Med Decis Making. 2018;38(2):174-188.
#*
#* Run from the repository root (open VOI-Gaussian-Approximation.Rproj).
#* Original version: MDM_2017_Appendix/EVSI_GA_Box1.R
#**************************************************************************

#### Step 1: Load packages, functions, and the PSA ####
library(mgcv)
library(matrixStats)
source("R/GA_functions.R")

psa <- read.csv("data/psa_demo_GA.csv")
head(psa)

## Net monetary benefit (NMB) of each strategy (A, B, C)
nmb <- psa[, 5:7]
## Parameter of interest: mean number of hospital visits under A (column 1)
theta_I <- psa[, 1]

#### Step 2: Optimal strategy (d*) with current information ####
d.star <- which.max(colMeans(nmb))

#### Step 3: Opportunity loss of each strategy relative to d* ####
loss <- nmb - nmb[, d.star]

#### Step 4: Metamodel of each strategy's loss on the parameter of interest ####
lmm1 <- gam(loss[, 1] ~ s(theta_I))
lmm2 <- gam(loss[, 2] ~ s(theta_I))
lmm3 <- gam(loss[, 3] ~ s(theta_I))

#### Step 5: EVPPI ####
Lhat  <- cbind(lmm1$fitted, lmm2$fitted, lmm3$fitted)
evppi <- mean(rowMaxs(Lhat))
evppi

#### Step 6: Prior (n0) and proposed study (n) sample sizes ####
n0 <- 10
n  <- 100

#### Step 7: Expected conditional loss from the preposterior (GA) ####
Ltilde1 <- predict.ga(lmm1, n = n, n0 = n0)
Ltilde2 <- predict.ga(lmm2, n = n, n0 = n0)
Ltilde3 <- predict.ga(lmm3, n = n, n0 = n0)
loss.predicted <- cbind(Ltilde1, Ltilde2, Ltilde3)

#### Step 8: EVSI ####
evsi <- mean(rowMaxs(loss.predicted))
evsi
