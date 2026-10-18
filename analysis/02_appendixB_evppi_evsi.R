#**************************************************************************
#* Appendix B: EVPI, EVPPI, and EVSI with the Gaussian approximation (GA)
#* for one and two parameters, with balanced and unbalanced designs.
#*
#* Jalal H, Alarid-Escudero F. A Gaussian Approximation Approach for Value
#* of Information Analysis. Med Decis Making. 2018;38(2):174-188.
#*
#* Run from the repository root (open VOI-Gaussian-Approximation.Rproj).
#* Original version: MDM_2017_Appendix/EVSI_GA_Appendix.R
#**************************************************************************

### Load packages and functions
library(matrixStats)  # row and column statistics
library(mgcv)         # spline (GAM) metamodels
source("R/GA_functions.R")

### Load PSA
psa <- read.csv("data/psa_demo_GA.csv")
head(psa)

## Net monetary benefit (NMB) matrix (strategies A, B, C)
nmb <- psa[, 5:7]
## Model parameter inputs
theta <- psa[, 1:4]
## Number of simulations and strategies
n.sim        <- nrow(nmb)
n.strategies <- ncol(nmb)

### Optimal strategy (d*) based on the highest expected NMB
d.star <- which.max(colMeans(nmb))
d.star

### Opportunity loss matrix
loss <- nmb - nmb[, d.star]

### EVPI
evpi <- mean(rowMaxs(as.matrix(loss)))
evpi

#========================#
#### Single parameter ####
#========================#
## Parameter of interest: probability of failing A (column 3)
sel.param <- 3
f1 <- paste("~ s(", colnames(theta)[sel.param], ")")
lmm1 <- gam(as.formula(paste("loss[, 1]", f1)), data = theta)
lmm2 <- gam(as.formula(paste("loss[, 2]", f1)), data = theta)
lmm3 <- gam(as.formula(paste("loss[, 3]", f1)), data = theta)

#### EVPPI on one parameter ####
loss.hat <- cbind(lmm1$fitted, lmm2$fitted, lmm3$fitted)
evppi <- mean(rowMaxs(loss.hat))
evppi

#### EVSI on one parameter ####
n0 <- 10   # effective prior sample size
n  <- 100  # sample size of the proposed study

loss.tilde <- cbind(predict.ga(lmm1, n = n, n0 = n0),
                    predict.ga(lmm2, n = n, n0 = n0),
                    predict.ga(lmm3, n = n, n0 = n0))
evsi <- mean(rowMaxs(loss.tilde))
evsi

#======================#
#### Two parameters ####
#======================#
## Parameters of interest: mean visits under A (1) and probability of failing A (3)
sel.params <- c(1, 3)
p <- colnames(theta)[sel.params]
f2 <- paste0("~ s(", p[1], ") + s(", p[2], ") + ti(", p[1], ", ", p[2], ")")
lmm1 <- gam(as.formula(paste("loss[, 1]", f2)), data = theta)
lmm2 <- gam(as.formula(paste("loss[, 2]", f2)), data = theta)
lmm3 <- gam(as.formula(paste("loss[, 3]", f2)), data = theta)

#### EVPPI on two parameters ####
loss.hat <- cbind(lmm1$fitted, lmm2$fitted, lmm3$fitted)
evppi <- mean(rowMaxs(loss.hat))
evppi

## EVSI for a given design (n and n0 can be named by parameter)
evsi_ga <- function(n, n0) {
  loss.tilde <- cbind(predict.ga(lmm1, n = n, n0 = n0),
                      predict.ga(lmm2, n = n, n0 = n0),
                      predict.ga(lmm3, n = n, n0 = n0))
  mean(rowMaxs(loss.tilde))
}

#### Case 1: same n0 and n ####
evsi_ga(n = c(MeanVisitsA = 100, ProbFailingA = 100),
        n0 = c(MeanVisitsA = 10,  ProbFailingA = 10))

#### Case 2: different n0, same n ####
evsi_ga(n = c(MeanVisitsA = 100, ProbFailingA = 100),
        n0 = c(MeanVisitsA = 10,  ProbFailingA = 50))

#### Case 3: different n0 and n ####
evsi_ga(n = c(MeanVisitsA = 100, ProbFailingA = 1000),
        n0 = c(MeanVisitsA = 10,  ProbFailingA = 50))
