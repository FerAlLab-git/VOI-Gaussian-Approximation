#**************************************************************************
#* Checks for R/GA_functions.R. Run from the repository root:
#*   Rscript tests/test_GA_functions.R
#* Stops with an error if any check fails.
#**************************************************************************

suppressMessages({
  library(mgcv)
  library(matrixStats)
})

## Load the original (published) and maintained functions in separate environments
orig <- new.env(); sys.source("MDM_2017_Appendix/GA_functions.R", envir = orig)
new  <- new.env(); sys.source("R/GA_functions.R", envir = new)

psa   <- read.csv("data/psa_demo_GA.csv")
nmb   <- psa[, 5:7]
theta <- psa[, 1:4]
loss  <- nmb - nmb[, which.max(colMeans(nmb))]
tol   <- 1e-6

check <- function(ok, msg) {
  if (!isTRUE(ok)) stop("FAILED: ", msg, call. = FALSE)
  cat("ok   ", msg, "\n")
}

ga  <- function(env, m, n, n0) as.vector(env$predict.ga(m, n = n, n0 = n0, verbose = FALSE))
fit <- function(f, j) gam(as.formula(paste0("loss[, ", j, "] ", f)), data = theta)
f.one <- "~ s(ProbFailingA)"
f.two <- "~ s(MeanVisitsA) + s(ProbFailingA) + ti(MeanVisitsA, ProbFailingA)"

## 1. Data reproduce Table D.2 of the supplement
check(all(round(colMeans(nmb)) == c(445179, 449543, 445305)),
      "mean NMB of A, B, C match Table D.2")
check(all(tabulate(max.col(as.matrix(nmb), ties.method = "first"), 3) == c(2847, 4586, 2567)),
      "times each strategy is optimal match Table D.2")

## 2. Maintained functions give the same results as the published ones
m1 <- fit(f.one, 1)
m2 <- fit(f.two, 1)
check(max(abs(ga(orig, m1, 100, 10) - ga(new, m1, 100, 10))) < tol,
      "single parameter: same as published code")
for (d in list(c(100, 100, 10, 10), c(100, 100, 10, 50), c(100, 1000, 10, 50))) {
  check(max(abs(ga(orig, m2, d[1:2], d[3:4]) - ga(new, m2, d[1:2], d[3:4]))) < tol,
        paste0("two parameters (n = ", d[1], ", ", d[2], "; n0 = ", d[3], ", ", d[4],
               "): same as published code"))
}

## 3. Limits: n -> Inf gives the EVPPI metamodel; n = 0 gives no information
check(max(abs(ga(new, m1, 1e12, 1) - fitted(m1))) < 1e-4,
      "n -> Inf reproduces fitted values (EVSI -> EVPPI), one parameter")
check(max(abs(ga(new, m2, 1e12, 1) - fitted(m2))) < 1e-4,
      "n -> Inf reproduces fitted values (EVSI -> EVPPI), two parameters")
check(sd(ga(new, m1, 0, 10)) < tol,
      "n = 0 gives a constant expected loss (EVSI = 0)")

## 4. Named n and n0 give the same answer whatever the order of terms in the formula
f.rev <- "~ s(ProbFailingA) + s(MeanVisitsA) + ti(MeanVisitsA, ProbFailingA)"
m.rev <- fit(f.rev, 1)
n  <- c(MeanVisitsA = 100, ProbFailingA = 1000)
n0 <- c(MeanVisitsA = 10,  ProbFailingA = 50)
check(max(abs(ga(new, m2, n, n0) - ga(new, m.rev, n, n0))) < 1e-4,
      "named n, n0 are invariant to formula order")

## 5. Results match results/expected_results.csv
expected <- read.csv("results/expected_results.csv")
lm1 <- lapply(1:3, function(j) fit(f.one, j))
lm2 <- lapply(1:3, function(j) fit(f.two, j))
evsi <- function(ms, n, n0) mean(rowMaxs(sapply(ms, function(m) ga(new, m, n, n0))))
got <- c(evsi(lm1, 100, 10),
         evsi(lm2, c(100, 100), c(10, 10)),
         evsi(lm2, c(100, 100), c(10, 50)),
         evsi(lm2, c(100, 1000), c(10, 50)))
want <- expected$value[expected$source %in% c("Appendix B", "Appendix B case 1",
                                              "Appendix B case 2", "Appendix B case 3") &
                         expected$quantity == "EVSI"]
check(max(abs(got - want)) < 1e-2, "EVSI values match results/expected_results.csv")

cat("\nAll checks passed (mgcv ", as.character(packageVersion("mgcv")), ").\n", sep = "")
