#**************************************************************************
#* Reproduce all results. Run from the repository root:
#*   Rscript analysis/00_run_all.R
#* Set run_mcmc <- TRUE to also run the MCMC examples (needs JAGS, ~3 min).
#**************************************************************************

run_mcmc <- FALSE

scripts <- c("analysis/01_box1_evsi_single_parameter.R",
             "analysis/02_appendixB_evppi_evsi.R",
             "analysis/03_appendixC_n0_summary_statistic.R")
if (run_mcmc) scripts <- c(scripts, "analysis/04_appendixC_n0_mcmc.R")

for (s in scripts) {
  cat("\n==========", s, "==========\n")
  source(s, echo = FALSE, print.eval = TRUE)
}

cat("\n========== tests ==========\n")
source("tests/test_GA_functions.R")

cat("\n========== session info ==========\n")
print(sessionInfo())
