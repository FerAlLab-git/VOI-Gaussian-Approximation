#**************************************************************************
#*
#* Original work Copyright (C) 2017  Fernando Alarid-Escudero
#* Original work Copyright (C) 2017  Hawre Jalal
#*
#* Released under the MIT License. See the LICENSE file in the root of this
#* repository or <https://opensource.org/licenses/MIT>.
#*
#* Functions to compute the expected value of sample information (EVSI) with
#* the Gaussian approximation (GA) of:
#*   Jalal H, Alarid-Escudero F. A Gaussian Approximation Approach for Value
#*   of Information Analysis. Med Decis Making. 2018;38(2):174-188.
#*   https://doi.org/10.1177/0272989X17715627
#*
#* This file is the maintained version of MDM_2017_Appendix/GA_functions.R.
#* For the models used in the paper it returns the same results as the
#* original file (checked in tests/test_GA_functions.R). Changes (v2.0.0):
#*   - n and n0 are matched to parameters by name, so results no longer
#*     depend on the order of terms in the GAM formula. Unnamed vectors
#*     follow the order in which parameters first appear in the formula,
#*     as in the original code.
#*   - Tensor smooths (te, ti) are detected by class instead of by the
#*     first letter of their label.
#*   - Informative errors for unsupported models (parametric terms other
#*     than the intercept, t2 smooths).
#**************************************************************************

#' Expected conditional loss under the Gaussian approximation
#'
#' Computes the preposterior expectation of the opportunity loss for one
#' strategy by shrinking each basis function of a fitted GAM metamodel toward
#' its mean by the factor sqrt(n / (n + n0)).
#'
#' @param object A `gam` object (package mgcv) regressing the opportunity loss
#'   of one strategy on the parameter(s) of interest, using only smooth terms
#'   (`s()`, `te()`, `ti()`) plus an intercept.
#' @param n Sample size of the proposed study. A scalar (same n for all
#'   parameters) or a vector with one value per parameter. Named vectors are
#'   matched by parameter name.
#' @param n0 Effective prior sample size. Same format as `n`.
#' @param verbose If `TRUE`, prints the factor sqrt(n / (n + n0)) for each
#'   parameter.
#' @return A one-column matrix with the expected conditional loss for each PSA
#'   sample.
predict.ga <- function(object, n, n0, verbose = TRUE) {
  ### Parameters used in the metamodel (all model frame columns but the response)
  data <- object$model[, -1, drop = FALSE]
  param.names <- colnames(data)
  n.params <- length(param.names)

  ### Only an intercept plus smooth terms are supported
  if (object$nsdf != 1) {
    stop("predict.ga() supports GAMs with an intercept and smooth terms only ",
         "(s(), te(), ti()). Remove parametric terms from the formula.")
  }

  ### Match n and n0 to parameters
  n  <- match_to_params(n,  param.names, "n")
  n0 <- match_to_params(n0, param.names, "n0")

  ### Factor that shrinks each basis function toward its mean, sqrt(n/(n+n0))
  v.ga <- sqrt(n / (n + n0))
  names(v.ga) <- param.names
  if (verbose) {
    print(paste("Variance reduction factor =", round(v.ga, 3)))
  }

  ### Initialize matrix for preposterior of all basis functions
  n.smooth <- length(object$smooth)
  X <- matrix(NA, nrow(object$model), length(object$coefficients))
  X[, 1] <- 1

  for (k in seq_len(n.smooth)) {
    sm <- object$smooth[[k]]
    if (inherits(sm, "t2.smooth")) {
      stop("t2() smooths are not supported; use te() or ti() instead.")
    }
    if (inherits(sm, "tensor.smooth")) {
      ## One factor per margin, matched by the margin's parameter name
      v.margin <- vapply(sm$margin, function(mg) v.ga[[mg$term[1]]], numeric(1))
      Xfrag <- Predict.matrix.tensor.smooth.ga(sm, data, v.margin)
    } else {
      Xfrag <- Predict.smooth.ga(sm, data, v.ga[[sm$term[1]]])
    }
    X[, sm$first.para:sm$last.para] <- Xfrag
  }

  ### Expected conditional loss
  Beta <- coef(object)
  Ltilde <- X %*% Beta

  return(Ltilde)
}

match_to_params <- function(x, param.names, arg) {
  #### Expands a scalar or reorders a named vector so that it has one value
  #### per parameter, in the order of param.names
  n.params <- length(param.names)
  if (!is.null(names(x))) {
    missing <- setdiff(param.names, names(x))
    if (length(missing) > 0) {
      stop("'", arg, "' is named but has no value for: ",
           paste(missing, collapse = ", "))
    }
    return(unname(x[param.names]))
  }
  if (length(x) == 1) {
    return(rep(x, n.params))
  }
  if (length(x) != n.params) {
    stop("'", arg, "' should be either a scalar or a vector with one value ",
         "per parameter (", n.params, ": ", paste(param.names, collapse = ", "),
         ")")
  }
  return(x)
}

Predict.smooth.ga <- function(object, data, v.ga = 1) {
  #### Preposterior of each basis function of a single-parameter smooth

  ### Basis functions evaluated at the PSA samples
  X <- PredictMat(object, data)
  n.obs <- nrow(X)

  ### Shrink each basis function toward its mean
  ones <- matrix(1, n.obs, 1)
  X.ga <- v.ga * X + (1 - v.ga) * (ones %*% colMeans(X))

  return(X.ga)
}

Predict.matrix.tensor.smooth.ga <- function(object,
                                            data,
                                            v.ga = rep(1, length(object$margin))) {
  #### Preposterior of the basis functions of a tensor-product smooth.
  #### v.ga has one value per margin, in margin order.
  #### (Heavily based on Predict.matrix.tensor.smooth from package mgcv)

  m <- length(object$margin)
  X <- list()
  for (i in 1:m) {
    term <- object$margin[[i]]$term
    dat <- list()
    for (j in seq_along(term)) {
      dat[[term[j]]] <- data[[term[j]]]
    }
    X[[i]] <- if (!is.null(object$mc[i]))
      PredictMat(object$margin[[i]], dat, n = length(dat[[1]]))
    else Predict.matrix(object$margin[[i]], dat)
    n.obs <- nrow(X[[i]])
  }
  mxp <- length(object$XP)
  if (mxp > 0)
    for (i in 1:mxp) if (!is.null(object$XP[[i]]))
      X[[i]] <- X[[i]] %*% object$XP[[i]]

  ### Shrink the basis functions of each margin toward their means
  ones <- matrix(1, n.obs, 1)
  X.ga <- list()
  for (i in 1:m) {
    X.ga[[i]] <- v.ga[i] * X[[i]] + (1 - v.ga[i]) * (ones %*% colMeans(X[[i]]))
  }

  ### Tensor product of the shrunken margins
  T.ga <- tensor.prod.model.matrix(X.ga)

  return(T.ga)
}
