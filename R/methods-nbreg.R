#' Methods for nbreg objects
#'
#' Standard S3 methods for fitted negative binomial regression models
#' returned by [nbreg()].
#'
#' @param x,object A fitted `"nbreg"` object (or, for
#'   `print.summary.nbreg()`, a `"summary.nbreg"` object).
#' @param digits Number of significant digits to print.
#' @param type For `residuals()`: `"response"` (`y - mu`) or `"pearson"`
#'   (`(y - mu) / sqrt(mu + mu^2 / theta)`). For `predict()`: `"link"`
#'   (linear predictor) or `"response"` (`mu`).
#' @param newdata Optional data frame for prediction. If omitted, the fitted
#'   values from the original data are returned.
#' @param ... Further arguments (ignored).
#'
#' @return `summary()` returns an object of class `"summary.nbreg"` with a
#'   coefficient table (estimates, standard errors, z values, p-values).
#'   `logLik()` returns a `"logLik"` object with a `df` attribute (number of
#'   coefficients plus one for theta), so `AIC()` and `BIC()` work.
#'   The other methods return the corresponding vectors or matrices.
#' @examplesIf requireNamespace("pscl", quietly = TRUE)
#' data("bioChemists", package = "pscl")
#' fit <- nbreg(art ~ fem + ment, data = bioChemists)
#' coef(fit)
#' sqrt(diag(vcov(fit)))
#' logLik(fit)
#' BIC(fit)
#' head(residuals(fit, type = "pearson"))
#' predict(fit, newdata = data.frame(fem = "Women", ment = c(0, 10)),
#'         type = "response")
#' @name nbreg-methods
NULL

#' @rdname nbreg-methods
#' @export
print.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Negative binomial regression (IRLS + Newton-Raphson for theta)\n\n")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\nTheta:", format(x$theta, digits = digits),
      " Log-likelihood:", format(x$loglik, digits = digits + 2L),
      "on", length(x$coefficients) + 1L, "Df\n")
  if (!x$converged) cat("Warning: algorithm did not converge\n")
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
summary.nbreg <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(diag(object$vcov))
  zval <- est / se
  coefs <- cbind(Estimate = est, `Std. Error` = se, `z value` = zval,
                 `Pr(>|z|)` = 2 * stats::pnorm(-abs(zval)))
  structure(list(call = object$call,
                 coefficients = coefs,
                 theta = object$theta,
                 SE.theta = object$SE.theta,
                 loglik = object$loglik,
                 df = length(est) + 1L,
                 aic = stats::AIC(object),
                 iter = object$iter,
                 converged = object$converged,
                 pearson = stats::residuals(object, type = "pearson")),
            class = "summary.nbreg")
}

#' @rdname nbreg-methods
#' @export
print.summary.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Pearson residuals:\n")
  q <- stats::quantile(x$pearson)
  names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(q, digits = digits)
  cat("\nCoefficients (negative binomial with log link):\n")
  stats::printCoefmat(x$coefficients, digits = digits)
  cat("\nTheta:", format(x$theta, digits = digits),
      " Std. Err.:", format(x$SE.theta, digits = digits), "\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits + 2L), "on",
      x$df, "Df,  AIC:", format(x$aic, digits = digits + 2L), "\n")
  cat("Outer iterations (IRLS / Newton-Raphson):", x$iter,
      if (x$converged) "(converged)" else "(NOT converged)", "\n")
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
coef.nbreg <- function(object, ...) object$coefficients

#' @rdname nbreg-methods
#' @export
vcov.nbreg <- function(object, ...) object$vcov

#' @rdname nbreg-methods
#' @export
logLik.nbreg <- function(object, ...) {
  structure(object$loglik,
            df = length(object$coefficients) + 1L,
            nobs = stats::nobs(object),
            class = "logLik")
}

#' @rdname nbreg-methods
#' @export
nobs.nbreg <- function(object, ...) sum(object$prior.weights != 0)

#' @rdname nbreg-methods
#' @export
fitted.nbreg <- function(object, ...) object$fitted.values

#' @rdname nbreg-methods
#' @export
residuals.nbreg <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  r <- object$y - object$fitted.values
  if (type == "pearson") {
    mu <- object$fitted.values
    r <- r * sqrt(object$prior.weights) / sqrt(mu + mu^2 / object$theta)
  }
  r
}

#' @rdname nbreg-methods
#' @export
predict.nbreg <- function(object, newdata, type = c("link", "response"),
                          ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    eta <- object$linear.predictors
  } else {
    X <- new_design(object$terms, newdata, object$xlevels, object$contrasts)
    off <- new_offset(object, newdata)
    eta <- drop(X$X %*% object$coefficients) + off
  }
  if (type == "response") exp(eta) else eta
}

# Design matrix for new data, using the factor levels and contrasts of the
# original fit.
new_design <- function(terms, newdata, xlevels, contrasts) {
  tt <- stats::delete.response(terms)
  mf <- stats::model.frame(tt, newdata, xlev = xlevels,
                           na.action = stats::na.pass)
  X <- stats::model.matrix(tt, mf, contrasts.arg = contrasts)
  list(X = X, mf = mf)
}

# Offset for new data: offset() terms in the formula plus any offset given
# through the 'offset' argument of the original call.
new_offset <- function(object, newdata, terms = object$terms) {
  tt <- stats::delete.response(terms)
  mf <- stats::model.frame(tt, newdata, xlev = object$xlevels,
                           na.action = stats::na.pass)
  off <- stats::model.offset(mf)
  if (is.null(off)) off <- rep(0, nrow(mf))
  if (!is.null(object$call$offset)) {
    off <- off + eval(object$call$offset, newdata,
                      environment(object$terms))
  }
  off
}
