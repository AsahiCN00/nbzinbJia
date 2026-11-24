#' Methods for zinbreg objects
#'
#' Standard S3 methods for fitted zero-inflated negative binomial models
#' returned by [zinbreg()]. Coefficient names follow `pscl::zeroinfl()`:
#' `count_*` for the count model and `zero_*` for the zero-inflation model.
#'
#' @param x,object A fitted `"zinbreg"` object (or, for
#'   `print.summary.zinbreg()`, a `"summary.zinbreg"` object).
#' @param digits Number of significant digits to print.
#' @param model For `coef()` and `vcov()`: `"full"` (both parts), `"count"`,
#'   or `"zero"`.
#' @param type For `residuals()`: `"response"` (`y - (1 - pi) mu`) or
#'   `"pearson"` (response residual divided by
#'   `sqrt((1 - pi) mu (1 + pi mu + mu / theta))`). For `predict()`:
#'   `"response"` (`(1 - pi) mu`, the mean of `Y`), `"count"` (`mu`, the
#'   mean of the NB component), `"zero"` (`pi`, the probability of a
#'   structural zero), or `"prob"` (a matrix of `P(Y = k)` for `k` in `at`).
#' @param newdata Optional data frame for prediction. If omitted, the
#'   original data are used.
#' @param at For `type = "prob"`: the counts at which to evaluate
#'   probabilities (default `0:max(y)`).
#' @param ... Further arguments (ignored).
#'
#' @return `summary()` returns an object of class `"summary.zinbreg"` with
#'   coefficient tables (estimates, standard errors, z values, p-values) for
#'   the count part (including `Log(theta)`) and the zero part.
#'   `logLik()` returns a `"logLik"` object with a `df` attribute (all
#'   coefficients plus one for theta), so `AIC()` and `BIC()` work.
#' @examplesIf requireNamespace("pscl", quietly = TRUE)
#' data("bioChemists", package = "pscl")
#' fit <- zinbreg(art ~ fem + kid5 + ment | ment, data = bioChemists)
#' coef(fit, model = "zero")
#' sqrt(diag(vcov(fit)))
#' AIC(fit)
#' nd <- data.frame(fem = "Women", kid5 = 0, ment = c(0, 10, 30))
#' predict(fit, newdata = nd, type = "zero")
#' predict(fit, newdata = nd, type = "response")
#' # expected number of zeros
#' sum(predict(fit, type = "prob", at = 0))
#' @name zinbreg-methods
NULL

#' @rdname zinbreg-methods
#' @export
print.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Zero-inflated negative binomial regression (EM with IRLS M-steps)\n\n")
  cat("Count model coefficients (negbin with log link):\n")
  print.default(format(x$coefficients$count, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  print.default(format(x$coefficients$zero, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\nTheta:", format(x$theta, digits = digits),
      " Log-likelihood:", format(x$loglik, digits = digits + 2L), "on",
      attr(stats::logLik(x), "df"), "Df\n")
  if (!x$converged) cat("Warning: EM did not converge\n")
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
summary.zinbreg <- function(object, ...) {
  kx <- length(object$coefficients$count)
  kz <- length(object$coefficients$zero)
  se <- sqrt(diag(object$vcov))
  make_table <- function(est, se) {
    z <- est / se
    cbind(Estimate = est, `Std. Error` = se, `z value` = z,
          `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  }
  count_tab <- make_table(c(object$coefficients$count, log(object$theta)),
                          se[c(seq_len(kx), kx + kz + 1L)])
  rownames(count_tab) <- c(names(object$coefficients$count), "Log(theta)")
  zero_tab <- make_table(object$coefficients$zero, se[kx + seq_len(kz)])
  rownames(zero_tab) <- names(object$coefficients$zero)
  structure(list(call = object$call,
                 coefficients = list(count = count_tab, zero = zero_tab),
                 theta = object$theta,
                 SE.theta = object$SE.theta,
                 loglik = object$loglik,
                 df = attr(stats::logLik(object), "df"),
                 aic = stats::AIC(object),
                 iter = object$iter,
                 converged = object$converged,
                 pearson = stats::residuals(object, type = "pearson")),
            class = "summary.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
print.summary.zinbreg <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Pearson residuals:\n")
  q <- stats::quantile(x$pearson)
  names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(q, digits = digits)
  cat("\nCount model coefficients (negbin with log link):\n")
  stats::printCoefmat(x$coefficients$count, digits = digits)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  stats::printCoefmat(x$coefficients$zero, digits = digits)
  cat("\nTheta =", format(x$theta, digits = digits),
      " (Std. Err.:", format(x$SE.theta, digits = digits), ")\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits + 2L), "on",
      x$df, "Df,  AIC:", format(x$aic, digits = digits + 2L), "\n")
  cat("EM iterations:", x$iter,
      if (x$converged) "(converged)" else "(NOT converged)", "\n")
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
coef.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
         count = cf$count,
         zero = cf$zero,
         full = stats::setNames(c(cf$count, cf$zero),
                                c(paste0("count_", names(cf$count)),
                                  paste0("zero_", names(cf$zero)))))
}

#' @rdname zinbreg-methods
#' @export
vcov.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  kx <- length(object$coefficients$count)
  kz <- length(object$coefficients$zero)
  idx <- switch(model,
                count = seq_len(kx),
                zero = kx + seq_len(kz),
                full = seq_len(kx + kz))
  V <- object$vcov[idx, idx, drop = FALSE]
  if (model != "full") {
    nm <- names(object$coefficients[[model]])
    dimnames(V) <- list(nm, nm)
  }
  V
}

#' @rdname zinbreg-methods
#' @export
logLik.zinbreg <- function(object, ...) {
  structure(object$loglik,
            df = length(object$coefficients$count) +
              length(object$coefficients$zero) + 1L,
            nobs = stats::nobs(object),
            class = "logLik")
}

#' @rdname zinbreg-methods
#' @export
nobs.zinbreg <- function(object, ...) sum(object$prior.weights != 0)

#' @rdname zinbreg-methods
#' @export
fitted.zinbreg <- function(object, ...) object$fitted.values

#' @rdname zinbreg-methods
#' @export
residuals.zinbreg <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  r <- object$y - object$fitted.values
  if (type == "pearson") {
    # Var(Y) = (1 - pi) mu [1 + mu (pi + 1 / theta)]
    mu <- object$mu
    p <- object$pi
    v <- (1 - p) * mu * (1 + mu * (p + 1 / object$theta))
    r <- r * sqrt(object$prior.weights) / sqrt(v)
  }
  r
}

#' @rdname zinbreg-methods
#' @export
predict.zinbreg <- function(object, newdata,
                            type = c("response", "count", "zero", "prob"),
                            at = NULL, ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    mu <- object$mu
    p <- object$pi
  } else {
    dx <- new_design(object$terms$count, newdata, object$xlevels$count,
                     object$contrasts$count)
    offx <- stats::model.offset(dx$mf)
    if (is.null(offx)) offx <- 0
    if (!is.null(object$call$offset)) {
      offx <- offx + eval(object$call$offset, newdata,
                          environment(object$terms$count))
    }
    dz <- new_design(object$terms$zero, newdata, object$xlevels$zero,
                     object$contrasts$zero)
    offz <- stats::model.offset(dz$mf)
    if (is.null(offz)) offz <- 0
    mu <- exp(drop(dx$X %*% object$coefficients$count) + offx)
    p <- stats::plogis(drop(dz$X %*% object$coefficients$zero) + offz)
  }
  switch(type,
         response = (1 - p) * mu,
         count = mu,
         zero = p,
         prob = {
           if (is.null(at)) at <- 0:max(object$y)
           # P(Y = k) = pi 1(k = 0) + (1 - pi) f_NB(k; mu, theta)
           pr <- vapply(at, function(k) {
             (k == 0) * p + (1 - p) *
               stats::dnbinom(k, size = object$theta, mu = mu)
           }, numeric(length(mu)))
           pr <- matrix(pr, nrow = length(mu))
           colnames(pr) <- at
           pr
         })
}
