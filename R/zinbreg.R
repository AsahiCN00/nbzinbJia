#' Control parameters for zinbreg()
#'
#' @param epsilon Convergence tolerance for the EM algorithm: iterations stop
#'   when the observed-data log-likelihood changes by less than `epsilon`.
#'   EM converges linearly (slowly), so a small change in the log-likelihood
#'   can still leave the parameters noticeably away from the MLE; the default
#'   `1e-10` gives coefficients accurate to about `1e-5`.
#' @param maxit Maximum number of EM iterations.
#' @param inner_maxit Maximum number of iterations of each IRLS /
#'   Newton-Raphson loop inside the M-steps.
#' @param trace Logical; if `TRUE`, print the observed-data log-likelihood at
#'   each EM iteration.
#' @param theta_max Upper bound for theta (see [nbreg_control()]).
#' @param theta_warn A warning is issued when the estimated theta exceeds this
#'   value.
#'
#' @return A list of control parameters.
#' @examples
#' zinbreg_control(epsilon = 1e-10, trace = TRUE)
#' @export
zinbreg_control <- function(epsilon = 1e-10, maxit = 5000, inner_maxit = 100,
                            trace = FALSE, theta_max = 1e8, theta_warn = 1e5) {
  if (!is.numeric(epsilon) || length(epsilon) != 1 || epsilon <= 0) {
    stop("'epsilon' must be a single positive number", call. = FALSE)
  }
  if (!is.numeric(maxit) || length(maxit) != 1 || maxit < 1 ||
      !is.numeric(inner_maxit) || length(inner_maxit) != 1 ||
      inner_maxit < 1) {
    stop("'maxit' and 'inner_maxit' must be single numbers >= 1",
         call. = FALSE)
  }
  if (!is.numeric(theta_max) || theta_max <= 0 ||
      !is.numeric(theta_warn) || theta_warn <= 0) {
    stop("'theta_max' and 'theta_warn' must be positive", call. = FALSE)
  }
  list(epsilon = epsilon, maxit = as.integer(maxit),
       inner_maxit = as.integer(inner_maxit), trace = isTRUE(trace),
       theta_max = theta_max, theta_warn = theta_warn)
}

#' Zero-inflated negative binomial regression by EM with IRLS M-steps
#'
#' Fits the ZINB model: with probability `pi` an observation is a structural
#' zero, and otherwise it follows an NB2 distribution with mean `mu`, where
#' `log(mu) = X beta + offset` and `logit(pi) = Z gamma`. The observed-data
#' log-likelihood is maximized by the EM algorithm. The E-step computes the
#' posterior probability that each zero is structural; the M-step fits a
#' logistic regression with that fractional response (by IRLS) and a
#' weighted NB regression (IRLS for beta, Newton-Raphson for theta). No
#' general-purpose optimizer is used. Standard errors come from the analytic
#' observed information matrix of `(beta, gamma, log(theta))`. The results
#' reproduce `pscl::zeroinfl(..., dist = "negbin")`.
#'
#' @param formula A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model
#'   covariates before `|`, zero-model covariates after it. Without `|`, the
#'   same covariates are used in both parts.
#' @param data A data frame containing the variables in `formula`.
#' @param offset Optional offset for the count model, evaluated in `data`.
#'   Offsets can also be given with `offset()` in either part of the formula.
#' @param control A list from [zinbreg_control()].
#' @param weights Optional non-negative prior weights.
#' @param subset Optional expression selecting the rows of `data` to use.
#' @param na.action How to handle missing values.
#'
#' @return An object of class `"zinbreg"`: a list containing `coefficients`
#'   (a list with elements `count` and `zero`), `vcov` (covariance matrix of
#'   all coefficients and `log(theta)`), `theta`, `SE.theta`, `SE.logtheta`,
#'   `loglik`, `iter`, `converged`, `fitted.values` (`(1 - pi) mu`), `mu`,
#'   `pi`, `call`, and the model information needed by [predict.zinbreg()].
#' @seealso [zinbreg-methods], [nbreg()]
#' @examplesIf requireNamespace("pscl", quietly = TRUE)
#' data("bioChemists", package = "pscl")
#' fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
#'                data = bioChemists)
#' summary(fit)
#' head(predict(fit, type = "zero"))
#' @export
zinbreg <- function(formula, data, offset, control = zinbreg_control(),
                    weights, subset, na.action) {
  cl <- match.call()
  if (missing(data)) data <- environment(formula)

  # Parse the two-part formula. Without '|', copy the count covariates to
  # the zero part (as zeroinfl() does).
  f <- Formula::as.Formula(formula)
  if (length(f)[1L] != 1L) stop("the formula needs one response", call. = FALSE)
  if (length(f)[2L] > 2L) {
    stop("the formula can have at most two parts: y ~ x | z", call. = FALSE)
  }
  if (length(f)[2L] == 1L) {
    f <- Formula::as.Formula(stats::formula(f, rhs = 1L),
                             stats::formula(f, lhs = 0L, rhs = 1L))
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- f
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mtX <- stats::terms(f, data = data, rhs = 1L)
  mtZ <- stats::terms(f, data = data, lhs = 0L, rhs = 2L)
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)

  y <- stats::model.response(mf, "numeric")
  if (is.null(y)) stop("the formula must have a response", call. = FALSE)
  check_counts(y)
  if (all(y > 0)) {
    stop("the response has no zeros; a zero-inflated model is not ",
         "appropriate", call. = FALSE)
  }
  n <- NROW(y)
  a <- as.vector(stats::model.weights(mf))
  if (is.null(a)) a <- rep(1, n)
  if (any(a < 0)) stop("weights must be non-negative", call. = FALSE)
  offx <- terms_offset(mf, mtX, include_arg = TRUE)
  offz <- terms_offset(mf, mtZ, include_arg = FALSE)

  fit <- zinb_fit(X, Z, y, offx, offz, a, control)

  fit$call <- cl
  fit$formula <- f
  fit$terms <- list(count = mtX, zero = mtZ)
  fit$model <- mf
  fit$xlevels <- list(count = stats::.getXlevels(mtX, mf),
                      zero = stats::.getXlevels(mtZ, mf))
  fit$contrasts <- list(count = attr(X, "contrasts"),
                        zero = attr(Z, "contrasts"))
  class(fit) <- "zinbreg"
  fit
}

# ---------------------------------------------------------------------------
# EM algorithm (Section 3.2)
# ---------------------------------------------------------------------------
zinb_fit <- function(X, Z, y, offx, offz, a, control) {
  kx <- ncol(X)
  kz <- ncol(Z)
  nb_ctrl <- nbreg_control(epsilon = control$epsilon * 1e-2,
                           maxit = control$inner_maxit,
                           theta_max = control$theta_max,
                           theta_warn = control$theta_warn)

  # Starting values (Section 3.2):
  #  - (beta, theta) from an NB fit to all observations;
  #  - gamma from an IRLS logistic regression of I(y_i = 0) on z_i.
  nb0 <- nb_fit(X, y, offx, a, nb_ctrl, warn = FALSE)
  beta <- nb0$coefficients
  theta <- nb0$theta
  gamma <- logit_irls(Z, as.numeric(y == 0), offz, a,
                      gamma = stats::setNames(rep(0, kz), colnames(Z)),
                      control = nb_ctrl)

  mu <- exp(drop(X %*% beta) + offx)
  eta_z <- drop(Z %*% gamma) + offz
  ll <- zinb_loglik(y, mu, theta, eta_z, a)
  ll_path <- ll
  converged <- FALSE
  decreased <- FALSE
  iter <- 0L
  for (iter in seq_len(control$maxit)) {
    ll_old <- ll

    # E-step: tau_i = P(D_i = 1 | y_i) at the current (beta, gamma, theta)
    tau <- zinb_estep(y, mu, theta, eta_z)

    # M-step, zero part: logistic IRLS with fractional response tau
    gamma <- logit_irls(Z, tau, offz, a, gamma, nb_ctrl)
    eta_z <- drop(Z %*% gamma) + offz

    # M-step, count part: the NB fit of Section 3.1 with prior weights
    # a_i (1 - tau_i), warm-started at the current (beta, theta)
    nb <- nb_fit(X, y, offx, a * (1 - tau), nb_ctrl, start = beta,
                 theta_start = theta, warn = FALSE)
    beta <- nb$coefficients
    theta <- nb$theta
    mu <- nb$fitted.values

    # Convergence check on the observed-data log-likelihood.
    # EM must never decrease it; a decrease beyond rounding error signals
    # a bug, so it is recorded and reported.
    ll <- zinb_loglik(y, mu, theta, eta_z, a)
    ll_path <- c(ll_path, ll)
    if (ll < ll_old - 1e-8 * max(1, abs(ll_old))) decreased <- TRUE
    if (control$trace) {
      cat(sprintf("zinbreg EM iter %4d: logLik = %.10f, theta = %.6f\n",
                  iter, ll, theta))
    }
    if (abs(ll - ll_old) < control$epsilon) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    warning("zinbreg: EM did not converge in ", control$maxit, " iterations",
            call. = FALSE)
  }
  if (decreased) {
    warning("zinbreg: the log-likelihood decreased during EM", call. = FALSE)
  }
  if (theta > control$theta_warn) {
    warning("zinbreg: theta is very large (", signif(theta, 4), "); the ",
            "count part shows little overdispersion (zero-inflated Poisson ",
            "limit)", call. = FALSE)
  }

  # Standard errors: invert the observed information -H of the
  # observed-data log-likelihood at the MLE, in (beta, gamma, log theta).
  sh <- zinb_score_hessian(X, Z, y, a, beta, gamma, theta, offx, offz)
  V <- tryCatch(solve(-sh$hessian), error = function(e) {
    warning("zinbreg: observed information is singular; standard errors ",
            "are not available", call. = FALSE)
    matrix(NA_real_, kx + kz + 1, kx + kz + 1)
  })
  nm <- c(paste0("count_", colnames(X)), paste0("zero_", colnames(Z)),
          "Log(theta)")
  dimnames(V) <- list(nm, nm)
  se_logtheta <- sqrt(V[kx + kz + 1, kx + kz + 1])

  pi_hat <- stats::plogis(eta_z)
  # Boundary case: if every pi_i is essentially 0, the MLE is at
  # gamma -> -Inf (no zero inflation). The log-likelihood is then flat in
  # gamma, so gamma and its SE are not identified.
  if (max(pi_hat) < 1e-6) {
    warning("zinbreg: the estimated zero-inflation probabilities are ",
            "essentially zero (max pi = ", signif(max(pi_hat), 3), "); the ",
            "zero-model coefficients are not identified and their standard ",
            "errors are unreliable. The NB model (nbreg) may suffice.",
            call. = FALSE)
  }
  list(coefficients = list(count = beta, zero = gamma),
       vcov = V,
       theta = theta,
       SE.logtheta = se_logtheta,
       # delta method: SE(theta) = theta * SE(log theta)
       SE.theta = theta * se_logtheta,
       loglik = ll,
       loglik_path = ll_path,
       iter = iter,
       converged = converged,
       gradient = stats::setNames(sh$gradient, nm),
       fitted.values = (1 - pi_hat) * mu,
       mu = mu,
       pi = pi_hat,
       y = y,
       prior.weights = a,
       offset = list(count = offx, zero = offz),
       df.residual = sum(a > 0) - (kx + kz + 1),
       control = control)
}
