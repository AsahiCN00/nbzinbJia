#' Control parameters for nbreg()
#'
#' @param epsilon Convergence tolerance: the outer loop stops when the
#'   log-likelihood changes by less than `epsilon` between iterations.
#' @param maxit Maximum number of outer iterations (IRLS for beta alternated
#'   with Newton-Raphson for theta). Also caps the inner IRLS and
#'   Newton-Raphson loops.
#' @param trace Logical; if `TRUE`, print the log-likelihood and theta at each
#'   outer iteration.
#' @param theta_max Upper bound for theta. Estimates that reach it indicate
#'   the Poisson limit (no overdispersion).
#' @param theta_warn A warning is issued when the estimated theta exceeds this
#'   value.
#'
#' @return A list of control parameters.
#' @examples
#' nbreg_control(epsilon = 1e-10, trace = TRUE)
#' @export
nbreg_control <- function(epsilon = 1e-8, maxit = 100, trace = FALSE,
                          theta_max = 1e8, theta_warn = 1e5) {
  if (!is.numeric(epsilon) || length(epsilon) != 1 || epsilon <= 0) {
    stop("'epsilon' must be a single positive number", call. = FALSE)
  }
  if (!is.numeric(maxit) || length(maxit) != 1 || maxit < 1) {
    stop("'maxit' must be a single number >= 1", call. = FALSE)
  }
  if (!is.numeric(theta_max) || theta_max <= 0 ||
      !is.numeric(theta_warn) || theta_warn <= 0) {
    stop("'theta_max' and 'theta_warn' must be positive", call. = FALSE)
  }
  list(epsilon = epsilon, maxit = as.integer(maxit), trace = isTRUE(trace),
       theta_max = theta_max, theta_warn = theta_warn)
}

#' Negative binomial regression by IRLS and Newton-Raphson
#'
#' Fits the NB2 regression model `log(mu) = X beta + offset`,
#' `Var(Y) = mu + mu^2 / theta`. The regression coefficients are updated by
#' iteratively reweighted least squares (Fisher scoring) with theta held
#' fixed, and `log(theta)` is updated by Newton-Raphson with beta held fixed.
#' The two steps alternate until the log-likelihood converges. No
#' general-purpose optimizer is used. The results reproduce
#' `MASS::glm.nb()`.
#'
#' @param formula A model formula, e.g. `y ~ x1 + x2`. Offsets may be given
#'   in the formula with `offset()`.
#' @param data A data frame containing the variables in `formula`.
#' @param offset Optional offset for the linear predictor (e.g. log exposure),
#'   evaluated in `data`.
#' @param control A list from [nbreg_control()].
#' @param weights Optional non-negative prior weights.
#' @param subset Optional expression selecting the rows of `data` to use.
#' @param na.action How to handle missing values (default from
#'   `getOption("na.action")`).
#'
#' @return An object of class `"nbreg"`: a list containing
#'   `coefficients`, `vcov` (variance-covariance matrix of the coefficients,
#'   treating theta as known), `theta`, `SE.theta`, `loglik`, `iter`,
#'   `converged`, `fitted.values`, `linear.predictors`, `call`, and the
#'   model information needed by [predict.nbreg()].
#' @seealso [nbreg-methods] for summary(), predict() and other methods.
#' @examplesIf requireNamespace("pscl", quietly = TRUE)
#' data("bioChemists", package = "pscl")
#' fit <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#' summary(fit)
#' AIC(fit)
#' @export
nbreg <- function(formula, data, offset, control = nbreg_control(),
                  weights, subset, na.action) {
  cl <- match.call()
  if (missing(data)) data <- environment(formula)

  # Build the model frame the same way glm() does, so factor coding,
  # subsets, weights, and offsets match the reference functions.
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y <- stats::model.response(mf, "numeric")
  if (is.null(y)) stop("the formula must have a response", call. = FALSE)
  check_counts(y)
  X <- stats::model.matrix(mt, mf)
  n <- NROW(y)

  off <- as.vector(stats::model.offset(mf))
  if (is.null(off)) off <- rep(0, n)
  a <- as.vector(stats::model.weights(mf))
  if (is.null(a)) a <- rep(1, n)
  if (any(a < 0)) stop("weights must be non-negative", call. = FALSE)

  fit <- nb_fit(X, y, offset = off, weights = a, control = control)

  fit$call <- cl
  fit$terms <- mt
  fit$model <- mf
  fit$xlevels <- stats::.getXlevels(mt, mf)
  fit$contrasts <- attr(X, "contrasts")
  class(fit) <- "nbreg"
  fit
}

# ---------------------------------------------------------------------------
# Internal fitting engine. Also used for the count-part M-step of zinbreg(),
# where the prior weights are a_i = 1 - tau_i.
# ---------------------------------------------------------------------------
nb_fit <- function(X, y, offset = rep(0, NROW(y)), weights = rep(1, NROW(y)),
                   control = nbreg_control(), start = NULL, theta_start = NULL,
                   warn = TRUE) {
  a <- weights
  eps <- control$epsilon

  # (c) Starting values: a Poisson IRLS fit for beta, and the moment-type
  #     estimate theta0 = sum(a) / sum(a (y/mu - 1)^2) (as in MASS).
  beta <- if (is.null(start)) pois_irls(X, y, offset, a, control) else start
  mu <- exp(drop(X %*% beta) + offset)
  theta <- if (is.null(theta_start)) {
    sum(a) / sum(a * (y / mu - 1)^2)
  } else {
    theta_start
  }
  # The moment estimate can be extreme (e.g. Inf when y/mu - 1 is ~0 for
  # every row, i.e. no overdispersion), so keep it inside [1e-4, theta_max].
  theta <- min(max(theta, 1e-4), control$theta_max)

  ll <- nb_loglik(y, mu, theta, a)
  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(control$maxit)) {
    ll_old <- ll

    # (a) IRLS for beta with theta fixed
    b <- nb_irls_beta(X, y, offset, a, theta, beta, control)
    beta <- b$beta
    mu <- b$mu

    # (b) Newton-Raphson for phi = log(theta) with beta fixed
    theta <- nb_theta_nr(y, mu, a, theta, control)

    ll <- nb_loglik(y, mu, theta, a)
    if (control$trace) {
      cat(sprintf("nbreg iter %3d: logLik = %.10f, theta = %.6f\n",
                  iter, ll, theta))
    }
    # (c) Outer convergence: stop when |l(beta, theta) - l_old| < epsilon.
    if (abs(ll - ll_old) < eps) {
      converged <- TRUE
      break
    }
  }
  if (warn && !converged) {
    warning("nbreg: alternation did not converge in ", control$maxit,
            " iterations", call. = FALSE)
  }
  if (warn && theta > control$theta_warn) {
    warning("nbreg: theta is very large (", signif(theta, 4), "); the data ",
            "show little overdispersion and a Poisson model may suffice",
            call. = FALSE)
  }

  # (d) Standard errors at convergence.
  # Var(beta) = (X' W X)^{-1}, W = diag(a theta mu / (theta + mu)),
  # treating theta as known (as glm.nb does).
  eta <- drop(X %*% beta) + offset
  w <- a * theta * mu / (theta + mu)
  XtWX <- crossprod(X * sqrt(w))
  V <- solve(XtWX)
  dimnames(V) <- list(colnames(X), colnames(X))
  # SE(theta) = {-H(theta)}^{-1/2}
  H <- nb_theta_hess(theta, y, mu, a)
  se_theta <- if (H < 0) 1 / sqrt(-H) else NA_real_

  list(coefficients = beta,
       vcov = V,
       theta = theta,
       SE.theta = se_theta,
       loglik = ll,
       iter = iter,
       converged = converged,
       fitted.values = mu,
       linear.predictors = eta,
       residuals = y - mu,
       y = y,
       prior.weights = a,
       offset = offset,
       rank = ncol(X),
       df.residual = sum(a > 0) - ncol(X) - 1,
       control = control)
}

# IRLS for beta with theta fixed (Section 3.1a).
#   w_i = a_i * theta mu_i / (theta + mu_i)
#   z_i = (eta_i - o_i) + (y_i - mu_i) / mu_i
#   beta <- (X'WX)^{-1} X'W z
# Step halving guards against a (rare) decrease of the log-likelihood.
# theta is fixed here, so only the mu-dependent kernel of the log-likelihood
# (nb_kernel) is needed to compare candidate steps.
nb_irls_beta <- function(X, y, offset, a, theta, beta, control) {
  eta <- drop(X %*% beta) + offset
  mu <- exp(eta)
  ll <- nb_kernel(y, mu, theta, a)
  for (it in seq_len(control$maxit)) {
    w <- a * theta * mu / (theta + mu)
    z <- (eta - offset) + (y - mu) / mu
    beta_new <- wls_solve(X, w, z)
    eta_new <- drop(X %*% beta_new) + offset
    mu_new <- exp(eta_new)
    ll_new <- nb_kernel(y, mu_new, theta, a)
    k <- 0L
    while ((!is.finite(ll_new) || ll_new < ll - 1e-12) && k < 30L) {
      beta_new <- (beta + beta_new) / 2
      eta_new <- drop(X %*% beta_new) + offset
      mu_new <- exp(eta_new)
      ll_new <- nb_kernel(y, mu_new, theta, a)
      k <- k + 1L
    }
    # Inner convergence uses a tolerance 100x tighter than the outer loop, so
    # beta is fully converged for this theta before theta is updated. We
    # also stop once beta stops moving: for a large |loglik| a very small
    # absolute tolerance can be below floating-point resolution and would
    # otherwise never be met.
    done <- abs(ll_new - ll) < control$epsilon * 1e-2 ||
      max(abs(beta_new - beta)) < 1e-10
    beta <- beta_new
    eta <- eta_new
    mu <- mu_new
    ll <- ll_new
    if (done) break
  }
  list(beta = beta, mu = mu)
}

# Newton-Raphson for phi = log(theta) with mu fixed (Section 3.1b).
#   U_phi = theta U(theta),  H_phi = theta^2 H(theta) + theta U(theta)
#   phi <- phi - U_phi / H_phi
# Working on the log scale keeps theta > 0. If H_phi is not negative
# (log-likelihood not locally concave) we take a unit step uphill instead,
# and step halving guarantees the log-likelihood does not decrease.
nb_theta_nr <- function(y, mu, a, theta, control) {
  phi <- log(theta)
  phi_max <- log(control$theta_max)
  ll <- nb_loglik(y, mu, theta, a)
  for (it in seq_len(control$maxit)) {
    U <- nb_theta_score(theta, y, mu, a)
    H <- nb_theta_hess(theta, y, mu, a)
    U_phi <- theta * U
    H_phi <- theta^2 * H + theta * U
    step <- if (H_phi < 0) -U_phi / H_phi else sign(U_phi)
    # Cap the step at |5| on the log scale (theta changes by at most a factor
    # e^5 ~ 148 per step), so one bad Newton step cannot overshoot wildly.
    step <- max(min(step, 5), -5)
    phi_new <- min(phi + step, phi_max)
    ll_new <- nb_loglik(y, mu, exp(phi_new), a)
    k <- 0L
    while ((!is.finite(ll_new) || ll_new < ll - 1e-12) && k < 30L) {
      phi_new <- (phi + phi_new) / 2
      ll_new <- nb_loglik(y, mu, exp(phi_new), a)
      k <- k + 1L
    }
    done <- abs(phi_new - phi) < 1e-10 || phi_new >= phi_max
    phi <- phi_new
    theta <- exp(phi)
    ll <- ll_new
    if (done) break
  }
  theta
}

# Poisson IRLS used only for starting values:
#   w_i = a_i mu_i,  z_i = (eta_i - o_i) + (y_i - mu_i) / mu_i
pois_irls <- function(X, y, offset, a, control) {
  # Start from mu = y + 0.1 (as glm() does for Poisson): adding 0.1 avoids
  # log(0) and division by zero in z_i when y_i = 0.
  mu <- y + 0.1
  eta <- log(mu)
  ll <- -Inf
  beta <- NULL
  for (it in seq_len(control$maxit)) {
    w <- a * mu
    z <- (eta - offset) + (y - mu) / mu
    beta <- wls_solve(X, w, z)
    eta <- drop(X %*% beta) + offset
    mu <- exp(eta)
    ll_new <- pois_loglik(y, mu, a)
    if (abs(ll_new - ll) < control$epsilon) break
    ll <- ll_new
  }
  beta
}
