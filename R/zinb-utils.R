# Internal helpers for the ZINB model (Sections 2.2 and 3.2).
# Notation: mu_i = exp(x_i'beta + o_i) (count part),
#           pi_i = plogis(z_i'gamma + o_zi) (zero-inflation part),
#           phi = log(theta), a_i = prior weights.

# log(exp(u) + exp(v)) computed without overflow/underflow.
log_add_exp <- function(u, v) {
  m <- pmax(u, v)
  m + log1p(exp(-abs(u - v)))
}

# Per-observation observed-data log-likelihood (Section 2.2):
#   y_i = 0: log[ pi_i + (1 - pi_i) f_NB(0; mu_i, theta) ]
#   y_i > 0: log(1 - pi_i) + log f_NB(y_i; mu_i, theta)
# log(pi) and log(1 - pi) come from plogis(..., log.p = TRUE) so that pi
# close to 0 or 1 does not lose precision.
zinb_loglik_i <- function(y, mu, theta, eta_z) {
  log_pi <- stats::plogis(eta_z, log.p = TRUE)
  log_1mpi <- stats::plogis(-eta_z, log.p = TRUE)
  log_f <- stats::dnbinom(y, size = theta, mu = mu, log = TRUE)
  ifelse(y == 0, log_add_exp(log_pi, log_1mpi + log_f), log_1mpi + log_f)
}

zinb_loglik <- function(y, mu, theta, eta_z, a) {
  sum(a * zinb_loglik_i(y, mu, theta, eta_z))
}

# E-step (Section 3.2): posterior probability of a structural zero
#   tau_i = pi_i / [ pi_i + (1 - pi_i) f_NB(0; mu_i, theta) ]  if y_i = 0
#   tau_i = 0                                                 if y_i > 0
# computed on the log scale as exp(log pi_i - log[denominator]).
zinb_estep <- function(y, mu, theta, eta_z) {
  log_pi <- stats::plogis(eta_z, log.p = TRUE)
  log_1mpi <- stats::plogis(-eta_z, log.p = TRUE)
  log_f0 <- stats::dnbinom(0, size = theta, mu = mu, log = TRUE)
  tau <- exp(log_pi - log_add_exp(log_pi, log_1mpi + log_f0))
  tau[y > 0] <- 0
  tau
}

# M-step, zero part (Section 3.2): maximize
#   sum_i a_i [ tau_i log pi_i + (1 - tau_i) log(1 - pi_i) ]
# over gamma: a logistic regression with fractional response tau_i, fitted
# by IRLS with
#   w_i = a_i pi_i (1 - pi_i),
#   z_i = (eta_i - o_i) + (tau_i - pi_i) / [ pi_i (1 - pi_i) ].
# Also used (with tau = I(y = 0)) for the starting value of gamma.
logit_irls <- function(Z, tau, offset, a, gamma, control) {
  obj <- function(eta) {
    sum(a * (tau * stats::plogis(eta, log.p = TRUE) +
               (1 - tau) * stats::plogis(-eta, log.p = TRUE)))
  }
  eta <- drop(Z %*% gamma) + offset
  ll <- obj(eta)
  for (it in seq_len(control$maxit)) {
    p <- stats::plogis(eta)
    # pi (1 - pi) is floored at machine epsilon so z_i stays finite when
    # pi_i is numerically 0 or 1 (those rows get ~zero weight anyway).
    v <- pmax(p * (1 - p), .Machine$double.eps)
    w <- a * v
    z <- (eta - offset) + (tau - p) / v
    gamma_new <- wls_solve(Z, w, z)
    eta_new <- drop(Z %*% gamma_new) + offset
    ll_new <- obj(eta_new)
    # step halving: the objective must not decrease
    k <- 0L
    while ((!is.finite(ll_new) || ll_new < ll - 1e-12) && k < 30L) {
      gamma_new <- (gamma + gamma_new) / 2
      eta_new <- drop(Z %*% gamma_new) + offset
      ll_new <- obj(eta_new)
      k <- k + 1L
    }
    done <- abs(ll_new - ll) < control$epsilon * 1e-2
    gamma <- gamma_new
    eta <- eta_new
    ll <- ll_new
    if (done) break
  }
  gamma
}

# Offset of one model part, taken from a model frame: the offset() terms
# of that part's formula, plus the 'offset' argument when include_arg is
# TRUE (zeroinfl() applies the argument to the count part only).
terms_offset <- function(mf, terms, include_arg) {
  off <- rep(0, nrow(mf))
  if (include_arg && !is.null(mf[["(offset)"]])) off <- off + mf[["(offset)"]]
  idx <- attr(terms, "offset")
  if (length(idx)) {
    vars <- attr(terms, "variables")
    for (i in idx) {
      off <- off + mf[[paste(deparse(vars[[i + 1L]], width.cutoff = 500L),
                              collapse = "")]]
    }
  }
  off
}
