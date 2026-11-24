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

# Derivatives of g_i = log f_NB(y_i; mu_i, theta) with respect to the count
# linear predictor eta_c (mu = exp(eta_c)) and phi = log(theta):
#   g_c   = theta (y - mu) / (theta + mu)
#   g_cc  = -theta mu (theta + y) / (theta + mu)^2
#   g_phi = theta U_i(theta)
#   g_pp  = theta^2 H_i(theta) + theta U_i(theta)
#   g_cp  = theta mu (y - mu) / (theta + mu)^2
# where U_i and H_i are the per-observation terms of Section 3.1(b).
nb_log_density_derivs <- function(y, mu, theta) {
  tm <- theta + mu
  U <- digamma(y + theta) - digamma(theta) + log(theta) + 1 - log(tm) -
    (y + theta) / tm
  H <- trigamma(y + theta) - trigamma(theta) + 1 / theta - 2 / tm +
    (y + theta) / tm^2
  list(c = theta * (y - mu) / tm,
       cc = -theta * mu * (theta + y) / tm^2,
       p = theta * U,
       pp = theta^2 * H + theta * U,
       cp = theta * mu * (y - mu) / tm^2)
}

# Score and analytic Hessian of the observed-data ZINB log-likelihood with
# respect to (beta, gamma, phi = log theta).
#
# Write l_i = log[ pi_i 1(y_i = 0) + (1 - pi_i) f(y_i) ] and g_i = log f(y_i).
# Differentiating, with tau_i the E-step posterior (tau_i = 0 when y_i > 0),
# gives one formula for every observation:
#   dl/d eta_c = (1 - tau) g_c          dl/d eta_z = tau - pi
#   dl/d phi   = (1 - tau) g_phi
#   d2l/(d eta_c)^2   = (1 - tau) g_cc + tau (1 - tau) g_c^2
#   d2l/d eta_c d phi = (1 - tau) g_cp + tau (1 - tau) g_c g_phi
#   d2l/d phi^2       = (1 - tau) g_pp + tau (1 - tau) g_phi^2
#   d2l/(d eta_z)^2   = tau (1 - tau) - pi (1 - pi)
#   d2l/d eta_c d eta_z = -tau (1 - tau) g_c
#   d2l/d eta_z d phi   = -tau (1 - tau) g_phi
# and the chain rule (d eta_c / d beta = x_i, d eta_z / d gamma = z_i)
# assembles the gradient and Hessian below.
zinb_score_hessian <- function(X, Z, y, a, beta, gamma, theta, offx, offz) {
  mu <- exp(drop(X %*% beta) + offx)
  eta_z <- drop(Z %*% gamma) + offz
  p <- stats::plogis(eta_z)
  tau <- zinb_estep(y, mu, theta, eta_z)
  d <- nb_log_density_derivs(y, mu, theta)
  tt <- tau * (1 - tau)

  s_c <- a * (1 - tau) * d$c
  s_z <- a * (tau - p)
  s_p <- a * (1 - tau) * d$p
  h_cc <- a * ((1 - tau) * d$cc + tt * d$c^2)
  h_cp <- a * ((1 - tau) * d$cp + tt * d$c * d$p)
  h_pp <- a * ((1 - tau) * d$pp + tt * d$p^2)
  h_zz <- a * (tt - p * (1 - p))
  h_cz <- a * (-tt * d$c)
  h_zp <- a * (-tt * d$p)

  grad <- c(colSums(X * s_c), colSums(Z * s_z), sum(s_p))
  Hbb <- crossprod(X, X * h_cc)
  Hgg <- crossprod(Z, Z * h_zz)
  Hbg <- crossprod(X, Z * h_cz)
  Hbp <- colSums(X * h_cp)
  Hgp <- colSums(Z * h_zp)
  Hpp <- sum(h_pp)
  H <- rbind(cbind(Hbb, Hbg, Hbp),
             cbind(t(Hbg), Hgg, Hgp),
             c(Hbp, Hgp, Hpp))
  list(gradient = unname(grad), hessian = unname(H))
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
