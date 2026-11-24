# Internal numerical helpers shared by nbreg() and zinbreg().
# None of these call an optimizer: every update is an explicit
# weighted-least-squares solve or a Newton-Raphson step.

# Weighted least squares: beta = (X'WX)^{-1} X'W z.
# Computed through the QR decomposition of W^{1/2} X, which is
# numerically more stable than forming X'WX directly.
wls_solve <- function(X, w, z) {
  sw <- sqrt(w)
  qrx <- qr(X * sw)
  if (qrx$rank < ncol(X)) {
    stop("design matrix is rank deficient; remove collinear covariates",
         call. = FALSE)
  }
  beta <- qr.coef(qrx, z * sw)
  names(beta) <- colnames(X)
  beta
}

# NB2 log-likelihood (Section 2.1), with prior weights a_i:
# l = sum_i a_i [ lgamma(y+theta) - lgamma(theta) - log y! + theta log theta
#                 + y log mu - (theta + y) log(theta + mu) ]
# dnbinom(size = theta, mu = mu) uses exactly this parameterization, and
# it evaluates lgamma(y + theta) - lgamma(theta) accurately even for huge
# theta (a direct lgamma difference loses precision there).
nb_loglik <- function(y, mu, theta, a) {
  sum(a * stats::dnbinom(y, size = theta, mu = mu, log = TRUE))
}

# The part of the NB log-likelihood that depends on mu:
#   theta log theta + y log mu - (theta + y) log(theta + mu) + const
#     = -theta log1p(mu / theta) + y log(mu / (theta + mu)) + const.
# With theta fixed, the lgamma terms are constants, so this kernel is all
# the beta-IRLS step needs to check that the log-likelihood increases. It
# is much cheaper than dnbinom(), which matters because it is evaluated at
# every IRLS step. (The y log(...) term is set to 0 when y = 0, avoiding
# 0 * -Inf.)
nb_kernel_i <- function(y, mu, theta) {
  ylog <- y * log(mu / (theta + mu))
  ylog[y == 0] <- 0
  -theta * log1p(mu / theta) + ylog
}

nb_kernel <- function(y, mu, theta, a) sum(a * nb_kernel_i(y, mu, theta))

# Poisson log-likelihood (used only for the Poisson starting values).
pois_loglik <- function(y, mu, a) {
  sum(a * stats::dpois(y, lambda = mu, log = TRUE))
}

# Score U(theta) of the NB log-likelihood with beta held fixed (Section 3.1b):
# U = sum a [ psi(y+theta) - psi(theta) + log theta + 1 - log(theta+mu)
#             - (y+theta)/(theta+mu) ]
nb_theta_score <- function(theta, y, mu, a) {
  sum(a * (digamma(y + theta) - digamma(theta) + log(theta) + 1 -
             log(theta + mu) - (y + theta) / (theta + mu)))
}

# Second derivative H(theta) (Section 3.1b):
# H = sum a [ psi'(y+theta) - psi'(theta) + 1/theta - 2/(theta+mu)
#             + (y+theta)/(theta+mu)^2 ]
nb_theta_hess <- function(theta, y, mu, a) {
  sum(a * (trigamma(y + theta) - trigamma(theta) + 1 / theta -
             2 / (theta + mu) + (y + theta) / (theta + mu)^2))
}

# Check that a response is a vector of non-negative integer counts.
check_counts <- function(y) {
  if (!is.numeric(y)) stop("response must be numeric counts", call. = FALSE)
  if (anyNA(y)) stop("response contains missing values", call. = FALSE)
  if (any(y < 0)) stop("response must be non-negative counts", call. = FALSE)
  if (any(abs(y - round(y)) > 1e-8)) {
    stop("response must be integer counts", call. = FALSE)
  }
  invisible(TRUE)
}
