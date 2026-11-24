# Compare zinbreg() with pscl::zeroinfl(dist = "negbin") using the
# tolerances of Table 4: coefficients, SEs, log(theta) and the
# log-likelihood all within 1e-3 (maximum absolute difference).
max_abs_diff <- function(a, b) max(abs(unname(a) - unname(b)))

# All coefficients of a zinbreg fit, named like coef(zeroinfl(...)).
zcoef <- function(fit) {
  cf <- fit$coefficients
  stats::setNames(c(cf$count, cf$zero),
                  c(paste0("count_", names(cf$count)),
                    paste0("zero_", names(cf$zero))))
}

# SEs of the coefficients (fit$vcov also has a row for log(theta)).
zse <- function(fit) {
  k <- length(fit$coefficients$count) + length(fit$coefficients$zero)
  sqrt(diag(fit$vcov))[seq_len(k)]
}

expect_matches_zeroinfl <- function(fit, ref) {
  expect_true(fit$converged)
  expect_identical(names(zcoef(fit)), names(coef(ref)))
  expect_lte(max_abs_diff(zcoef(fit), coef(ref)), 1e-3)
  expect_lte(max_abs_diff(zse(fit), sqrt(diag(vcov(ref)))), 1e-3)
  expect_lte(max_abs_diff(log(fit$theta), log(ref$theta)), 1e-3)
  expect_lte(max_abs_diff(fit$SE.logtheta, ref$SE.logtheta), 1e-3)
  expect_lte(max_abs_diff(fit$loglik, logLik(ref)), 1e-3)
}

ref_zinb <- function(f, data, ...) {
  pscl::zeroinfl(f, data = data, dist = "negbin",
                 control = pscl::zeroinfl.control(reltol = 1e-12), ...)
}

test_that("zinbreg matches zeroinfl on bioChemists", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  f <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
  expect_matches_zeroinfl(zinbreg(f, data = bioChemists),
                          ref_zinb(f, bioChemists))
})

test_that("zinbreg matches zeroinfl on NMES1988", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER")
  f <- visits ~ hospital + health + chronic + gender + school + insurance |
    hospital + chronic + gender + school + insurance
  expect_matches_zeroinfl(zinbreg(f, data = NMES1988), ref_zinb(f, NMES1988))
})

test_that("EM never decreases the observed-data log-likelihood", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- zinbreg(art ~ fem + kid5 + ment | ment, data = bioChemists)
  expect_true(all(diff(fit$loglik_path) > -1e-8))
  # at the MLE the score of the observed-data log-likelihood is ~0
  expect_lt(max(abs(fit$gradient)), 1e-2)
})

test_that("analytic Hessian agrees with a numerical Hessian", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("numDeriv")
  data("bioChemists", package = "pscl")
  d <- bioChemists
  fit <- zinbreg(art ~ fem + kid5 + ment | kid5 + ment, data = d)
  X <- model.matrix(~ fem + kid5 + ment, d)
  Z <- model.matrix(~ kid5 + ment, d)
  kx <- ncol(X)
  kz <- ncol(Z)
  a <- rep(1, nrow(d))
  o <- rep(0, nrow(d))
  llfun <- function(p) {
    zinb_loglik(d$art, exp(drop(X %*% p[1:kx])), exp(p[kx + kz + 1]),
                drop(Z %*% p[kx + 1:kz]), a)
  }
  par <- c(fit$coefficients$count, fit$coefficients$zero, log(fit$theta))
  sh <- zinb_score_hessian(X, Z, d$art, a, fit$coefficients$count,
                           fit$coefficients$zero, fit$theta, o, o)
  Hn <- numDeriv::hessian(llfun, par)
  expect_lt(max(abs(sh$hessian - Hn)) / max(abs(Hn)), 1e-6)
  expect_lt(max(abs(sh$gradient - numDeriv::grad(llfun, par))), 1e-4)
  # and away from the MLE too (the formulas must hold everywhere)
  par2 <- par + 0.1
  sh2 <- zinb_score_hessian(X, Z, d$art, a, par2[1:kx], par2[kx + 1:kz],
                            exp(par2[kx + kz + 1]), o, o)
  Hn2 <- numDeriv::hessian(llfun, par2)
  expect_lt(max(abs(sh2$hessian - Hn2)) / max(abs(Hn2)), 1e-6)
})

test_that("a formula without '|' uses the same covariates in both parts", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  f1 <- zinbreg(art ~ fem + ment, data = bioChemists)
  f2 <- zinbreg(art ~ fem + ment | fem + ment, data = bioChemists)
  expect_equal(zcoef(f1), zcoef(f2))
  expect_equal(names(f1$coefficients$zero), c("(Intercept)", "femWomen",
                                              "ment"))
})

test_that("offsets match zeroinfl", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  d <- bioChemists
  d$expo <- 1 + d$kid5 + (d$mar == "Married")
  ref <- pscl::zeroinfl(art ~ fem + ment | ment, data = d, dist = "negbin",
                        offset = log(expo),
                        control = pscl::zeroinfl.control(reltol = 1e-12))
  fit1 <- zinbreg(art ~ fem + ment | ment, data = d, offset = log(expo))
  fit2 <- zinbreg(art ~ fem + ment + offset(log(expo)) | ment, data = d)
  expect_matches_zeroinfl(fit1, ref)
  expect_matches_zeroinfl(fit2, ref)
})

test_that("single-covariate and intercept-only zero models work", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("AER")
  data("bioChemists", package = "pscl")
  data("NMES1988", package = "AER")
  f1 <- art ~ ment | ment
  expect_matches_zeroinfl(zinbreg(f1, data = bioChemists),
                          ref_zinb(f1, bioChemists))
  f2 <- visits ~ hospital + chronic | 1
  expect_matches_zeroinfl(zinbreg(f2, data = NMES1988),
                          ref_zinb(f2, NMES1988))
})

test_that("no zero inflation (pi -> 0 boundary) gives a warning", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  # With an intercept-only zero part these data show no zero inflation:
  # the MLE has gamma -> -Inf, but the log-likelihood and count part are
  # still those of zeroinfl().
  f <- art ~ ment | 1
  expect_warning(fit <- zinbreg(f, data = bioChemists), "essentially zero")
  ref <- ref_zinb(f, bioChemists)
  expect_lte(max_abs_diff(fit$coefficients$count,
                          coef(ref, model = "count")), 1e-3)
  expect_lte(max_abs_diff(fit$loglik, logLik(ref)), 1e-3)
})

test_that("bad input gives informative errors", {
  d <- data.frame(y = c(0, 1, 2, 0, 5, 1), x = c(1, 2, 3, 4, 5, 6))
  expect_error(zinbreg(y ~ x, data = transform(d, y = y + 1)), "no zeros")
  expect_error(zinbreg(y ~ x, data = transform(d, y = y - 1)),
               "non-negative")
  expect_error(zinbreg(y ~ x | x | x, data = d), "at most two parts")
  expect_error(zinbreg_control(maxit = 0), "maxit")
})

test_that("hitting maxit gives a warning", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  expect_warning(zinbreg(art ~ ment | ment, data = bioChemists,
                         control = zinbreg_control(maxit = 2)),
                 "did not converge")
})
