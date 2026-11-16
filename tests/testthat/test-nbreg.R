# Compare nbreg() with MASS::glm.nb() using the tolerances of Table 4:
# coefficients and SEs <= 1e-5, theta and log-likelihood <= 1e-4.
max_abs_diff <- function(a, b) max(abs(unname(a) - unname(b)))

expect_matches_glm_nb <- function(fit, ref) {
  expect_true(fit$converged)
  expect_identical(names(coef(fit)), names(coef(ref)))
  expect_lte(max_abs_diff(coef(fit), coef(ref)), 1e-5)
  expect_lte(max_abs_diff(sqrt(diag(vcov(fit))), sqrt(diag(vcov(ref)))), 1e-5)
  expect_lte(max_abs_diff(fit$theta, ref$theta), 1e-4)
  expect_lte(max_abs_diff(fit$SE.theta, ref$SE.theta), 1e-4)
  expect_lte(max_abs_diff(logLik(fit), logLik(ref)), 1e-4)
}

ref_control <- function() glm.control(epsilon = 1e-10, maxit = 100)

test_that("nbreg matches glm.nb on bioChemists", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  library(MASS)
  data("bioChemists", package = "pscl")
  f <- art ~ fem + mar + kid5 + phd + ment
  fit <- nbreg(f, data = bioChemists)
  ref <- glm.nb(f, data = bioChemists, control = ref_control())
  expect_matches_glm_nb(fit, ref)
})

test_that("nbreg matches glm.nb on NMES1988", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("AER")
  library(MASS)
  data("NMES1988", package = "AER")
  f <- visits ~ hospital + health + chronic + gender + school + insurance
  fit <- nbreg(f, data = NMES1988)
  ref <- glm.nb(f, data = NMES1988, control = ref_control())
  expect_matches_glm_nb(fit, ref)
})

test_that("offset argument and offset() in the formula agree with glm.nb", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  library(MASS)
  data("bioChemists", package = "pscl")
  d <- bioChemists
  d$expo <- 1 + d$kid5 + (d$mar == "Married")
  ref <- glm.nb(art ~ fem + ment + offset(log(expo)), data = d,
                control = ref_control())
  fit1 <- nbreg(art ~ fem + ment, data = d, offset = log(expo))
  fit2 <- nbreg(art ~ fem + ment + offset(log(expo)), data = d)
  expect_matches_glm_nb(fit1, ref)
  expect_matches_glm_nb(fit2, ref)
  # predict() must re-apply the offset for new data
  expect_equal(predict(fit1, newdata = d, type = "response"),
               fitted(fit1), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(predict(fit2, newdata = d, type = "response"),
               fitted(fit2), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("single-covariate and intercept-only models work", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  library(MASS)
  data("bioChemists", package = "pscl")
  expect_matches_glm_nb(nbreg(art ~ ment, data = bioChemists),
                        glm.nb(art ~ ment, data = bioChemists,
                               control = ref_control()))
  fit0 <- nbreg(art ~ 1, data = bioChemists)
  expect_equal(unname(exp(coef(fit0))), mean(bioChemists$art),
               tolerance = 1e-6)
})

test_that("prior weights give the same estimates as replicated rows", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  d <- bioChemists[1:200, ]
  w <- rep(1:2, length.out = nrow(d))
  fit_w <- nbreg(art ~ fem + ment, data = d, weights = w)
  fit_r <- nbreg(art ~ fem + ment, data = d[rep(seq_len(nrow(d)), w), ])
  expect_equal(coef(fit_w), coef(fit_r), tolerance = 1e-6)
  expect_equal(fit_w$theta, fit_r$theta, tolerance = 1e-6)
})

test_that("S3 methods return sensible output", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- nbreg(art ~ fem + ment, data = bioChemists)
  ll <- logLik(fit)
  expect_equal(attr(ll, "df"), 4L)
  expect_equal(nobs(fit), nrow(bioChemists))
  expect_equal(AIC(fit), -2 * as.numeric(ll) + 2 * 4)
  expect_equal(BIC(fit), -2 * as.numeric(ll) + log(nrow(bioChemists)) * 4)
  expect_equal(residuals(fit), bioChemists$art - fitted(fit))
  mu <- fitted(fit)
  expect_equal(residuals(fit, type = "pearson"),
               (bioChemists$art - mu) / sqrt(mu + mu^2 / fit$theta))
  expect_equal(predict(fit, newdata = bioChemists, type = "response"),
               fitted(fit), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(predict(fit), log(fitted(fit)), ignore_attr = TRUE)
  s <- summary(fit)
  expect_s3_class(s, "summary.nbreg")
  expect_equal(colnames(s$coefficients),
               c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_output(print(fit), "Theta")
  expect_output(print(s), "Pearson residuals")
})

test_that("bad input gives informative errors", {
  d <- data.frame(y = c(0, 1, 2, 3, 5, 1), x = c(1, 2, 3, 4, 5, 6))
  expect_error(nbreg(y ~ x, data = transform(d, y = y - 1)), "non-negative")
  expect_error(nbreg(y ~ x, data = transform(d, y = y + 0.5)), "integer")
  expect_error(nbreg(y ~ x + I(2 * x), data = d), "rank deficient")
  expect_error(nbreg(y ~ x, data = d, weights = rep(-1, 6)), "non-negative")
  expect_error(nbreg_control(epsilon = -1), "epsilon")
})

test_that("near-Poisson data trigger the large-theta warning", {
  set.seed(1)
  x <- rnorm(500)
  y <- rpois(500, exp(0.5 + 0.3 * x))
  expect_warning(nbreg(y ~ x, data = data.frame(x, y)), "theta is very large")
})

test_that("hitting maxit gives a warning", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  expect_warning(nbreg(art ~ fem + ment, data = bioChemists,
                       control = nbreg_control(maxit = 1)),
                 "did not converge")
})
