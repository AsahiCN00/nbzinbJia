# Simulation study (Section 5.1): does zinbreg() recover known parameters?
#
# Design: R = 200 data sets of size n = 500 from a ZINB model with two
# covariates in each part.
#   count part:  log(mu) = b0 + b1 x1 + b2 x2,  x1 ~ N(0, 1), x2 ~ Bern(0.5)
#   zero part:   logit(pi) = g0 + g1 z1 + g2 x2,  z1 ~ N(0, 1)
#   true values: beta = (0.8, 0.5, -0.6), gamma = (-1.0, 0.7, 0.5),
#                theta = 1.5 (log theta = 0.405)
# For each parameter we report the bias, the mean of the estimated SEs vs
# the Monte Carlo SD of the estimates, and the empirical coverage of the
# 95% Wald interval (estimate +/- 1.96 SE).
#
# Run from the package root after installing the package:
#   Rscript inst/benchmark/simulation.R
# Results are saved to inst/benchmark/simulation_results.rds.

library(nbzinbJia)

set.seed(651)
R <- 200
n <- 500
beta <- c("(Intercept)" = 0.8, x1 = 0.5, x2 = -0.6)
gamma <- c("(Intercept)" = -1.0, z1 = 0.7, x2 = 0.5)
theta <- 1.5
truth <- c(setNames(beta, paste0("count_", names(beta))),
           setNames(gamma, paste0("zero_", names(gamma))),
           "log(theta)" = log(theta))

simulate_zinb <- function(n) {
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  z1 <- rnorm(n)
  mu <- exp(beta[1] + beta[2] * x1 + beta[3] * x2)
  p <- plogis(gamma[1] + gamma[2] * z1 + gamma[3] * x2)
  structural_zero <- runif(n) < p
  y <- ifelse(structural_zero, 0, rnbinom(n, size = theta, mu = mu))
  data.frame(y, x1, x2, z1)
}

est <- se <- matrix(NA_real_, R, length(truth),
                    dimnames = list(NULL, names(truth)))
iter <- integer(R)
converged <- logical(R)
n_warn <- integer(R)
prop_zero <- numeric(R)

t0 <- proc.time()
for (r in seq_len(R)) {
  d <- simulate_zinb(n)
  prop_zero[r] <- mean(d$y == 0)
  warns <- character(0)
  fit <- withCallingHandlers(
    zinbreg(y ~ x1 + x2 | z1 + x2, data = d),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  est[r, ] <- c(coef(fit), log(fit$theta))
  se[r, ] <- c(sqrt(diag(vcov(fit))), fit$SE.logtheta)
  iter[r] <- fit$iter
  converged[r] <- fit$converged
  n_warn[r] <- length(warns)
}
elapsed <- (proc.time() - t0)[["elapsed"]]

lo <- est - 1.96 * se
hi <- est + 1.96 * se
covered <- sweep(lo, 2, truth, "<=") & sweep(hi, 2, truth, ">=")

summary_table <- data.frame(
  parameter = names(truth),
  true = unname(truth),
  mean_estimate = colMeans(est),
  bias = colMeans(est) - truth,
  mean_se = colMeans(se),
  mc_sd = apply(est, 2, sd),
  coverage = colMeans(covered),
  row.names = NULL)

out <- list(summary = summary_table, est = est, se = se, iter = iter,
            converged = converged, n_warnings = n_warn,
            prop_zero = prop_zero, truth = truth,
            design = list(R = R, n = n, beta = beta, gamma = gamma,
                          theta = theta, seed = 651),
            elapsed_sec = elapsed, date = Sys.Date())
saveRDS(out, file.path("inst", "benchmark", "simulation_results.rds"))

print(summary_table, digits = 3)
cat(sprintf(paste0("\nconverged: %d/%d, fits with warnings: %d, ",
                   "median EM iterations: %d, mean proportion of zeros: ",
                   "%.3f, time: %.1f s\n"),
            sum(converged), R, sum(n_warn > 0), as.integer(median(iter)),
            mean(prop_zero), elapsed))
# Monte Carlo standard error of a 95% coverage estimate with R = 200:
cat(sprintf("MC SE of coverage: %.3f\n", sqrt(0.95 * 0.05 / R)))
