# Benchmark nbzinbJia against the reference implementations (Section 5).
#
#   nbreg()   vs MASS::glm.nb()
#   zinbreg() vs pscl::zeroinfl(dist = "negbin")
#
# on pscl::bioChemists and AER::NMES1988: side-by-side estimates and
# standard errors, the maximum absolute differences against the tolerances
# of Table 4, and computation time.
#
# Run from the package root after installing the package:
#   Rscript inst/benchmark/benchmark.R
# Results are saved to inst/benchmark/benchmark_results.rds.

library(nbzinbJia)
library(MASS)
library(pscl)
data("bioChemists", package = "pscl")
data("NMES1988", package = "AER")

# Reference fits use tight tolerances (as recommended in Section 5).
ref_nb <- function(f, d) {
  glm.nb(f, data = d, control = glm.control(epsilon = 1e-10, maxit = 100))
}
ref_zinb <- function(f, d) {
  zeroinfl(f, data = d, dist = "negbin",
           control = zeroinfl.control(reltol = 1e-12))
}

# Side-by-side table of estimates and SEs (ours vs reference).
side_by_side <- function(est, se, ref_est, ref_se) {
  data.frame(term = names(ref_est),
             estimate = unname(est), ref_estimate = unname(ref_est),
             abs_diff_est = abs(unname(est) - unname(ref_est)),
             se = unname(se), ref_se = unname(ref_se),
             abs_diff_se = abs(unname(se) - unname(ref_se)),
             row.names = NULL)
}

compare_nb <- function(f, d) {
  fit <- nbreg(f, data = d)
  ref <- ref_nb(f, d)
  tab <- side_by_side(coef(fit), sqrt(diag(vcov(fit))),
                      coef(ref), sqrt(diag(vcov(ref))))
  tab <- rbind(tab, data.frame(term = "theta",
                               estimate = fit$theta, ref_estimate = ref$theta,
                               abs_diff_est = abs(fit$theta - ref$theta),
                               se = fit$SE.theta, ref_se = ref$SE.theta,
                               abs_diff_se = abs(fit$SE.theta - ref$SE.theta)))
  summary <- data.frame(
    quantity = c("Coefficients", "Standard errors", "theta", "Log-likelihood"),
    max_abs_diff = c(max(abs(coef(fit) - coef(ref))),
                     max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(ref))))),
                     abs(fit$theta - ref$theta),
                     abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref)))),
    tolerance = c(1e-5, 1e-5, 1e-4, 1e-4))
  summary$pass <- summary$max_abs_diff <= summary$tolerance
  list(table = tab, summary = summary,
       loglik = c(ours = as.numeric(logLik(fit)),
                  reference = as.numeric(logLik(ref))),
       iter = fit$iter)
}

compare_zinb <- function(f, d) {
  fit <- zinbreg(f, data = d)
  ref <- ref_zinb(f, d)
  tab <- side_by_side(coef(fit), sqrt(diag(vcov(fit))),
                      coef(ref), sqrt(diag(vcov(ref))))
  tab <- rbind(tab, data.frame(term = "log(theta)",
                               estimate = log(fit$theta),
                               ref_estimate = log(ref$theta),
                               abs_diff_est = abs(log(fit$theta) -
                                                    log(ref$theta)),
                               se = fit$SE.logtheta, ref_se = ref$SE.logtheta,
                               abs_diff_se = abs(fit$SE.logtheta -
                                                   ref$SE.logtheta)))
  summary <- data.frame(
    quantity = c("Coefficients", "Standard errors", "log(theta)",
                 "Log-likelihood"),
    max_abs_diff = c(max(abs(coef(fit) - coef(ref))),
                     max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(ref))))),
                     abs(log(fit$theta) - log(ref$theta)),
                     abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref)))),
    tolerance = rep(1e-3, 4))
  summary$pass <- summary$max_abs_diff <= summary$tolerance
  list(table = tab, summary = summary,
       loglik = c(ours = as.numeric(logLik(fit)),
                  reference = as.numeric(logLik(ref))),
       iter = fit$iter)
}

models <- list(
  bioChemists = list(
    data = bioChemists,
    nb = art ~ fem + mar + kid5 + phd + ment,
    zinb = art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment),
  NMES1988 = list(
    data = NMES1988,
    nb = visits ~ hospital + health + chronic + gender + school + insurance,
    zinb = visits ~ hospital + health + chronic + gender + school + insurance |
      hospital + chronic + gender + school + insurance))

results <- list()
for (nm in names(models)) {
  m <- models[[nm]]
  message("Agreement: ", nm)
  results[[nm]] <- list(nb = compare_nb(m$nb, m$data),
                        zinb = compare_zinb(m$zinb, m$data))
}

# Computation time (median of several runs). Our code is not expected to be
# faster: EM converges linearly and every M-step is a full IRLS fit.
timing <- list()
for (nm in names(models)) {
  m <- models[[nm]]
  message("Timing: ", nm)
  n_it <- if (nm == "NMES1988") 3 else 10
  bm <- bench::mark(
    nbreg = nbreg(m$nb, data = m$data),
    glm.nb = ref_nb(m$nb, m$data),
    zinbreg = zinbreg(m$zinb, data = m$data),
    zeroinfl = ref_zinb(m$zinb, m$data),
    iterations = n_it, check = FALSE, memory = FALSE, filter_gc = FALSE)
  timing[[nm]] <- data.frame(data = nm, fn = as.character(bm$expression),
                             median_sec = as.numeric(bm$median),
                             iterations = n_it)
}
timing <- do.call(rbind, timing)
rownames(timing) <- NULL

out <- list(results = results, timing = timing,
            date = Sys.Date(), session = utils::sessionInfo())
saveRDS(out, file.path("inst", "benchmark", "benchmark_results.rds"))

for (nm in names(results)) {
  cat("\n=====", nm, "=====\n")
  cat("-- NB vs glm.nb\n"); print(results[[nm]]$nb$summary, digits = 3)
  cat("-- ZINB vs zeroinfl\n"); print(results[[nm]]$zinb$summary, digits = 3)
}
cat("\n-- Timing (median seconds)\n"); print(timing, digits = 3)
