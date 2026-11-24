# STA 651 final project: NB and ZINB regression R package (nbzinbJia)

## Rules (must follow)

- Fitting code in R/ must NOT call optim, optimize, optimise, nlm,
  nlminb, uniroot, glm, glm.fit, lm, lm.fit, lm.wfit, MASS::glm.nb,
  MASS::theta.ml, pscl::zeroinfl, countreg, glmmTMB, gamlss, VGAM,
  optimx, maxLik, stats4::mle, bbmle, or any other optimizer package.
- Allowed in fitting code: solve, qr, chol, crossprod, %\*%, lgamma,
  digamma, trigamma, dnbinom, plogis, qlogis, model.frame, model.matrix,
  model.response, model.offset, the Formula package, and
  numDeriv::hessian (standard errors only).
- NB: IRLS for beta alternated with Newton-Raphson on log(theta).
- ZINB: EM algorithm; M-steps are IRLS (weighted logistic with
  fractional response tau + NB fit with prior weights 1 - tau). EM must
  never decrease the observed-data log-likelihood.
- MASS::glm.nb, pscl::zeroinfl and glm may be used ONLY in tests/,
  inst/benchmark/, and the report, for comparison.
- Follow the equations in project-docs/STA651_Final_Project.pdf, Section
  3.
- Comment estimation code with the equation it implements.
- Explain your reasoning to me; I must be able to explain every line.

## Required agreement (max abs difference)

- NB vs glm.nb: coef, SE \<= 1e-5; theta, logLik \<= 1e-4.
- ZINB vs zeroinfl: coef, SE, log(theta), logLik \<= 1e-3.

## Workflow

- After editing roxygen comments: Rscript -e ‘devtools::document()’
- Run tests: Rscript -e ‘devtools::test()’
- Full check (must be 0 errors, 0 warnings): Rscript -e
  ‘devtools::check()’
- Report PDF: quarto render vignettes/articles/report.qmd –to pdf
- Website: Rscript -e ‘pkgdown::build_site()’
- Make small git commits with clear messages after each working step.

## Deadline

- Submit by 2026-11-23 (extended); oral walkthrough before finals week.
