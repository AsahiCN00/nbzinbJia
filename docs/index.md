# nbzinbJia

**nbzinbJia** fits

- **negative binomial (NB2) regression** with
  [`nbreg()`](https://asahicn00.github.io/nbzinbJia/reference/nbreg.md),
  which reproduces
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html), and
- **zero-inflated negative binomial (ZINB) regression** with
  [`zinbreg()`](https://asahicn00.github.io/nbzinbJia/reference/zinbreg.md),
  which reproduces `pscl::zeroinfl(..., dist = "negbin")`.

Both models are estimated **without any general-purpose optimizer** (no
[`optim()`](https://rdrr.io/r/stats/optim.html),
[`nlm()`](https://rdrr.io/r/stats/nlm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html), …):

| Model | Algorithm |
|----|----|
| NB | IRLS (Fisher scoring) for β with θ fixed, alternated with Newton–Raphson for log θ with β fixed |
| ZINB | EM algorithm: the E-step computes the posterior probability τᵢ that each zero is structural; the M-steps are a logistic IRLS fit with fractional response τᵢ (zero part) and the NB fit above with prior weights 1 − τᵢ (count part) |

Standard errors for the ZINB model come from the **analytic** observed
information matrix of (β, γ, log θ).

This package was written for the STA 651 final project (Miami
University, Fall 2026).

**Function reference website** <https://asahicn00.github.io/nbzinbJia/>

The first implementation is preserved on the `first-pass` branch. The
commits on `main` reconstruct the development steps from saved snapshots
for the walkthrough.

## Installation

\
`# install.packages(c("devtools", "remotes"))`\
`devtools``::`[`install_github`](https://devtools.r-lib.org/reference/install-deprecated.html)`(``"AsahiCN00/nbzinbJia"``)`

With devtools 2.5.0 or later, `install_github()` needs the **remotes**
package (devtools no longer installs it automatically), so install it
first if R asks for it. `remotes::install_github("AsahiCN00/nbzinbJia")`
or `pak::pak("AsahiCN00/nbzinbJia")` work as well.

The package needs only `stats` and `Formula`. The examples and tests
also use `MASS`, `pscl`, `AER` and `numDeriv`:

\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"MASS"``, ``"pscl"``, ``"AER"``, ``"numDeriv"``)``)`

## Usage

\
[`library`](https://rdrr.io/r/base/library.html)`(``nbzinbJia``)`\
[`data`](https://rdrr.io/r/utils/data.html)`(``"bioChemists"``, package ``=`` ``"pscl"``)`\
\
`# Negative binomial regression`\
`nb`` ``<-`` `[`nbreg`](https://asahicn00.github.io/nbzinbJia/reference/nbreg.md)`(``art`` ``~`` ``fem`` ``+`` ``mar`` ``+`` ``kid5`` ``+`` ``phd`` ``+`` ``ment``, data ``=`` ``bioChemists``)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``nb``)`\
`#>             Estimate Std. Error z value Pr(>|z|)`\
`#> (Intercept)   0.2561     0.1373  1.8649   0.0622`\
`#> femWomen     -0.2164     0.0726 -2.9795   0.0029`\
`#> marMarried    0.1505     0.0821  1.8331   0.0668`\
`#> kid5         -0.1764     0.0528 -3.3404   0.0008`\
`#> phd           0.0153     0.0359  0.4257   0.6703`\
`#> ment          0.0291     0.0032  9.0477   0.0000`\
`#> Theta: 2.2644   Log-likelihood: -1560.958 on 7 Df`\
\
`# Zero-inflated NB: count covariates | zero-inflation covariates`\
`zinb`` ``<-`` `[`zinbreg`](https://asahicn00.github.io/nbzinbJia/reference/zinbreg.md)`(``art`` ``~`` ``fem`` ``+`` ``mar`` ``+`` ``kid5`` ``+`` ``phd`` ``+`` ``ment`` ``|`` ``ment``,`\
`                data ``=`` ``bioChemists``)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``zinb``)``$``coefficients``$``zero`\
`#>             Estimate Std. Error z value Pr(>|z|)`\
`#> (Intercept)  -0.8067     0.3532 -2.2840   0.0224`\
`#> ment         -0.6095     0.2458 -2.4797   0.0131`\
\
`# Model comparison and prediction`\
[`c`](https://rdrr.io/r/base/c.html)`(``NB ``=`` `[`AIC`](https://rdrr.io/r/stats/AIC.html)`(``nb``)``, ZINB ``=`` `[`AIC`](https://rdrr.io/r/stats/AIC.html)`(``zinb``)``)`\
`#>     NB   ZINB`\
`#> 3135.9 3124.5`\
\
`nd`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``fem ``=`` ``"Women"``, mar ``=`` ``"Married"``, kid5 ``=`` ``0``, phd ``=`` ``3``,`\
`                 ment ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``10``, ``30``)``)`\
[`predict`](https://rdrr.io/r/stats/predict.html)`(``zinb``, newdata ``=`` ``nd``, type ``=`` ``"zero"``)``   ``# P(structural zero)`\
`#>      1      2      3`\
`#> 0.3086 0.0010 0.0000`

Without `|`,
[`zinbreg()`](https://asahicn00.github.io/nbzinbJia/reference/zinbreg.md)
uses the same covariates in both parts, as `zeroinfl()` does. Offsets
can be given with the `offset` argument (count part) or with
[`offset()`](https://rdrr.io/r/stats/offset.html) in either part of the
formula.

### Main functions

| Function | Purpose |
|----|----|
| [`nbreg()`](https://asahicn00.github.io/nbzinbJia/reference/nbreg.md), [`nbreg_control()`](https://asahicn00.github.io/nbzinbJia/reference/nbreg_control.md) | NB regression and its tolerances / iteration limits / trace |
| [`zinbreg()`](https://asahicn00.github.io/nbzinbJia/reference/zinbreg.md), [`zinbreg_control()`](https://asahicn00.github.io/nbzinbJia/reference/zinbreg_control.md) | ZINB regression and its tolerances / iteration limits / trace |

Both fitted classes (`"nbreg"`, `"zinbreg"`) have `print`, `summary`,
`coef`, `vcov`, `logLik` (so [`AIC()`](https://rdrr.io/r/stats/AIC.html)
and [`BIC()`](https://rdrr.io/r/stats/AIC.html) work), `fitted`,
`residuals` (`"response"`, `"pearson"`), `nobs` and `predict` (with
`newdata`) methods.
[`predict.zinbreg()`](https://asahicn00.github.io/nbzinbJia/reference/zinbreg-methods.md)
supports `type = "response"` ((1 − π)μ), `"count"` (μ), `"zero"` (π) and
`"prob"` (P(Y = k)).

## Agreement with the reference implementations

Maximum absolute differences on
[`pscl::bioChemists`](https://rdrr.io/pkg/pscl/man/bioChemists.html) and
[`AER::NMES1988`](https://rdrr.io/pkg/AER/man/NMES1988.html) (checked in
`tests/testthat/`):

|                       | NB vs `glm.nb()` | ZINB vs `zeroinfl()` |
|-----------------------|------------------|----------------------|
| Coefficients          | ≤ 5e-8           | ≤ 2e-5               |
| Standard errors       | ≤ 4e-10          | ≤ 9e-6               |
| θ (NB) / log θ (ZINB) | ≤ 2e-8           | ≤ 2e-6               |
| Log-likelihood        | ≤ 7e-11          | ≤ 2e-9               |

## R CMD check

[`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
finishes with 0 errors, 0 warnings and 0 notes.
