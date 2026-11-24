# Methods for zinbreg objects

Standard S3 methods for fitted zero-inflated negative binomial models
returned by
[`zinbreg()`](https://asahicn00.github.io/nbzinbJia/reference/zinbreg.md).
Coefficient names follow
[`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html):
`count_*` for the count model and `zero_*` for the zero-inflation model.

## Usage

``` r
# S3 method for class 'zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
summary(object, ...)

# S3 method for class 'summary.zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
coef(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
vcov(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
logLik(object, ...)

# S3 method for class 'zinbreg'
nobs(object, ...)

# S3 method for class 'zinbreg'
fitted(object, ...)

# S3 method for class 'zinbreg'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'zinbreg'
predict(
  object,
  newdata,
  type = c("response", "count", "zero", "prob"),
  at = NULL,
  ...
)
```

## Arguments

- x, object:

  A fitted `"zinbreg"` object (or, for `print.summary.zinbreg()`, a
  `"summary.zinbreg"` object).

- digits:

  Number of significant digits to print.

- ...:

  Further arguments (ignored).

- model:

  For [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html): `"full"` (both parts),
  `"count"`, or `"zero"`.

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"response"` (`y - (1 - pi) mu`) or `"pearson"` (response residual
  divided by `sqrt((1 - pi) mu (1 + pi mu + mu / theta))`). For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"response"`
  (`(1 - pi) mu`, the mean of `Y`), `"count"` (`mu`, the mean of the NB
  component), `"zero"` (`pi`, the probability of a structural zero), or
  `"prob"` (a matrix of `P(Y = k)` for `k` in `at`).

- newdata:

  Optional data frame for prediction. If omitted, the original data are
  used.

- at:

  For `type = "prob"`: the counts at which to evaluate probabilities
  (default `0:max(y)`).

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `"summary.zinbreg"` with coefficient tables (estimates, standard
errors, z values, p-values) for the count part (including `Log(theta)`)
and the zero part. [`logLik()`](https://rdrr.io/r/stats/logLik.html)
returns a `"logLik"` object with a `df` attribute (all coefficients plus
one for theta), so [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) work.

## Examples

``` r
data("bioChemists", package = "pscl")
fit <- zinbreg(art ~ fem + kid5 + ment | ment, data = bioChemists)
coef(fit, model = "zero")
#> (Intercept)        ment 
#>  -0.7949061  -0.6193810 
sqrt(diag(vcov(fit)))
#> count_(Intercept)    count_femWomen        count_kid5        count_ment 
#>       0.071295206       0.071503972       0.047968403       0.003456912 
#>  zero_(Intercept)         zero_ment 
#>       0.350163301       0.249061969 
AIC(fit)
#> [1] 3123.505
nd <- data.frame(fem = "Women", kid5 = 0, ment = c(0, 10, 30))
predict(fit, newdata = nd, type = "zero")
#>            1            2            3 
#> 3.111162e-01 9.213803e-04 3.845603e-09 
predict(fit, newdata = nd, type = "response")
#>         1         2         3 
#> 0.8988515 1.6630515 2.7091566 
# expected number of zeros
sum(predict(fit, type = "prob", at = 0))
#> [1] 286.3748
```
