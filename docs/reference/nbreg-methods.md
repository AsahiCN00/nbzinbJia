# Methods for nbreg objects

Standard S3 methods for fitted negative binomial regression models
returned by
[`nbreg()`](https://asahicn00.github.io/nbzinbJia/reference/nbreg.md).

## Usage

``` r
# S3 method for class 'nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
summary(object, ...)

# S3 method for class 'summary.nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
coef(object, ...)

# S3 method for class 'nbreg'
vcov(object, ...)

# S3 method for class 'nbreg'
logLik(object, ...)

# S3 method for class 'nbreg'
nobs(object, ...)

# S3 method for class 'nbreg'
fitted(object, ...)

# S3 method for class 'nbreg'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'nbreg'
predict(object, newdata, type = c("link", "response"), ...)
```

## Arguments

- x, object:

  A fitted `"nbreg"` object (or, for `print.summary.nbreg()`, a
  `"summary.nbreg"` object).

- digits:

  Number of significant digits to print.

- ...:

  Further arguments (ignored).

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"response"` (`y - mu`) or `"pearson"`
  (`(y - mu) / sqrt(mu + mu^2 / theta)`). For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"link"` (linear
  predictor) or `"response"` (`mu`).

- newdata:

  Optional data frame for prediction. If omitted, the fitted values from
  the original data are returned.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `"summary.nbreg"` with a coefficient table (estimates, standard
errors, z values, p-values).
[`logLik()`](https://rdrr.io/r/stats/logLik.html) returns a `"logLik"`
object with a `df` attribute (number of coefficients plus one for
theta), so [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) work. The other methods
return the corresponding vectors or matrices.

## Examples

``` r
data("bioChemists", package = "pscl")
fit <- nbreg(art ~ fem + ment, data = bioChemists)
coef(fit)
#> (Intercept)    femWomen        ment 
#>  0.30153194 -0.17439736  0.02929749 
sqrt(diag(vcov(fit)))
#> (Intercept)    femWomen        ment 
#> 0.057422210 0.069671041 0.003126916 
logLik(fit)
#> 'log Lik.' -1566.722 (df=4)
BIC(fit)
#> [1] 3160.72
head(residuals(fit, type = "pearson"))
#>          1          2          3          4          5          6 
#> -0.9708665 -0.9139980 -0.9139980 -0.9383033 -1.0724182 -0.8808327 
predict(fit, newdata = data.frame(fem = "Women", ment = c(0, 10)),
        type = "response")
#>        1        2 
#> 1.135570 1.522128 
```
