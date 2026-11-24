# Zero-inflated negative binomial regression by EM with IRLS M-steps

Fits the ZINB model: with probability `pi` an observation is a
structural zero, and otherwise it follows an NB2 distribution with mean
`mu`, where `log(mu) = X beta + offset` and `logit(pi) = Z gamma`. The
observed-data log-likelihood is maximized by the EM algorithm. The
E-step computes the posterior probability that each zero is structural;
the M-step fits a logistic regression with that fractional response (by
IRLS) and a weighted NB regression (IRLS for beta, Newton-Raphson for
theta). No general-purpose optimizer is used. Standard errors come from
the analytic observed information matrix of `(beta, gamma, log(theta))`.
The results reproduce `pscl::zeroinfl(..., dist = "negbin")`.

## Usage

``` r
zinbreg(
  formula,
  data,
  offset,
  control = zinbreg_control(),
  weights,
  subset,
  na.action
)
```

## Arguments

- formula:

  A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model covariates
  before `|`, zero-model covariates after it. Without `|`, the same
  covariates are used in both parts.

- data:

  A data frame containing the variables in `formula`.

- offset:

  Optional offset for the count model, evaluated in `data`. Offsets can
  also be given with [`offset()`](https://rdrr.io/r/stats/offset.html)
  in either part of the formula.

- control:

  A list from
  [`zinbreg_control()`](https://asahicn00.github.io/nbzinbJia/reference/zinbreg_control.md).

- weights:

  Optional non-negative prior weights.

- subset:

  Optional expression selecting the rows of `data` to use.

- na.action:

  How to handle missing values.

## Value

An object of class `"zinbreg"`: a list containing `coefficients` (a list
with elements `count` and `zero`), `vcov` (covariance matrix of all
coefficients and `log(theta)`), `theta`, `SE.theta`, `SE.logtheta`,
`loglik`, `iter`, `converged`, `fitted.values` (`(1 - pi) mu`), `mu`,
`pi`, `call`, and the model information needed by
[`predict.zinbreg()`](https://asahicn00.github.io/nbzinbJia/reference/zinbreg-methods.md).

## See also

[zinbreg-methods](https://asahicn00.github.io/nbzinbJia/reference/zinbreg-methods.md),
[`nbreg()`](https://asahicn00.github.io/nbzinbJia/reference/nbreg.md)

## Examples

``` r
data("bioChemists", package = "pscl")
fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
               data = bioChemists)
summary(fit)
#> 
#> Call:
#> zinbreg(formula = art ~ fem + mar + kid5 + phd + ment | ment, 
#>     data = bioChemists)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.3041 -0.7685 -0.2632  0.4671  6.3765 
#> 
#> Count model coefficients (negbin with log link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.404005   0.141715   2.851  0.00436 ** 
#> femWomen    -0.211906   0.071922  -2.946  0.00322 ** 
#> marMarried   0.139462   0.081193   1.718  0.08586 .  
#> kid5        -0.167624   0.052457  -3.195  0.00140 ** 
#> phd          0.001963   0.035586   0.055  0.95601    
#> ment         0.024393   0.003518   6.934 4.08e-12 ***
#> Log(theta)   1.002892   0.142829   7.022 2.19e-12 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  -0.8067     0.3532  -2.284   0.0224 *
#> ment         -0.6095     0.2458  -2.480   0.0131 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Theta = 2.726  (Std. Err.: 0.3894 )
#> Log-likelihood: -1553.27 on 9 Df,  AIC: 3124.54 
#> EM iterations: 139 (converged) 
head(predict(fit, type = "zero"))
#>            1            2            3            4            5            6 
#> 6.223342e-03 1.138847e-02 1.138847e-02 6.690806e-02 5.852532e-08 1.165330e-01 
```
