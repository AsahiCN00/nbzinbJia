# Negative binomial regression by IRLS and Newton-Raphson

Fits the NB2 regression model `log(mu) = X beta + offset`,
`Var(Y) = mu + mu^2 / theta`. The regression coefficients are updated by
iteratively reweighted least squares (Fisher scoring) with theta held
fixed, and `log(theta)` is updated by Newton-Raphson with beta held
fixed. The two steps alternate until the log-likelihood converges. No
general-purpose optimizer is used. The results reproduce
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html).

## Usage

``` r
nbreg(
  formula,
  data,
  offset,
  control = nbreg_control(),
  weights,
  subset,
  na.action
)
```

## Arguments

- formula:

  A model formula, e.g. `y ~ x1 + x2`. Offsets may be given in the
  formula with [`offset()`](https://rdrr.io/r/stats/offset.html).

- data:

  A data frame containing the variables in `formula`.

- offset:

  Optional offset for the linear predictor (e.g. log exposure),
  evaluated in `data`.

- control:

  A list from
  [`nbreg_control()`](https://asahicn00.github.io/nbzinbJia/reference/nbreg_control.md).

- weights:

  Optional non-negative prior weights.

- subset:

  Optional expression selecting the rows of `data` to use.

- na.action:

  How to handle missing values (default from `getOption("na.action")`).

## Value

An object of class `"nbreg"`: a list containing `coefficients`, `vcov`
(variance-covariance matrix of the coefficients, treating theta as
known), `theta`, `SE.theta`, `loglik`, `iter`, `converged`,
`fitted.values`, `linear.predictors`, `call`, and the model information
needed by
[`predict.nbreg()`](https://asahicn00.github.io/nbzinbJia/reference/nbreg-methods.md).

## See also

[nbreg-methods](https://asahicn00.github.io/nbzinbJia/reference/nbreg-methods.md)
for summary(), predict() and other methods.

## Examples

``` r
data("bioChemists", package = "pscl")
fit <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
summary(fit)
#> 
#> Call:
#> nbreg(formula = art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.2677 -0.8754 -0.2613  0.4983  6.6573 
#> 
#> Coefficients (negative binomial with log link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.256144   0.137348   1.865 0.062191 .  
#> femWomen    -0.216418   0.072636  -2.979 0.002887 ** 
#> marMarried   0.150489   0.082097   1.833 0.066791 .  
#> kid5        -0.176415   0.052813  -3.340 0.000837 ***
#> phd          0.015271   0.035873   0.426 0.670326    
#> ment         0.029082   0.003214   9.048  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Theta: 2.264  Std. Err.: 0.2712 
#> Log-likelihood: -1560.96 on 7 Df,  AIC: 3135.92 
#> Outer iterations (IRLS / Newton-Raphson): 4 (converged) 
AIC(fit)
#> [1] 3135.917
```
