# Control parameters for zinbreg()

Control parameters for zinbreg()

## Usage

``` r
zinbreg_control(
  epsilon = 1e-10,
  maxit = 5000,
  inner_maxit = 100,
  trace = FALSE,
  theta_max = 1e+08,
  theta_warn = 1e+05
)
```

## Arguments

- epsilon:

  Convergence tolerance for the EM algorithm: iterations stop when the
  observed-data log-likelihood changes by less than `epsilon`. EM
  converges linearly (slowly), so a small change in the log-likelihood
  can still leave the parameters noticeably away from the MLE; the
  default `1e-10` gives coefficients accurate to about `1e-5`.

- maxit:

  Maximum number of EM iterations.

- inner_maxit:

  Maximum number of iterations of each IRLS / Newton-Raphson loop inside
  the M-steps.

- trace:

  Logical; if `TRUE`, print the observed-data log-likelihood at each EM
  iteration.

- theta_max:

  Upper bound for theta (see
  [`nbreg_control()`](https://asahicn00.github.io/nbzinbJia/reference/nbreg_control.md)).

- theta_warn:

  A warning is issued when the estimated theta exceeds this value.

## Value

A list of control parameters.

## Examples

``` r
zinbreg_control(epsilon = 1e-10, trace = TRUE)
#> $epsilon
#> [1] 1e-10
#> 
#> $maxit
#> [1] 5000
#> 
#> $inner_maxit
#> [1] 100
#> 
#> $trace
#> [1] TRUE
#> 
#> $theta_max
#> [1] 1e+08
#> 
#> $theta_warn
#> [1] 1e+05
#> 
```
