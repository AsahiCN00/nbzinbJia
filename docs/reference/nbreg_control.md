# Control parameters for nbreg()

Control parameters for nbreg()

## Usage

``` r
nbreg_control(
  epsilon = 1e-08,
  maxit = 100,
  trace = FALSE,
  theta_max = 1e+08,
  theta_warn = 1e+05
)
```

## Arguments

- epsilon:

  Convergence tolerance: the outer loop stops when the log-likelihood
  changes by less than `epsilon` between iterations.

- maxit:

  Maximum number of outer iterations (IRLS for beta alternated with
  Newton-Raphson for theta). Also caps the inner IRLS and Newton-Raphson
  loops.

- trace:

  Logical; if `TRUE`, print the log-likelihood and theta at each outer
  iteration.

- theta_max:

  Upper bound for theta. Estimates that reach it indicate the Poisson
  limit (no overdispersion).

- theta_warn:

  A warning is issued when the estimated theta exceeds this value.

## Value

A list of control parameters.

## Examples

``` r
nbreg_control(epsilon = 1e-10, trace = TRUE)
#> $epsilon
#> [1] 1e-10
#> 
#> $maxit
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
