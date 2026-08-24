# Fit the multi-outcome restraint-falls model

Compiles and fits the final Bayesian hierarchical model (M12) with NUTS
via cmdstanr: multinomial fall outcomes (total / minor / major) linked
to the observed restraint count, with per-pair quadratic restraint
trends and hierarchical fall rates and outcome shares.

## Usage

``` r
fit_model(
  data = restraints_falls,
  prepared = prepare_model_data(data),
  seed = 919214260,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  quick = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame with the structure of
  [restraints_falls](https://upipa.github.io/restraints-falls-nh/reference/restraints_falls.md).
  Ignored if `prepared` is supplied.

- prepared:

  Optional output of
  [`prepare_model_data()`](https://upipa.github.io/restraints-falls-nh/reference/prepare_model_data.md);
  computed from `data` if not supplied.

- seed:

  Random seed passed to the sampler. The default reproduces the
  published fit.

- chains:

  Number of MCMC chains.

- iter_warmup:

  Warmup iterations per chain.

- iter_sampling:

  Sampling iterations per chain.

- adapt_delta:

  Target acceptance probability. The default (0.95) reproduces the
  published fit exactly. The Stan default (0.8) is likely sufficient and
  considerably faster — worth trying if you do not need an exact
  reproduction.

- quick:

  If `TRUE`, runs a single short chain (200 warmup + 100 sampling
  iterations) to verify that the model compiles and samples, before
  committing to the full fit. Overrides `chains`, `iter_warmup` and
  `iter_sampling`.

- ...:

  Further arguments passed to the `sample` method of the compiled model.

## Value

A `CmdStanMCMC` object with the posterior draws.

## Details

Fitting the full model takes several hours (about 6 h for the published
fit) and requires a C++ toolchain (Rtools on Windows) to compile the
Stan model. cmdstanr is not on CRAN; install it with
`install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_model(quick = TRUE) # quick check (minutes)
fit <- fit_model() # full fit, reproduces the published analysis (~6 h)
} # }
```
